//! Utilities for designing IIR polyphase coefficients.
//!
//! These functions calculate coefficients for all-pass filters
//! based on desired attenuation and transition bandwidth.

use std::f64::consts::PI;
const TWO_PI: f64 = std::f64::consts::TAU;

/// Calculates the phase delay of the IIR Polyphase filter at a given frequency.
///
/// * `coef_arr`: The coefficients used in the filter.
/// * `f_fs`: Frequency relative to the high sample rate [0.0 to 0.5].
///   e.g., if you want the delay at 1kHz for a 88.2kHz high-rate signal,
///   f_fs = 1000. / 88200.
pub fn compute_phase_delay(coef_arr: &[f32], f_fs: f64) -> f64 {
    assert!(
        (0.0..0.5).contains(&f_fs),
        "Frequency must be in range [0, 0.5)"
    );

    // Starting delay for the two branches
    let mut delay = [0.0, 1.0];

    for (k, &c) in coef_arr.iter().enumerate() {
        let dly = compute_unit_phase_delay(f64::from(c), f_fs);

        // Accumulate into even/odd branch
        delay[k % 2] += dly;
    }

    let delay_even = delay[0];
    let delay_odd = delay[1];

    if f_fs > 0.0 {
        let w = TWO_PI * f_fs;
        // Phase difference between the two paths
        let phi = (delay_odd - delay_even) * w;

        // Sum the complex vectors of both paths:
        let re = phi.cos() + 1.0;
        let im = phi.sin();

        // The resulting angle is the relative delay contribution
        let dif = f64::atan2(im, re) / w;
        delay_even + dif
    } else {
        // At DC, the delay is the average of the two paths
        (delay_even + delay_odd) * 0.5
    }
}

// Computes the phase delay introduced by a single filter stage
fn compute_unit_phase_delay(a: f64, f_fs: f64) -> f64 {
    let w = TWO_PI * f_fs;

    if w > 0.0 {
        let ac = (a + 1.0) * w.cos();
        let as_val = (a - 1.0) * w.sin();

        let x = ac * ac - as_val * as_val;
        let y = 2.0 * ac * as_val;

        let mut ph = f64::atan2(-y, x);

        // Ensure phase is wrapped correctly
        if ph < 0.0 {
            ph += TWO_PI;
        }

        ph / w
    } else {
        // Limit as frequency approaches 0: -2 * (a - 1) / (a + 1)
        -2.0 * (a - 1.0) / (a + 1.0)
    }
}

/// Computes IIR coefficients for a given stopband attenuation and transition bandwidth.
///
/// * `attenuation`: Target stopband attenuation in decibels (e.g., 96.0).
/// * `transition`: Transition bandwidth relative to the high-rate Nyquist (0.0 to 0.5).
pub fn compute_coefs(attenuation: f64, transition: f64) -> Vec<f32> {
    assert!(attenuation > 0.0);
    assert!(transition > 0.0);
    assert!(transition < 0.5);

    let (k, q) = compute_transition_param(transition);

    let order = compute_order(attenuation, q);
    let n_coefs = (order - 1) / 2;

    let mut coefs = Vec::with_capacity(n_coefs);

    for i in 0..n_coefs {
        coefs.push(compute_coef(i, k, q, order));
    }

    coefs.iter().map(|x| *x as f32).collect()
}

/// Computes IIR coefficients for a fixed number of stages and transition bandwidth.
///
/// * `n_coefs`: Number of coefficients (total filter order is `n_coefs * 2 + 1`).
/// * `transition`: Transition bandwidth relative to the high-rate Nyquist (0.0 to 0.5).
pub fn compute_coefs_tbw(n_coefs: usize, transition: f64) -> Vec<f32> {
    assert!(n_coefs > 0);
    assert!(transition > 0.0);
    assert!(transition < 0.5);

    let (k, q) = compute_transition_param(transition);

    let mut coefs = Vec::with_capacity(n_coefs);
    let order = n_coefs * 2 + 1;
    for i in 0..n_coefs {
        coefs.push(compute_coef(i, k, q, order));
    }

    coefs.iter().map(|x| *x as f32).collect()
}

fn compute_order(attenuation: f64, q: f64) -> usize {
    assert!(attenuation > 0.0);
    assert!(q > 0.0);

    let attn_p2 = 10.0_f64.powf(-attenuation / 10.0);

    let a = attn_p2 / (1.0 - attn_p2);

    let mut order = ((a * a / 16.0).ln() / q.ln()).ceil() as usize;
    if order.is_multiple_of(2) {
        order += 1;
    }
    if order == 1 {
        order = 3;
    }

    order
}

fn compute_transition_param(transition: f64) -> (f64, f64) {
    assert!(transition > 0.);
    assert!(transition < 0.5);

    let mut k = ((1. - transition * 2.) * PI / 4.).tan();
    k *= k;
    assert!(k < 1.);
    assert!(k > 0.);
    let kksqrt: f64 = (1. - k * k).powf(0.25);
    let e = 0.5 * (1. - kksqrt) / (1. + kksqrt);
    let e4 = e.powi(4);
    let q = e * (1. + e4 * (2. + e4 * (15. + 150. * e4)));
    assert!(q > 0.);
    (k, q)
}

fn compute_coef(index: usize, k: f64, q: f64, order: usize) -> f64 {
    assert!(index * 2 < order);

    let c = index + 1;
    let num: f64 = compute_acc_num(q, order, c) * q.powf(0.25);
    let den: f64 = compute_acc_den(q, order, c) + 0.5;
    let ww = num / den;
    let wwsq = ww * ww;

    let x = ((1. - wwsq * k) * (1. - wwsq / k)).sqrt() / (1. + wwsq);
    assert!(!x.is_nan());

    (1. - x) / (1. + x)
}

fn compute_acc_num(q: f64, order: usize, c: usize) -> f64 {
    assert!(c >= 1);
    assert!(c < order * 2);

    let mut i = 0;
    let mut j = 1;
    let mut acc: f64 = 0.;
    let mut q_ii1;
    loop {
        q_ii1 = q.powi((i * (i + 1)).try_into().unwrap());
        q_ii1 *= ((((i * 2 + 1) * c) as f64 * PI) / (order as f64)).sin() * f64::from(j);
        acc += q_ii1;

        j = -j;
        i += 1;

        if q_ii1.abs() <= 1e-100 {
            break;
        }
    }

    assert!(!acc.is_nan());
    acc
}

fn compute_acc_den(q: f64, order: usize, c: usize) -> f64 {
    assert!(c >= 1);
    assert!(c < order * 2);
    let mut i = 1;
    let mut j = -1;
    let mut acc: f64 = 0.;
    let mut q_i2;
    loop {
        q_i2 = q.powi((i * i).try_into().unwrap());
        q_i2 *= (((i * 2 * c) as f64 * PI) / (order as f64)).cos() * f64::from(j);
        acc += q_i2;

        j = -j;
        i += 1;

        if q_i2.abs() <= 1e-100 {
            break;
        }
    }

    assert!(!acc.is_nan());
    acc
}

#[cfg(test)]
mod tests {
    use crate::iir::design::*;
    use approx::assert_relative_eq;

    // Expected output from HIIR
    const EXPECTED: [f32; 8] = [
        0.0771150813,
        0.265968531,
        0.482070625,
        0.665104151,
        0.796820462,
        0.88410151,
        0.941251457,
        0.982005417,
    ];

    #[test]
    fn test_compute_coefs_tbw() {
        let coefs = compute_coefs_tbw(8, 0.01);

        for (actual, expected) in coefs.iter().zip(EXPECTED.iter()) {
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_compute_coefs() {
        let coefs = compute_coefs(64.0, 0.01);

        for (actual, expected) in coefs.iter().zip(EXPECTED.iter()) {
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_phase_delay() {
        let coefs = compute_coefs_tbw(8, 0.0343747);
        let sample_rate = 44100.0;
        let freq_hz = 4000.0;

        let f_relative = freq_hz / (2.0 * sample_rate);
        let delay = compute_phase_delay(&coefs, f_relative.into()) - 0.5;

        assert_relative_eq!(delay, 3.0, epsilon = 1e-6);
    }
}
