//! Utilities for designing FIR windowed sinc filters.

use crate::fir::{Downsampler, Upsampler};
use std::f64::consts::PI;
use windowfunctions::{Symmetry, WindowFunction, window};

/// Computes FIR coefficients by estimating the best N and Beta for a Kaiser window.
///
/// * `attenuation_db`: Target stopband attenuation in decibels (e.g., 96.0).
/// * `transition`: Transition bandwidth relative to the high-rate Nyquist (0.0 to 0.5).
pub fn kaiser_spec(attenuation_db: f64, transition: f64) -> Vec<f32> {
    assert!(attenuation_db > 0.0);
    assert!(transition > 0.0);
    assert!(transition < 0.5);

    let n_coefs = estimate_n(attenuation_db, transition);
    let beta = estimate_beta(attenuation_db);
    kaiser_beta(n_coefs, beta)
}

/// Computes FIR coefficients for a fixed number of stages and stopband attenuation.
///
/// * `n_coefs`: Number of coefficients.
/// * `attenuation_db`: Target stopband attenuation in decibels (e.g., 96.0).
pub fn kaiser_attenuation(n_coefs: usize, attenuation_db: f64) -> Vec<f32> {
    assert!(attenuation_db > 0.0);

    let beta = estimate_beta(attenuation_db);
    kaiser_beta(n_coefs, beta)
}

/// Computes FIR coefficients for a Kaiser window with parameter beta.
pub fn kaiser_beta(n_coefs: usize, beta: f32) -> Vec<f32> {
    fir_coefs(n_coefs, WindowFunction::Kaiser { beta })
}

/// Computes FIR coefficients for a Hamming window.
pub fn hamming(n_coefs: usize) -> Vec<f32> {
    fir_coefs(n_coefs, WindowFunction::Hamming)
}

/// Estimates the required Kaiser beta parameter for a given stopband attenuation.
pub fn estimate_beta(attenuation_db: f64) -> f32 {
    (if attenuation_db > 50.0 {
        0.1102 * (attenuation_db - 8.7)
    } else if attenuation_db > 21.0 {
        0.5842 * (attenuation_db - 21.0).powf(0.4) + 0.07886 * (attenuation_db - 21.0)
    } else {
        0.0
    }) as f32
}

/// Estimates the number of non-zero coefficients (N) needed to meet
/// a specific attenuation and transition bandwidth.
pub fn estimate_n(attenuation_db: f64, transition: f64) -> usize {
    // k is the filter order
    let k = (attenuation_db - 8.0) / (14.36 * transition);

    // k = 4 * N - 1
    // N = (k + 1) / 4
    ((k + 1.0) / 4.0).round() as usize
}

fn fir_coefs(n_coefs: usize, window_type: WindowFunction) -> Vec<f32> {
    let ntaps = 4 * n_coefs - 1;
    let center = (ntaps / 2) as isize;

    // Generate the window
    let w: Vec<f64> = window::<f64>(ntaps, window_type, Symmetry::Symmetric).collect();

    // Calculate windowed sinc and extract non-zero coefficients
    // We want indices 0, 2, 4 ... just before the center.
    let mut coefs = Vec::with_capacity(n_coefs);
    for i in 0..n_coefs {
        let idx = (i * 2) as isize;
        let x = (idx - center) as f64;

        // sinc(x/2) = sin(pi * x / 2) / (pi * x)
        // Note: x will never be 0 here because we are skipping the center tap.
        let sinc = (PI * x * 0.5).sin() / (PI * x);

        // We multiply by 2.0 to maintain unity gain in the polyphase structure
        // (since the 0.5 center tap is handled separately in your process_sample).
        coefs.push(2.0 * sinc * w[idx as usize]);
    }

    let s: f64 = coefs.iter().sum();

    let normalized = coefs.iter().map(|x| (0.5 * *x / s) as f32).collect();

    normalized
}

/// The default designs use a Kaiser window with 53dB of attenuation.
impl<const N: usize> Default for Downsampler<N> {
    fn default() -> Self {
        let coefs = kaiser_attenuation(N, 53.0);
        Self::new(&coefs)
    }
}

/// The default designs use a Kaiser window with 53dB of attenuation.
impl<const N: usize> Default for Upsampler<N> {
    fn default() -> Self {
        let coefs = kaiser_attenuation(N, 53.0);
        Self::new(&coefs)
    }
}

#[cfg(test)]
mod tests {
    use crate::fir::design::*;

    #[test]
    fn test_sum() {
        let coefs = hamming(8);

        let s: f32 = coefs.into_iter().sum();

        assert_eq!(s, 0.5);
    }

    #[test]
    fn test_transition_n() {
        assert_eq!(estimate_n(53.0, 0.1), 8);
    }
}
