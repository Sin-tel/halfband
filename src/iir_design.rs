use std::f64::consts::PI;

// Compute coefficients for given attenuation and transition bandwidth
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

// Compute coefficients given fixed N and transition bandwidth (TBW)
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
        q_ii1 *= ((((i * 2 + 1) * c) as f64 * PI) / (order as f64)).sin() * j as f64;
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
        q_i2 *= (((i * 2 * c) as f64 * PI) / (order as f64)).cos() * j as f64;
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
    // Expected output from HIIR
    use crate::iir_design::compute_coefs;
    use crate::iir_design::compute_coefs_tbw;

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
}
