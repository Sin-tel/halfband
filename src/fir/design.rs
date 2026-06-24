//! Utilities for designing FIR windowed sinc filters.
//!
//! These functions calculate coefficients for windowed sinc filters
//! based on desired attenuation and transition bandwidth.
//!
//! ## Example
//!
//! ```rust
//! use halfband::fir;
//! use halfband::fir::design::hamming;
//!
//! // Create a FIR upsampler using a 31-tap Hamming windowed sinc.
//! let coefs = hamming::<8>();
//! let up = fir::Upsampler8::new(&coefs);
//! ```

use crate::fir::{Downsampler, Upsampler};
use core::f64::consts::PI;
use heapless::Vec;
use windowfunctions::{Symmetry, WindowFunction, window};

/// Computes FIR coefficients for a fixed number of stages and stopband attenuation.
///
/// * `N`: Number of coefficients.
/// * `attenuation_db`: Target stopband attenuation in decibels (e.g., 96.0).
pub fn kaiser_attenuation<const N: usize>(attenuation_db: f64) -> Vec<f32, N> {
    assert!(attenuation_db > 0.0);

    let beta = estimate_beta(attenuation_db);
    kaiser_beta::<N>(beta)
}

/// Computes FIR coefficients for a Kaiser window with parameter beta.
pub fn kaiser_beta<const N: usize>(beta: f32) -> Vec<f32, N> {
    fir_coefs::<N>(WindowFunction::Kaiser { beta })
}

/// Computes FIR coefficients for a Hamming window.
pub fn hamming<const N: usize>() -> Vec<f32, N> {
    fir_coefs::<N>(WindowFunction::Hamming)
}

/// Estimates the required Kaiser beta parameter for a given stopband attenuation.
fn estimate_beta(attenuation_db: f64) -> f32 {
    (if attenuation_db > 50.0 {
        0.1102 * (attenuation_db - 8.7)
    } else if attenuation_db > 21.0 {
        0.5842 * (attenuation_db - 21.0).powf(0.4) + 0.07886 * (attenuation_db - 21.0)
    } else {
        0.0
    }) as f32
}

fn fir_coefs<const N: usize>(window_type: WindowFunction) -> Vec<f32, N> {
    let ntaps = 4 * N - 1;
    let center = (ntaps / 2) as isize;

    // worst case 64x oversampling, so taps: 4 * 64 - 1 = 255
    assert!(N <= 64);
    let w: Vec<f64, 255> = window::<f64>(ntaps, window_type, Symmetry::Symmetric).collect();

    // Calculate windowed sinc and extract non-zero coefficients
    // We want indices 0, 2, 4 ... just before the center.
    let mut coefs: Vec<f64, N> = Vec::new();
    for i in 0..N {
        let idx = (i * 2) as isize;
        let x = (idx - center) as f64;

        // x will never be 0 here because we are skipping the center tap.
        let sinc = (PI * x * 0.5).sin() / (PI * x);

        coefs.push(sinc * w[idx as usize]).ok();
    }

    let s: f64 = coefs.iter().sum();

    let normalized = coefs.iter().map(|x| (0.5 * *x / s) as f32).collect();

    normalized
}

/// The default designs use a Kaiser window with 53dB of attenuation.
impl<const N: usize> Default for Downsampler<N> {
    fn default() -> Self {
        let coefs = kaiser_attenuation::<N>(53.0);
        Self::new(&coefs)
    }
}

/// The default designs use a Kaiser window with 53dB of attenuation.
impl<const N: usize> Default for Upsampler<N> {
    fn default() -> Self {
        let coefs = kaiser_attenuation::<N>(53.0);
        Self::new(&coefs)
    }
}

#[cfg(test)]
mod tests {
    use crate::fir::design::*;

    #[test]
    fn test_sum() {
        let coefs = hamming::<8>();

        let s: f32 = coefs.into_iter().sum();

        assert_eq!(s, 0.5);
    }

    #[test]
    fn test_transition_n() {
        assert_eq!(estimate_n(53.0, 0.1), 8);
    }

    /// Estimates the number of non-zero coefficients (N) needed to meet
    /// a specific attenuation and transition bandwidth.
    const fn estimate_n(attenuation_db: f64, transition: f64) -> usize {
        // k is the filter order
        let k = (attenuation_db - 8.0) / (14.36 * transition);

        // k = 4 * N - 1
        // N = (k + 1) / 4
        ((k + 1.0) / 4.0) as usize
    }
}
