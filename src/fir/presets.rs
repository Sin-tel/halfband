//! Type aliases and `Default` implementations for standard FIR filter sizes.
//!
//! These presets use windowed sinc coefficients, using a Hamming window.
//! Tap counts refer to the total filter length (4N - 1).

use crate::fir::{Downsampler, Upsampler};

// Coefficients can be generated using `calc_fir.py`.
const COEF_19: [f32; 5] = [
    0.005627999,
    -0.016970186,
    0.05826522,
    -0.16250847,
    0.61558545,
];

const COEF_31: [f32; 8] = [
    -0.0034063,
    0.0058841747,
    -0.0134819755,
    0.028233415,
    -0.053656805,
    0.09835691,
    -0.19419056,
    0.63226116,
];

const COEF_47: [f32; 12] = [
    -0.0022190535,
    0.0029486136,
    -0.0049348283,
    0.008482312,
    -0.013966189,
    0.02190768,
    -0.033139557,
    0.049198344,
    -0.07339233,
    0.11443117,
    -0.20456062,
    0.6352444,
];

const COEF_63: [f32; 16] = [
    -0.0016455112,
    0.0019660334,
    -0.002769686,
    0.004143326,
    -0.006184378,
    0.009008168,
    -0.012760753,
    0.017641623,
    -0.023945283,
    0.03214175,
    -0.04304512,
    0.058206834,
    -0.08098234,
    0.12015598,
    -0.20806137,
    0.6361307,
];

/// A 2x Upsampler using a 19-tap FIR filter.
pub type Upsampler19 = Upsampler<{ COEF_19.len() }>;

impl Default for Upsampler19 {
    fn default() -> Self {
        Self::new(&COEF_19)
    }
}

/// A 2x Downsampler using a 19-tap FIR filter.
pub type Downsampler19 = Downsampler<{ COEF_19.len() }>;

impl Default for Downsampler19 {
    fn default() -> Self {
        Self::new(&COEF_19)
    }
}

/// A 2x Upsampler using a 31-tap FIR filter.
pub type Upsampler31 = Upsampler<{ COEF_31.len() }>;

impl Default for Upsampler31 {
    fn default() -> Self {
        Self::new(&COEF_31)
    }
}

/// A 2x Downsampler using a 31-tap FIR filter.
pub type Downsampler31 = Downsampler<{ COEF_31.len() }>;

impl Default for Downsampler31 {
    fn default() -> Self {
        Self::new(&COEF_31)
    }
}

/// A 2x Upsampler using a 47-tap FIR filter.
pub type Upsampler47 = Upsampler<{ COEF_47.len() }>;

impl Default for Upsampler47 {
    fn default() -> Self {
        Self::new(&COEF_47)
    }
}

/// A 2x Downsampler using a 47-tap FIR filter.
pub type Downsampler47 = Downsampler<{ COEF_47.len() }>;

impl Default for Downsampler47 {
    fn default() -> Self {
        Self::new(&COEF_47)
    }
}

/// A 2x Upsampler using a 63-tap FIR filter.
pub type Upsampler63 = Upsampler<{ COEF_63.len() }>;

impl Default for Upsampler63 {
    fn default() -> Self {
        Self::new(&COEF_63)
    }
}

/// A 2x Downsampler using a 63-tap FIR filter.
pub type Downsampler63 = Downsampler<{ COEF_63.len() }>;

impl Default for Downsampler63 {
    fn default() -> Self {
        Self::new(&COEF_63)
    }
}
