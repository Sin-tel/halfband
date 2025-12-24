//! Infinite Impulse Response (IIR) polyphase resamplers.
//!
//! These resamplers utilize all-pass filter chains to achieve very steep
//! attenuation with minimal CPU overhead. They are based on the designs
//! found in the HIIR C++ library.
//!
//! These filters internally use two parallel all-pass filter chains.
//! Each stage `N` consists of a first-order all-pass section.
//! For even filters, the number of coefficients is `2*N`.
//! For odd filters, the number of coefficients is `2*N-1`.
//!
//! Because of the generic parameters, it is often more convenient
//! to work with one of the aliases provided in [types].
//!
//! `Default` implementations are provided for the following types, which give an integer delay when used together as a pair:
//! |  Upsampler    | Downsampler     | total delay | attenuation |
//! |---------------|-----------------|-------------|-------------|
//! | `Upsampler12` | `Downsampler12` | 5.0 samples | 157 dB      |
//! | `Upsampler10` | `Downsampler10` | 4.0 samples | 125 dB      |
//! | `Upsampler8`  | `Downsampler8`  | 3.0 samples | 94 dB       |
//! | `Upsampler6`  | `Downsampler6`  | 2.0 samples | 64 dB       |

// TODO: Implement IIR Hilbert filters
// TODO: Add presets for 4x / 8x / 16x

pub mod design;
pub mod types;
pub use types::*;

use design::coefs_transition;
use design::phase_delay;

// Frequency to measure the phase delay at, since it varies.
// Equivalent to 4kHz when sample rate is 44.1kHz.
const MEASURE_F: f32 = 4000. / 44100.;

/// Core Polyphase trait implemented by [PolyphaseEven] and [PolyphaseOdd].
pub trait Polyphase {
    /// Creates a new polyphase filter from a slice of coefficients.
    fn new(coef_arr: &[f32]) -> Self;
    /// Processes a pair of samples through the parallel all-pass chains.
    fn process(&mut self, s0: f32, s1: f32) -> [f32; 2];
    /// Resets internal filter state.
    fn clear(&mut self);
    /// Get a reference to the coefficient vector.
    fn coefs(&self) -> &[f32];
}

/// A polyphase filter with an even number of coefficients.
#[derive(Debug)]
pub struct PolyphaseEven<const N: usize> {
    coef: [[f32; 2]; N],
    state: [[f32; 2]; N],
    state_last: [f32; 2],
    coef_vec: Vec<f32>,
}

impl<const N: usize> Polyphase for PolyphaseEven<N> {
    /// Creates a new polyphase filter from a slice of coefficients.
    fn new(coef_arr: &[f32]) -> Self {
        assert!(N > 0);
        let mut coef = [[0.0; 2]; N];

        assert_eq!(coef_arr.len(), N * 2, "Coefficient array size mismatch");
        for i in 0..N {
            coef[i] = [coef_arr[i * 2], coef_arr[i * 2 + 1]];
        }

        let coef_vec = coef_arr.to_vec();

        Self {
            coef,
            state: [[0.0; 2]; N],
            state_last: [0.0; 2],
            coef_vec,
        }
    }

    /// Processes a pair of samples through the parallel all-pass chains.
    #[inline]
    fn process(&mut self, s0: f32, s1: f32) -> [f32; 2] {
        // Note: swapped!
        let mut signal = [s1, s0];

        // First N-1 stages
        for i in 0..N - 1 {
            let tmp = [
                (signal[0] - self.state[i + 1][0]) * self.coef[i][0] + self.state[i][0],
                (signal[1] - self.state[i + 1][1]) * self.coef[i][1] + self.state[i][1],
            ];
            self.state[i] = signal;
            signal = tmp;
        }

        // Process last stage
        let i = N - 1;
        let tmp = [
            (signal[0] - self.state_last[0]) * self.coef[i][0] + self.state[i][0],
            (signal[1] - self.state_last[1]) * self.coef[i][1] + self.state[i][1],
        ];
        self.state[i] = signal;
        signal = tmp;

        self.state_last = signal;

        signal
    }

    /// Resets internal filter state.
    fn clear(&mut self) {
        self.state = [[0.0; 2]; N];
        self.state_last = [0.0; 2];
    }

    fn coefs(&self) -> &[f32] {
        &self.coef_vec
    }
}

/// A polyphase filter with an odd number of coefficients.
#[derive(Debug)]
pub struct PolyphaseOdd<const N: usize> {
    coef: [[f32; 2]; N],
    state: [[f32; 2]; N],
    state_last: f32,
    coef_vec: Vec<f32>,
}

impl<const N: usize> Polyphase for PolyphaseOdd<N> {
    /// Creates a new polyphase filter from a slice of coefficients.
    fn new(coef_arr: &[f32]) -> Self {
        assert!(N > 0);
        let mut coef = [[0.0; 2]; N];

        assert_eq!(coef_arr.len(), N * 2 - 1, "Coefficient array size mismatch");
        for i in 0..N - 1 {
            coef[i] = [coef_arr[i * 2], coef_arr[i * 2 + 1]];
        }
        coef[N - 1][0] = coef_arr[N * 2 - 2];

        let coef_vec = coef_arr.to_vec();

        Self {
            coef,
            state: [[0.0; 2]; N],
            state_last: 0.0,
            coef_vec,
        }
    }

    /// Processes a pair of samples through the parallel all-pass chains.
    #[inline]
    fn process(&mut self, s0: f32, s1: f32) -> [f32; 2] {
        // Note: swapped!
        let mut signal = [s1, s0];

        // First N-1 stages
        for i in 0..N - 1 {
            let tmp = [
                (signal[0] - self.state[i + 1][0]) * self.coef[i][0] + self.state[i][0],
                (signal[1] - self.state[i + 1][1]) * self.coef[i][1] + self.state[i][1],
            ];
            self.state[i] = signal;
            signal = tmp;
        }

        // Process last stage
        let i = N - 1;
        let tmp = (signal[0] - self.state_last) * self.coef[i][0] + self.state[i][0];
        self.state[i] = signal;
        signal[0] = tmp;

        self.state_last = signal[0];

        signal
    }

    /// Resets internal filter state.
    fn clear(&mut self) {
        self.state = [[0.0; 2]; N];
        self.state_last = 0.0;
    }

    fn coefs(&self) -> &[f32] {
        &self.coef_vec
    }
}

/// A 2x IIR downsampler.
///
/// Takes high-rate input samples and produces low-rate output.
#[derive(Debug)]
pub struct Downsampler<P: Polyphase> {
    filter: P,
}

impl<P: Polyphase> Downsampler<P> {
    /// Creates a new downsampler with the provided coefficients.
    ///
    /// # Panics
    /// Panics if `coef_arr.len()` does not match the number of coefficients specified.
    pub fn new(coef_arr: &[f32]) -> Self {
        Self {
            filter: P::new(coef_arr),
        }
    }

    /// Processes a pair of high-rate samples to produce one downsampled output.
    #[inline]
    pub fn process(&mut self, mut s0: f32, mut s1: f32) -> f32 {
        [s0, s1] = self.filter.process(s0, s1);
        0.5 * (s0 + s1)
    }

    /// Processes a block of high-rate samples into a low-rate output buffer.
    ///
    /// # Panics
    /// Panics if `input.len()` is not exactly `output.len() * 2`.
    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            output.len() * 2,
            input.len(),
            "Input must be twice the size of output."
        );

        for (i, chunk) in input.chunks_exact(2).enumerate() {
            output[i] = self.process(chunk[0], chunk[1]);
        }
    }

    /// Resets internal filter state.
    pub fn clear(&mut self) {
        self.filter.clear();
    }

    /// Compute the latency of this stage, at the original rate.
    /// `rate` is the ratio of high/low sample rates.
    /// E.g. for a 4x cascade, the outer pair runs at rate 2, the inner pair at 4.
    ///
    /// For IIR filters, the latency is frequency-dependent.
    /// This method always measures at a relative frequency equivalent to 4kHz when the sample rate is 44.1kHz.
    pub fn latency(&self, rate: u32) -> f32 {
        let rate = rate as f32;
        let f_relative = MEASURE_F / rate;

        let delay = phase_delay(self.filter.coefs(), f_relative.into());
        // Note extra -1 term here. Downsampler is always one step ahead since it receives a pair of samples,
        // the last one being one step "in the future".
        (delay - 1.0) as f32 / rate
    }
}

/// A 2x IIR upsampler.
///
/// Takes low-rate input samples and produces high-rate output samples.
#[derive(Debug)]
pub struct Upsampler<P: Polyphase> {
    filter: P,
}

impl<P: Polyphase> Upsampler<P> {
    /// Creates a new upsampler with the provided coefficients.
    ///
    /// # Panics
    /// Panics if `coef_arr.len()` does not match the number of coefficients specified.
    pub fn new(coef_arr: &[f32]) -> Self {
        Self {
            filter: P::new(coef_arr),
        }
    }

    /// Processes one low-rate sample to produce a pair of high-rate samples.
    ///
    /// Returns `[even_sample, odd_sample]`.
    #[inline]
    pub fn process(&mut self, input: f32) -> [f32; 2] {
        self.filter.process(input, input)
    }

    /// Processes a block of low-rate samples into a high-rate output buffer.
    ///
    /// # Panics
    /// Panics if `output.len()` is not exactly `input.len() * 2`.
    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            input.len() * 2,
            output.len(),
            "Output must be twice the size of input."
        );

        for (i, chunk) in output.chunks_exact_mut(2).enumerate() {
            let [even, odd] = self.process(input[i]);
            chunk[0] = even;
            chunk[1] = odd;
        }
    }

    /// Resets internal filter states.
    pub fn clear(&mut self) {
        self.filter.clear();
    }

    /// Compute the latency of this stage, at the original rate.
    /// `rate` is the ratio of high/low sample rates.
    /// E.g. for a 4x cascade, the outer pair runs at rate 2, the inner pair at 4.
    ///
    /// For IIR filters, the latency is frequency-dependent.
    /// This method always measures at a relative frequency equivalent to 4kHz when the sample rate is 44.1kHz.
    pub fn latency(&self, rate: u32) -> f32 {
        let rate = rate as f32;
        let f_relative = MEASURE_F / rate;

        let delay = phase_delay(self.filter.coefs(), f_relative.into());
        delay as f32 / rate
    }
}

#[cfg(test)]
mod tests {
    use crate::iir::design::*;
    use crate::iir::*;

    // Test coefficients from HIIR
    const TEST_COEFFICIENTS: [f32; 8] = [
        0.0771150813,
        0.265968531,
        0.482070625,
        0.665104151,
        0.796820462,
        0.88410151,
        0.941251457,
        0.982005417,
    ];

    // Test input
    const INPUT_DOWN: [f32; 16] = [
        1.0, 0.5, 0.0, -0.5, -1.0, -0.5, 0.0, 0.5, 1.0, 0.5, 0.0, -0.5, -1.0, -0.5, 0.0, 0.5,
    ];

    // Expected output from HIIR (downsampled)
    const EXPECTED_DOWN: [f32; 8] = [
        0.0837605968,
        0.452451766,
        0.183772057,
        -0.876450181,
        -0.0602405369,
        0.923236489,
        -0.0583606362,
        -0.784886956,
    ];

    // Test input
    const INPUT_UP: [f32; 8] = [1.0, 0.5, 0.0, -0.5, -1.0, -0.5, 0.0, 0.5];

    // Expected output from HIIR (upsampled)
    const EXPECTED_UP: [f32; 16] = [
        0.0278814975,
        0.153580442,
        0.433896393,
        0.785656631,
        0.955273688,
        0.726828396,
        0.191594124,
        -0.282384992,
        -0.452594221,
        -0.474814117,
        -0.649253607,
        -0.940763056,
        -0.981431603,
        -0.587845564,
        -0.0524485111,
        0.234299481,
    ];

    #[test]
    fn test_downsampler() {
        let mut downsampler = Downsampler8::new(&TEST_COEFFICIENTS);

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_DOWN.len()];
        downsampler.process_block(&INPUT_DOWN, &mut output);

        for (actual, expected) in EXPECTED_DOWN.iter().zip(output.iter()) {
            assert_eq!(actual, expected);
        }

        assert_eq!(downsampler.filter.coefs(), &TEST_COEFFICIENTS);
    }

    #[test]
    fn test_upsampler() {
        let mut upsampler = Upsampler8::new(&TEST_COEFFICIENTS);

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_UP.len()];
        upsampler.process_block(&INPUT_UP, &mut output);

        for (actual, expected) in output.iter().zip(EXPECTED_UP.iter()) {
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_single_sample_processing() {
        // Test that block processing and single sample processing give same results
        let mut downsampler_block = Downsampler8::new(&TEST_COEFFICIENTS);
        let mut downsampler_single = Downsampler8::new(&TEST_COEFFICIENTS);

        // Process as block
        let mut block_output = [0.0f32; EXPECTED_DOWN.len()];
        downsampler_block.process_block(&INPUT_DOWN, &mut block_output);

        // Process sample by sample
        let mut single_output = [0.0f32; EXPECTED_DOWN.len()];

        for i in 0..EXPECTED_DOWN.len() {
            single_output[i] = downsampler_single.process(INPUT_DOWN[i * 2], INPUT_DOWN[i * 2 + 1]);
        }

        for (block, single) in block_output.iter().zip(single_output.iter()) {
            assert_eq!(block, single);
        }
    }

    #[test]
    fn test_buffer_clearing() {
        let mut downsampler = Downsampler8::new(&TEST_COEFFICIENTS);

        // Process some data to populate internal state
        let test_input = [1.0, 1.0, -1.0, -1.0];
        let mut output1 = [0.0f32; 2];
        downsampler.process_block(&test_input, &mut output1);

        // Clear buffers and process same input again
        downsampler.clear();
        let mut output2 = [0.0f32; 2];
        downsampler.process_block(&test_input, &mut output2);

        // Results should be identical after clearing
        for i in 0..2 {
            assert_eq!(output1[i], output2[i]);
        }
    }

    #[test]
    fn test_phase_delay_cascade() {
        // From resampling.txt 4×, 120 dB
        let coefs1 = coefs_transition(10, 0.0367598);
        let coefs2 = coefs_transition(4, 0.261666);
        let downsampler1 = Downsampler10::new(&coefs1);
        let upsampler1 = Upsampler10::new(&coefs1);
        let downsampler2 = Downsampler4::new(&coefs2);
        let upsampler2 = Upsampler4::new(&coefs2);

        let delay1 = upsampler1.latency(2) + downsampler1.latency(2);
        let delay2 = upsampler2.latency(4) + downsampler2.latency(4);

        assert_eq!(delay1 + delay2, 5.);
    }

    // Tests for odd number of coefficients

    // Test coefficients from HIIR
    const TEST_COEFFICIENTS_ODD: [f32; 7] = [
        0.0976001173,
        0.324210763,
        0.561098695,
        0.741631627,
        0.85908848,
        0.931541979,
        0.979553342,
    ];

    // Expected output from HIIR (downsampled)
    const EXPECTED_DOWN_ODD: [f32; 8] = [
        0.12351343,
        0.516426027,
        -0.017348364,
        -0.923705459,
        0.252689898,
        0.831034064,
        -0.247902423,
        -0.78623116,
    ];

    // Expected output from HIIR (upsampled)
    const EXPECTED_UP_ODD: [f32; 16] = [
        0.0460845754,
        0.223984569,
        0.562953711,
        0.89793092,
        0.931892037,
        0.544287086,
        -0.015807271,
        -0.358057559,
        -0.42016384,
        -0.497870803,
        -0.785470068,
        -1.02987862,
        -0.869121671,
        -0.362747431,
        0.0737292171,
        0.219533384,
    ];

    #[test]
    fn test_downsampler_odd() {
        let mut downsampler = Downsampler7::new(&TEST_COEFFICIENTS_ODD);

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_DOWN_ODD.len()];
        downsampler.process_block(&INPUT_DOWN, &mut output);

        dbg!(output);

        for (actual, expected) in EXPECTED_DOWN_ODD.iter().zip(output.iter()) {
            assert_eq!(actual, expected);
        }

        assert_eq!(downsampler.filter.coefs(), &TEST_COEFFICIENTS_ODD);
    }

    #[test]
    fn test_upsampler_odd() {
        let mut upsampler = Upsampler7::new(&TEST_COEFFICIENTS_ODD);

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_UP_ODD.len()];
        upsampler.process_block(&INPUT_UP, &mut output);

        dbg!(output);

        for (actual, expected) in output.iter().zip(EXPECTED_UP_ODD.iter()) {
            assert_eq!(actual, expected);
        }
    }
}
