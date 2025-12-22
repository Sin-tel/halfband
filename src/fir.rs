//! Finite Impulse Response (FIR) polyphase resamplers.
//!
//! These resamplers provide linear-phase filtering, meaning all frequencies
//! are delayed by the same amount.

pub mod presets;

use bit_mask_ring_buf::BitMaskRB;

/// A 2x FIR downsampler.
///
/// Implements a symmetric half-band FIR filter. `N` is the number of
/// non-zero coefficients in the polyphase branch.
/// Total filter length is `4 * N - 1`.
#[derive(Debug)]
pub struct Downsampler<const N: usize> {
    buf1: BitMaskRB<f32>,
    buf2: BitMaskRB<f32>,
    pos: isize,
    coef: [f32; N],
}
impl<const N: usize> Downsampler<N> {
    /// Creates a new downsampler from windowed-sinc coefficients.
    ///
    /// # Panics
    /// Panics if `coef_arr.len()` is not exactly `N`.
    pub fn new(coef_arr: &[f32]) -> Self {
        assert_eq!(coef_arr.len(), N, "Coefficient array size mismatch");

        let mut coef = [0.0; N];
        coef[..N].copy_from_slice(&coef_arr[..N]);

        Self {
            buf1: BitMaskRB::<f32>::new(2 * N, 0.0),
            buf2: BitMaskRB::<f32>::new(2 * N, 0.0),
            pos: 0,
            coef,
        }
    }

    /// Processes two high-rate samples (`s1`, `s2`) to produce one low-rate sample.
    pub fn process_sample(&mut self, s1: f32, s2: f32) -> f32 {
        self.pos = self.buf1.constrain(self.pos + 1);
        self.buf1[self.pos] = s1;
        self.buf2[self.pos] = s2;

        let k = N * 2 - 1;

        let s1 = self
            .coef
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                c * (self.buf1[self.pos - i as isize] + self.buf1[self.pos - (k - i) as isize])
            })
            .sum::<f32>();

        let s2 = self.buf2[self.pos - (N as isize)];
        0.5 * (s1 + s2)
    }

    /// Processes a block of high-rate samples into a low-rate output buffer.
    ///
    /// # Panics
    /// Panics if `input.len()` is not exactly `output.len() * 2`.
    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            output.len() * 2,
            input.len(),
            "Output must be twice the size of input."
        );

        for (i, chunk) in input.chunks_exact(2).enumerate() {
            output[i] = self.process_sample(chunk[0], chunk[1]);
        }
    }

    /// Resets internal filter state.
    pub fn clear(&mut self) {
        self.pos = 0;
        self.buf1.raw_data_mut().fill(0.);
        self.buf2.raw_data_mut().fill(0.);
    }

    /// Compute the latency of this stage at the half-band rate.
    pub fn get_latency(&self) -> f32 {
        let n_taps = 4 * N - 1;
        0.25 * ((n_taps as f32) - 1.0)
    }
}

/// A 2x FIR upsampler.
///
/// Implements a transposed symmetric half-band FIR filter.
#[derive(Debug)]
pub struct Upsampler<const N: usize> {
    buf: BitMaskRB<f32>,
    pos: isize,
    coef: [f32; N],
}

impl<const N: usize> Upsampler<N> {
    /// Creates a new upsampler from windowed-sinc coefficients.
    pub fn new(coef_arr: &[f32]) -> Self {
        assert_eq!(coef_arr.len(), N, "Coefficient array size mismatch");

        let mut coef = [0.0; N];
        coef[..N].copy_from_slice(&coef_arr[..N]);
        Self {
            buf: BitMaskRB::<f32>::new(2 * N, 0.0),
            pos: 0,
            coef,
        }
    }

    /// Processes one low-rate sample to produce two high-rate samples.
    ///
    /// Returns `[even_sample, odd_sample]`.
    #[rustfmt::skip]
    pub fn process_sample(&mut self, s: f32) -> [f32; 2] {
        self.pos = self.buf.constrain(self.pos + 1);
        self.buf[self.pos] = s;

        let k = N * 2 - 1;
        let s1 = self
            .coef
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                c * (self.buf[self.pos - i as isize] + self.buf[self.pos - (k - i) as isize])
            })
            .sum::<f32>();

        let s2 = self.buf[self.pos + 1 - N as isize];

        [s1, s2]
    }

    /// Processes a block of low-rate samples into a high-rate output buffer.
    ///
    /// # Panics
    /// Panics if `output.len()` is not exactly `input.len() * 2`.
    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            input.len() * 2,
            output.len(),
            "Input must be twice the size of output."
        );

        for (i, chunk) in output.chunks_exact_mut(2).enumerate() {
            let [even, odd] = self.process_sample(input[i]);
            chunk[0] = even;
            chunk[1] = odd;
        }
    }

    /// Resets internal filter state.
    pub fn clear(&mut self) {
        self.pos = 0;
        self.buf.raw_data_mut().fill(0.);
    }

    /// Compute the latency of this stage at the half-band rate.
    pub fn get_latency(&self) -> f32 {
        let n_taps = 4 * N - 1;
        0.25 * ((n_taps as f32) - 1.0)
    }
}

#[cfg(test)]
mod tests {
    use crate::fir::presets::*;

    // Test input (upsampled signal)
    const INPUT_DOWN: [f32; 32] = [
        1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
        1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
    ];

    // Expected output downsampled
    const EXPECTED_DOWN: [f32; 16] = [
        -0.00170315,
        0.0046452372,
        -0.011386225,
        0.025502931,
        -0.052331336,
        0.101509795,
        -0.19860508,
        0.51473564,
        0.30139494,
        -0.3984902,
        0.44766867,
        -0.47449708,
        0.48861378,
        -0.49535477,
        0.49829686,
        -0.5,
    ];

    const INPUT_UP: [f32; 16] = [
        1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
    ];

    const EXPECTED_UP: [f32; 32] = [
        -0.0034063,
        0.0,
        0.0024778747,
        0.0,
        -0.004191501,
        0.0,
        0.012273565,
        0.0,
        -0.02123189,
        0.0,
        0.032426544,
        0.0,
        -0.07460176,
        0.0,
        0.40564406,
        1.0,
        1.3391241,
        1.0,
        0.032426544,
        -1.0,
        -1.4349577,
        -1.0,
        0.012273565,
        1.0,
        1.4095343,
        1.0,
        0.0024778747,
        -1.0,
        -1.4171321,
        -1.0,
        0.0,
        1.0,
    ];

    #[test]
    fn test_downsampler() {
        let mut downsampler = Downsampler31::default();

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_DOWN.len()];
        downsampler.process_block(&INPUT_DOWN, &mut output);

        dbg!(output);

        // dbg!(output);
        for (actual, expected) in EXPECTED_DOWN.iter().zip(output.iter()) {
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_upsampler() {
        let mut upsampler = Upsampler31::default();

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_UP.len()];
        upsampler.process_block(&INPUT_UP, &mut output);

        dbg!(output);
        for (actual, expected) in EXPECTED_UP.iter().zip(output.iter()) {
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_clear_down() {
        let mut downsampler = Downsampler47::default();

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
    fn test_clear_up() {
        let mut upsampler = Upsampler19::default();

        // Process some data to populate internal state
        let test_input = [1.0, 1.0, -1.0, -1.0];
        let mut output1 = [0.0f32; 8];
        upsampler.process_block(&test_input, &mut output1);

        // Clear buffers and process same input again
        upsampler.clear();
        let mut output2 = [0.0f32; 8];
        upsampler.process_block(&test_input, &mut output2);

        // Results should be identical after clearing
        for i in 0..2 {
            assert_eq!(output1[i], output2[i]);
        }
    }
}
