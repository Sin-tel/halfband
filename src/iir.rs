pub mod design;

pub struct Polyphase<const N: usize> {
    coef: [[f32; 2]; N],
    state: [[f32; 2]; N],
    state_last: [f32; 2],
}

impl<const N: usize> Polyphase<N> {
    pub fn new(coef_arr: &[f32]) -> Self {
        assert_eq!(coef_arr.len(), N * 2, "Coefficient array size mismatch");
        assert!(N > 0);

        let mut coef = [[0.0; 2]; N];
        for i in 0..N {
            coef[i] = [coef_arr[i * 2], coef_arr[i * 2 + 1]];
        }

        Self {
            coef,
            state: [[0.0; 2]; N],
            state_last: [0.0; 2],
        }
    }

    #[inline]
    pub fn process(&mut self, s0: f32, s1: f32) -> [f32; 2] {
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

    pub fn clear(&mut self) {
        self.state = [[0.0; 2]; N];
        self.state_last = [0.0; 2];
    }
}

pub struct Downsampler<const N: usize> {
    filter: Polyphase<N>,
}

impl<const N: usize> Downsampler<N> {
    pub fn new(coef_arr: &[f32]) -> Self {
        Self {
            filter: Polyphase::new(coef_arr),
        }
    }

    #[inline]
    pub fn process_sample(&mut self, mut s0: f32, mut s1: f32) -> f32 {
        [s0, s1] = self.filter.process(s0, s1);
        0.5 * (s0 + s1)
    }

    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            output.len() * 2,
            input.len(),
            "Input must be twice the size of output."
        );

        for (i, chunk) in input.chunks_exact(2).enumerate() {
            output[i] = self.process_sample(chunk[0], chunk[1]);
        }
    }

    pub fn clear(&mut self) {
        self.filter.clear();
    }
}

pub struct Upsampler<const N: usize> {
    filter: Polyphase<N>,
}

impl<const N: usize> Upsampler<N> {
    pub fn new(coef_arr: &[f32]) -> Self {
        Self {
            filter: Polyphase::new(coef_arr),
        }
    }

    #[inline]
    pub fn process_sample(&mut self, input: f32) -> [f32; 2] {
        self.filter.process(input, input)
    }

    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            input.len() * 2,
            output.len(),
            "Output must be twice the size of input."
        );

        for (i, chunk) in output.chunks_exact_mut(2).enumerate() {
            let [even, odd] = self.process_sample(input[i]);
            chunk[0] = even;
            chunk[1] = odd;
        }
    }

    pub fn clear(&mut self) {
        self.filter.clear();
    }
}

#[cfg(test)]
mod tests {
    use crate::iir::Downsampler;
    use crate::iir::Upsampler;

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
        let mut downsampler = Downsampler::<4>::new(&TEST_COEFFICIENTS);

        // Process the input in pairs
        let mut output = [0.0f32; EXPECTED_DOWN.len()];
        downsampler.process_block(&INPUT_DOWN, &mut output);

        for (actual, expected) in EXPECTED_DOWN.iter().zip(output.iter()) {
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn test_upsampler() {
        let mut upsampler = Upsampler::<4>::new(&TEST_COEFFICIENTS);

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
        let mut downsampler_block = Downsampler::<4>::new(&TEST_COEFFICIENTS);
        let mut downsampler_single = Downsampler::<4>::new(&TEST_COEFFICIENTS);

        // Process as block
        let mut block_output = [0.0f32; EXPECTED_DOWN.len()];
        downsampler_block.process_block(&INPUT_DOWN, &mut block_output);

        // Process sample by sample
        let mut single_output = [0.0f32; EXPECTED_DOWN.len()];

        for i in 0..EXPECTED_DOWN.len() {
            single_output[i] =
                downsampler_single.process_sample(INPUT_DOWN[i * 2], INPUT_DOWN[i * 2 + 1]);
        }

        for (block, single) in block_output.iter().zip(single_output.iter()) {
            assert_eq!(block, single);
        }
    }

    #[test]
    fn test_buffer_clearing() {
        let mut downsampler = Downsampler::<4>::new(&TEST_COEFFICIENTS);

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
}
