use bit_mask_ring_buf::BitMaskRB;

pub const COEF_19: [f32; 5] = [
    1.6543841e-04,
    -5.942591e-3,
    3.6400687e-02,
    -1.3846028e-01,
    6.078056e-1,
];

pub const COEF_31: [f32; 8] = [
    -9.926305e-5,
    1.2845506e-03,
    -5.469064e-3,
    1.604118e-2,
    -3.845655e-2,
    8.3076164e-02,
    -1.8245567e-01,
    6.261118e-1,
];

pub const COEF_51: [f32; 13] = [
    1.3438646e-06,
    -3.0550054e-05,
    1.7178609e-04,
    -6.266282e-4,
    1.7876519e-03,
    -4.317118e-3,
    9.225773e-3,
    -1.7981611e-02,
    3.2782678e-02,
    -5.7463527e-02,
    1.0096023e-01,
    -1.9530381e-01,
    6.307933e-1,
];

#[derive(Debug)]
pub struct Downsampler<const N: usize> {
    buf1: BitMaskRB<f32>,
    buf2: BitMaskRB<f32>,
    pos: isize,
    coef: [f32; N],
}
impl<const N: usize> Downsampler<N> {
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

    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            output.len() * 2,
            input.len(),
            "Output must be twice the size of input."
        );

        assert!(input.len().is_multiple_of(2), "Input length must be even.");

        for (i, chunk) in input.chunks_exact(2).enumerate() {
            output[i] = self.process_sample(chunk[0], chunk[1]);
        }
    }

    pub fn clear(&mut self) {
        self.pos = 0;
        self.buf1.raw_data_mut().fill(0.);
        self.buf2.raw_data_mut().fill(0.);
    }
}

#[derive(Debug)]
pub struct Upsampler<const N: usize> {
    buf: BitMaskRB<f32>,
    pos: isize,
    coef: [f32; N],
}

impl<const N: usize> Upsampler<N> {
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

    #[rustfmt::skip]
    pub fn process_sample(&mut self, s: f32) -> (f32, f32) {
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

        (s1, s2)
    }

    pub fn process_block(&mut self, input: &[f32], output: &mut [f32]) {
        assert_eq!(
            input.len() * 2,
            output.len(),
            "Input must be twice the size of output."
        );

        for (i, chunk) in output.chunks_exact_mut(2).enumerate() {
            let (even, odd) = self.process_sample(input[i]);
            chunk[0] = even;
            chunk[1] = odd;
        }
    }

    pub fn clear(&mut self) {
        self.pos = 0;
        self.buf.raw_data_mut().fill(0.);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Test input (upsampled signal)
    const INPUT_DOWN: [f32; 32] = [
        1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
        1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
    ];

    // Expected output downsampled
    const EXPECTED_DOWN: [f32; 16] = [
        -4.9631526e-5,
        0.00069190684,
        -0.003426439,
        0.011447029,
        -0.030675305,
        0.07221339,
        -0.16344123,
        0.4764971,
        0.33655876,
        -0.42778662,
        0.4693247,
        -0.48855296,
        0.49657357,
        -0.49930808,
        0.49995038,
        -0.5,
    ];

    const INPUT_UP: [f32; 16] = [
        1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0,
    ];

    const EXPECTED_UP: [f32; 32] = [
        -9.926305e-5,
        0.0,
        0.0011852876,
        0.0,
        -0.0040852507,
        0.0,
        0.009386828,
        0.0,
        -0.018330121,
        0.0,
        0.035232782,
        0.0,
        -0.08104939,
        0.0,
        0.40842336,
        1.0,
        1.333273,
        1.0,
        0.035232782,
        -1.0,
        -1.4326525,
        -1.0,
        0.009386828,
        1.0,
        1.4102371,
        1.0,
        0.0011852876,
        -1.0,
        -1.4144217,
        -1.0,
        0.0,
        1.0,
    ];

    #[test]
    fn test_downsampler() {
        let mut downsampler = Downsampler::<{ COEF_31.len() }>::new(&COEF_31);

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
        let mut upsampler = Upsampler::<{ COEF_31.len() }>::new(&COEF_31);

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
        let mut downsampler = Downsampler::<{ COEF_51.len() }>::new(&COEF_51);

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
        let mut upsampler = Upsampler::<{ COEF_19.len() }>::new(&COEF_19);

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
