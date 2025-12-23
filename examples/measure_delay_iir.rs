//! Example of how to calculate the phase delay of an IIR up/down pair

#![allow(dead_code)]
#![allow(unused)]
mod util;

use crate::util::{measure_fractional_delay, save_wav};
use halfband::iir;
use halfband::iir::design::{coefs_transition, phase_delay};
use std::f32::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample_rate = 44100.0;
    let freq_hz = 4000.0;
    let total_samples = 5_000;

    let input: Vec<f32> = (0..total_samples)
        .map(|i| 0.5 * (2.0 * PI * freq_hz * i as f32 / sample_rate).sin())
        .collect();

    // Setup stages
    // This pair has an expected latency of ~5.0 samples
    let coefs = coefs_transition(12, 0.0472053);
    let mut upsampler = iir::Upsampler12::new(&coefs);
    let mut downsampler = iir::Downsampler12::new(&coefs);

    let mut upsampled = vec![0.0; input.len() * 2];

    upsampler.process_block(&input, &mut upsampled);

    let mut output = vec![0.0; input.len()];
    downsampler.process_block(&upsampled, &mut output);

    let delay = measure_fractional_delay(freq_hz, sample_rate, &output);

    println!("Observed delay: {}", delay);

    // The filters operate relative to the full-band sample rate (2*fs).
    // Then there's an additional 1/2 sample advance that we have to compensate for.
    // Finally, the total at half-band is double that, but we have two stages (up and down) so that cancels out.
    let f_relative = freq_hz / (2.0 * sample_rate);
    let delay = phase_delay(&coefs, f_relative.into()) - 0.5;

    println!("Computed delay: {}", delay);

    // save_wav("input.wav", &input, sample_rate as u32)?;
    // save_wav("output.wav", &output, sample_rate as u32)?;

    Ok(())
}
