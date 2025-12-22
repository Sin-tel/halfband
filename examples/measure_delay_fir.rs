//! Example of how to calculate the phase delay of an FIR up/down pair

#![allow(dead_code)]
#![allow(unused)]
mod util;

use crate::util::{measure_fractional_delay, save_wav};
use halfband::fir::presets::*;
use std::f32::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample_rate = 44100.0;
    // Since FIR phase delay is flat we can use lower frequencies
    let freq_hz = 800.0;
    let total_samples = 5_000;

    let input: Vec<f32> = (0..total_samples)
        .map(|i| 0.5 * (2.0 * PI * freq_hz * i as f32 / sample_rate).sin())
        .collect();

    // Setup stages
    let mut upsampler = Upsampler31::default();
    let mut downsampler = Downsampler19::default();

    let mut upsampled = vec![0.0; input.len() * 2];

    upsampler.process_block(&input, &mut upsampled);

    let mut output = vec![0.0; input.len()];
    downsampler.process_block(&upsampled, &mut output);

    let delay = measure_fractional_delay(freq_hz, sample_rate, &output);

    println!("Observed delay: {}", delay);

    let up_delay = upsampler.get_latency();
    let down_delay = downsampler.get_latency();
    println!("Computed delay: {}", (up_delay + down_delay));

    // save_wav("input.wav", &input, sample_rate as u32)?;
    // save_wav("output.wav", &output, sample_rate as u32)?;

    Ok(())
}
