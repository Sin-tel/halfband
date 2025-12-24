//! Example of how to calculate the phase delay of an FIR up/down cascade

#![allow(dead_code)]
#![allow(unused)]
mod util;

use crate::util::{measure_fractional_delay, save_wav};
use halfband::fir;
use std::f32::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample_rate = 44100.0;
    let freq_hz = 800.0;
    let total_samples = 5_000;

    let input: Vec<f32> = (0..total_samples)
        .map(|i| 0.5 * (2.0 * PI * freq_hz * i as f32 / sample_rate).sin())
        .collect();

    // Mix some different presets to be sure.
    let mut downsampler1 = fir::Downsampler5::default();
    let mut upsampler1 = fir::Upsampler8::default();
    let mut downsampler2 = fir::Downsampler12::default();
    let mut upsampler2 = fir::Upsampler16::default();

    // Up
    let mut upsampled1 = vec![0.0; input.len() * 2];
    let mut upsampled2 = vec![0.0; input.len() * 4];
    upsampler1.process_block(&input, &mut upsampled1);
    upsampler2.process_block(&upsampled1, &mut upsampled2);

    // Down
    let mut downsampled = vec![0.0; input.len() * 2];
    let mut output = vec![0.0; input.len()];
    downsampler2.process_block(&upsampled2, &mut downsampled);
    downsampler1.process_block(&downsampled, &mut output);

    assert_eq!(input.len(), output.len());

    let delay = measure_fractional_delay(freq_hz, sample_rate, &output);

    println!("Observed delay: {}", delay);

    // First stage runs at 2x
    let delay1 = upsampler1.latency(2) + downsampler1.latency(2);
    // Second stage runs at 4x
    let delay2 = upsampler2.latency(4) + downsampler2.latency(4);

    println!("Computed delay: {}", delay1 + delay2);

    // save_wav("input.wav", &input, sample_rate as u32)?;
    // save_wav("output.wav", &output, sample_rate as u32)?;

    Ok(())
}
