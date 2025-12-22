//! Example of how to calculate the phase delay of an IIR up/down cascade

#![allow(dead_code)]
#![allow(unused)]
mod util;

use crate::util::{measure_fractional_delay, save_wav};
use halfband::iir::design::{compute_coefs_tbw, compute_phase_delay};
use halfband::iir::{Downsampler, Upsampler};
use std::f32::consts::PI;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample_rate = 44100.0;
    let freq_hz = 4000.0;
    let total_samples = 5_000;

    let input: Vec<f32> = (0..total_samples)
        .map(|i| 0.5 * (2.0 * PI * freq_hz * i as f32 / sample_rate).sin())
        .collect();

    // Setup stages
    // This cascade has an expected latency of ~5.0 samples
    const N1: usize = 5;
    const N2: usize = 2;
    let coefs1 = compute_coefs_tbw(2 * N1, 0.0367598);
    let coefs2 = compute_coefs_tbw(2 * N2, 0.261666);
    let mut downsampler1 = Downsampler::<N1>::new(&coefs1);
    let mut upsampler1 = Upsampler::<N1>::new(&coefs1);
    let mut downsampler2 = Downsampler::<N2>::new(&coefs2);
    let mut upsampler2 = Upsampler::<N2>::new(&coefs2);

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

    // First filter operates at 2*fs
    let f1 = freq_hz / (2.0 * sample_rate);
    let delay1 = compute_phase_delay(&coefs1, f1.into()) - 0.5;

    // Second filter operates at 4*fs
    let f2 = freq_hz / (4.0 * sample_rate);
    let delay2 = compute_phase_delay(&coefs2, f2.into()) - 0.5;

    // Second stage contributes only half
    println!("Computed delay: {}", delay1 + delay2 * 0.5);

    // save_wav("input.wav", &input, sample_rate as u32)?;
    // save_wav("output.wav", &output, sample_rate as u32)?;

    Ok(())
}
