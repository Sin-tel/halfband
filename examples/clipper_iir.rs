//! Sine sweep through a soft clipper with 4x oversampling, using IIR filters.
//! Run with:
//! `cargo run --example clipper_iir`
//!
//! Generates `clipped_naive.wav` vs `clipped_oversampled.wav`.
//! The oversampled version should have less aliasing.

#![allow(dead_code)]
#![allow(unused)]
mod util;

use crate::util::{generate_sine_sweep, save_wav};
use halfband::iir;
use halfband::iir::design::coefs_transition;

fn softclip(x: f32) -> f32 {
    (x * 4.0).tanh() / 2.0
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sample_rate = 44100;
    let duration = 5.0;
    let start_freq = 20.0;
    let end_freq = sample_rate as f64 / 2.0;

    let original_sweep = generate_sine_sweep(duration, sample_rate, start_freq, end_freq);

    let clipped: Vec<f32> = original_sweep.iter().map(|x| softclip(*x)).collect();

    save_wav("clipped_naive.wav", &clipped, sample_rate)?;

    // Setup stages
    let coefs1 = coefs_transition(10, 0.0367598);
    let coefs2 = coefs_transition(4, 0.261666);
    let mut downsampler1 = iir::Downsampler10::new(&coefs1);
    let mut upsampler1 = iir::Upsampler10::new(&coefs1);
    let mut downsampler2 = iir::Downsampler4::new(&coefs2);
    let mut upsampler2 = iir::Upsampler4::new(&coefs2);

    // Up
    let mut upsampled1 = vec![0.0; original_sweep.len() * 2];
    let mut upsampled2 = vec![0.0; original_sweep.len() * 4];
    upsampler1.process_block(&original_sweep, &mut upsampled1);
    upsampler2.process_block(&upsampled1, &mut upsampled2);

    // Run clipper
    let clipped_up: Vec<f32> = upsampled2.iter().map(|x| softclip(*x)).collect();

    // Down
    let mut downsampled2 = vec![0.0; original_sweep.len() * 2];
    let mut downsampled1 = vec![0.0; original_sweep.len()];
    downsampler2.process_block(&clipped_up, &mut downsampled2);
    downsampler1.process_block(&downsampled2, &mut downsampled1);

    // Save downsampled signal
    save_wav("clipped_oversampled.wav", &downsampled1, sample_rate)?;

    Ok(())
}
