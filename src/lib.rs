//! A high-performance 2x resampling library for audio signals using polyphase designs.
//!
//! This crate provides two main approaches to resampling:
//! - Finite Impulse Response (FIR): Linear-phase windowed-sinc filters. Ideal when
//!   phase preservation is critical. Includes presets for 19, 31, 47, and 63 taps.
//! - Infinite Impulse Response (IIR): All-pass based recursive filters based on
//!   the [HIIR library](https://ldesoras.fr/prod.html#src_audio). Ideal for steep
//!   transition bands.

#![warn(missing_docs, missing_debug_implementations, rust_2018_idioms)]
#![warn(clippy::cast_lossless)]

pub mod fir;
pub mod iir;
