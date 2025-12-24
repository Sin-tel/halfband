//! A high-performance 2x resampling library for audio signals using polyphase designs.
//!
//! This crate provides two main approaches to resampling:
//! - [Finite Impulse Response (FIR)](fir): Linear-phase windowed-sinc filters. These filters
//!   have a constant phase delay (latency) for all frequencies.
//! - [Infinite Impulse Response (IIR)](iir): All-pass based recursive filters based on
//!   the [HIIR library](https://ldesoras.fr/prod.html#src_audio). These generally provide
//!   steeper filtering, lower latency and better performance, but the phase delay is not
//!   constant across the frequency range.
//!
//! Both FIR and IIR provide `Upsampler` and `Downsampler` structs that have the same API.
//!
//! ## Example
//!
//! ```rust
//! use halfband::fir;
//!
//! let mut up = fir::Upsampler8::default();
//! let mut down = fir::Downsampler8::default();
//!
//! let mut samples = vec![0.1, 0.2, 0.3];
//!
//! for s in &mut samples {
//!     // Upsample: returns two high-rate samples
//!     let [mut s0, mut s1] = up.process(*s);
//!
//!     // Apply nonlinearity (e.g., hard clipping or tanh)
//!     s0 = s0.tanh();
//!     s1 = s1.tanh();
//!
//!     // Downsample: filters and returns one low-rate sample
//!     *s = down.process(s0, s1);
//! }
//! ```

#![warn(missing_docs, missing_debug_implementations, rust_2018_idioms)]
#![warn(clippy::cast_lossless)]

pub mod fir;
pub mod iir;
