[![build](https://github.com/Sin-tel/halfband/actions/workflows/rust.yml/badge.svg)](https://github.com/Sin-tel/halfband/actions/workflows/rust.yml)
[![Crates.io](https://img.shields.io/crates/v/halfband.svg)](https://crates.io/crates/halfband) 
[![docs.rs](https://docs.rs/halfband/badge.svg)](https://docs.rs/halfband/)

# halfband

A Rust library for resampling audio signals. Provides 2x up- and downsampling with an efficient polyphase design. These can be cascaded to get 4x, 8x, etc.

Both FIR (windowed sinc) and IIR implementation are available. The IIR implementation is based on the [HIIR library by Laurent De Soras](https://ldesoras.fr/prod.html#src_audio), *"an oversampling and Hilbert transform library in C++"*.

See [examples/clipper_fir.rs](./examples/clipper_fir.rs) and [examples/clipper_iir.rs](./examples/clipper_iir.rs) for a complete examples that do 4x oversampling.

Currently does not use any SIMD intrinsics, but code is written in a way that allows for decent auto-vectorization. Benchmarks show it is plenty fast: a single stage takes ~1 μs to process a buffer of 1024 samples.

## Example

```rust
use halfband::fir;

let mut up = fir::Upsampler8::default();
let mut down = fir::Downsampler8::default();

let mut samples = vec![0.1, 0.2, 0.3];

for s in &mut samples {
    // Upsample: returns two high-rate samples
    let [mut s0, mut s1] = up.process(*s);

    // Apply nonlinearity
    s0 = s0.tanh();
    s1 = s1.tanh();

    // Downsample: filters and returns one low-rate sample
    *s = down.process(s0, s1);
}
```
