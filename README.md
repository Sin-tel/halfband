[![build](https://github.com/Sin-tel/halfband/actions/workflows/rust.yml/badge.svg)](https://github.com/Sin-tel/halfband/actions/workflows/rust.yml)
[![Crates.io](https://img.shields.io/crates/v/halfband.svg)](https://crates.io/crates/halfband) 
[![docs.rs](https://docs.rs/halfband/badge.svg)](https://docs.rs/halfband/)

# halfband

A Rust library for resampling audio signals. Provides 2x up- and downsampling with an efficient polyphase design. These can be cascaded to get 4x, 8x, etc.

Both FIR (windowed sinc) and IIR implementation are available. The IIR implementation is based on the [HIIR library by Laurent De Soras](https://ldesoras.fr/prod.html#src_audio), *"an oversampling and Hilbert transform library in C++"*.

See `examples/` for how to do 4x oversampling.

Currently does not use any SIMD intrinsics, but code is written in a way that allows for decent auto-vectorization. Benchmarks show it is plenty fast: a single stage takes ~1 μs to process a buffer of 1024 samples.

## TODO

* Provide presets for IIR (see `hiir/oversampling.txt`)
* Implement IIR Hilbert filters
