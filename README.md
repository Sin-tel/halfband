[![build](https://github.com/Sin-tel/halfband/actions/workflows/rust.yml/badge.svg)](https://github.com/Sin-tel/halfband/actions/workflows/rust.yml)
[![Crates.io](https://img.shields.io/crates/v/halfband.svg)](https://crates.io/crates/halfband) 
[![docs.rs](https://docs.rs/halfband/badge.svg)](https://docs.rs/halfband/)

# halfband

A Rust library for resampling audio signals. Provides 2x up- and downsampling with an efficient polyphase design. These can be cascaded to get 4x, 8x etc.

Both FIR (windowed sinc) and IIR implementation are available. The IIR implementation is based on the [HIIR library by Laurent De Soras](https://ldesoras.fr/prod.html#src_audio), *"an oversampling and Hilbert transform library in C++"*.

See `examples/` for how to do 4x oversampling.

The FIR comes with preset coefficients for 19, 31, 47 and 63 taps.

## TODO

* Provide presets for IIR (see `hiir/oversampling.txt`)
* Provide methods to compute total phase delay in samples
* Currently IIR only works with an even number of coefficients
