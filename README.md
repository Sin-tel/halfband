# halfband

A Rust library for resampling audio signals. Provides 2x up- and downsampling with an efficient polyphase design. These can be cascaded to get 4x, 8x etc.

Both FIR (windowed sinc) and IIR implementation are available. The IIR implementation is based on the [HIIR library by Laurent De Soras](https://ldesoras.fr/prod.html#src_audio), *"an oversampling and Hilbert transform library in C++"*.

