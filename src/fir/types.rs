//! Type aliases resamplers of commonly used filter sizes.

use crate::fir::{Downsampler, Upsampler};

macro_rules! generate_type_aliases {
    ($($taps:literal : $n:literal),*) => {
        paste::paste! {
            $(
                #[doc = "Convenience type for a 2x downsampler with " $n " coefficients."]
                #[doc = "\n\nEquivalent to a " $taps "-tap half-band FIR filter."]
                pub type [<Downsampler $n>] = Downsampler<$n>;
                #[doc = "Convenience type for a 2x upsampler with " $n " coefficients."]
                #[doc = "\n\nEquivalent to a " $taps "-tap half-band FIR filter."]
                pub type [<Upsampler $n>] = Upsampler<$n>;
            )*
        }
    };
}

generate_type_aliases!(15:4, 19:5, 23:6, 27:7, 31:8, 35:9, 39:10, 43:11, 47:12, 63:16);
