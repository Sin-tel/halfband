use crate::iir::Downsampler;
use crate::iir::PolyphaseEven;
use crate::iir::PolyphaseOdd;
use crate::iir::Upsampler;

macro_rules! generate_type_aliases {
    // Even: count = 2 * N
    (even $($count:literal : $n:literal),*) => {
        paste::paste! {
            $(
                #[doc = "Convenience type for a 2x downsampler with " $count " coefficients."]
                pub type [<Downsampler $count>] = Downsampler<PolyphaseEven<$n>>;
                #[doc = "Convenience type for a 2x upsampler with " $count " coefficients."]
                pub type [<Upsampler $count>] = Upsampler<PolyphaseEven<$n>>;
            )*
        }
    };
    // Odd: count = 2 * N - 1
    (odd $($count:literal : $n:literal),*) => {
        paste::paste! {
            $(
                #[doc = "Convenience type for a 2x downsampler with " $count " coefficients."]
                pub type [<Downsampler $count>] = Downsampler<PolyphaseOdd<$n>>;
                #[doc = "Convenience type for a 2x upsampler with " $count " coefficients."]
                pub type [<Upsampler $count>] = Upsampler<PolyphaseOdd<$n>>;
            )*
        }
    };
}

generate_type_aliases!(even 2:1, 4:2, 6:3, 8:4, 10:5, 12:6, 14:7, 16:8);
generate_type_aliases!(odd 3:2, 5:3, 7:4, 9:5, 11:6, 13:7, 15:8);
