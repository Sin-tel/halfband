//! Type aliases resamplers of commonly used filter sizes.

use crate::iir::Downsampler;
use crate::iir::PolyphaseEven;
use crate::iir::PolyphaseOdd;
use crate::iir::Upsampler;
use crate::iir::coefs_transition;

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

// Default implementations

impl Default for Upsampler12 {
    fn default() -> Self {
        let coefs = coefs_transition(12, 0.0472053);
        Self::new(&coefs)
    }
}

impl Default for Downsampler12 {
    fn default() -> Self {
        let coefs = coefs_transition(12, 0.0472053);
        Self::new(&coefs)
    }
}

impl Default for Upsampler10 {
    fn default() -> Self {
        let coefs = coefs_transition(10, 0.0416368);
        Self::new(&coefs)
    }
}

impl Default for Downsampler10 {
    fn default() -> Self {
        let coefs = coefs_transition(10, 0.0416368);
        Self::new(&coefs)
    }
}

impl Default for Upsampler8 {
    fn default() -> Self {
        let coefs = coefs_transition(8, 0.0343747);
        Self::new(&coefs)
    }
}

impl Default for Downsampler8 {
    fn default() -> Self {
        let coefs = coefs_transition(8, 0.0343747);
        Self::new(&coefs)
    }
}

impl Default for Upsampler6 {
    fn default() -> Self {
        let coefs = coefs_transition(6, 0.0246851);
        Self::new(&coefs)
    }
}

impl Default for Downsampler6 {
    fn default() -> Self {
        let coefs = coefs_transition(6, 0.0246851);
        Self::new(&coefs)
    }
}

#[cfg(test)]
mod tests {
    use crate::iir::*;
    use approx::assert_relative_eq;
    // Default impls
    #[test]
    fn test_default_12() {
        let up = Downsampler12::default();
        let down = Upsampler12::default();
        assert_relative_eq!(up.latency(2) + down.latency(2), 5.0, epsilon = 1e-6);
    }

    #[test]
    fn test_default_10() {
        let up = Downsampler10::default();
        let down = Upsampler10::default();
        assert_relative_eq!(up.latency(2) + down.latency(2), 4.0, epsilon = 1e-6);
    }

    #[test]
    fn test_default_8() {
        let up = Downsampler8::default();
        let down = Upsampler8::default();
        assert_relative_eq!(up.latency(2) + down.latency(2), 3.0, epsilon = 1e-6);
    }

    #[test]
    fn test_default_6() {
        let up = Downsampler6::default();
        let down = Upsampler6::default();
        assert_relative_eq!(up.latency(2) + down.latency(2), 2.0, epsilon = 1e-6);
    }
}
