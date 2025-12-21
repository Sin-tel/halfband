# Compute FIR coefficients for sinc filter with Hamming window

import numpy as np

k_taps = 16
ntaps = k_taps * 4 - 1

print("n taps", ntaps)

N = ntaps - 1

# we do a divide by zero but fix it later
np.seterr(divide="ignore", invalid="ignore")

if ntaps % 2 == 0:
    print("ERROR ntaps must be odd")

w = np.hamming(ntaps)
x = np.linspace(-N / 2, N / 2, ntaps)

sinc = np.sin(np.pi * x / 2) / (np.pi * x)

# fix divide by zero
sinc[N // 2] = 1 / 2

out = np.multiply(w, sinc)

out_even = out[0 : N // 2 : 2]

out_even = 0.5 * out_even / np.sum(out_even)

print("sum:", np.sum(out_even))

out_even = out_even.astype(np.float32)

print(f"pub const COEF_{ntaps}: [f32; {k_taps}] = ")
print(
    np.array2string(out_even, separator=", ", floatmode="unique", max_line_width=5)
    + ";"
)
