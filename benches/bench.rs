use criterion::{Criterion, criterion_group, criterion_main};
use halfband::fir::COEF_31;
use std::hint::black_box;

const TEST_COEFFICIENTS: [f32; 8] = [
    0.07711508, 0.26596853, 0.48207062, 0.66510415, 0.79682046, 0.8841015, 0.94125146, 0.9820054,
];

const N: usize = 1024;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("downsample_iir", |b| {
        let input = vec![0.0; N];
        let mut output = vec![0.0; N / 2];

        let mut downsampler = halfband::iir::Downsampler::<4>::new(&TEST_COEFFICIENTS);

        b.iter(|| {
            downsampler.process_block(black_box(&input), &mut output);
        })
    });

    c.bench_function("upsample_iir", |b| {
        let input = vec![0.0; N / 2];
        let mut output = vec![0.0; N];

        let mut upsampler = halfband::iir::Upsampler::<4>::new(&TEST_COEFFICIENTS);

        b.iter(|| {
            upsampler.process_block(black_box(&input), &mut output);
        })
    });

    c.bench_function("downsample_fir", |b| {
        let input = vec![0.0; N];
        let mut output = vec![0.0; N / 2];

        let mut downsampler = halfband::fir::Downsampler::<8>::new(&COEF_31);

        b.iter(|| {
            downsampler.process_block(black_box(&input), &mut output);
        })
    });

    c.bench_function("upsample_fir", |b| {
        let input = vec![0.0; N / 2];
        let mut output = vec![0.0; N];

        let mut upsampler = halfband::fir::Upsampler::<8>::new(&COEF_31);

        b.iter(|| {
            upsampler.process_block(black_box(&input), &mut output);
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
