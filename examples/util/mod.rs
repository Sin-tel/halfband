use hound::{WavSpec, WavWriter};
use std::f64::consts::PI;

pub fn generate_sine_sweep(
    duration_secs: f64,
    sample_rate: u32,
    start_freq: f64,
    end_freq: f64,
) -> Vec<f32> {
    // Internal calculation in f64 for a cleaner result
    let num_samples = (duration_secs * sample_rate as f64) as usize;
    let mut sweep = Vec::with_capacity(num_samples);

    for i in 0..num_samples {
        let t = (i as f64) / (sample_rate as f64);
        let t_normalized = t / duration_secs;

        let phase = 2.0 * PI * (start_freq * t + (end_freq - start_freq) * t * t_normalized * 0.5);
        sweep.push((phase).sin() as f32);
    }

    sweep
}

pub fn save_wav(filename: &str, samples: &[f32], sample_rate: u32) -> Result<(), hound::Error> {
    let spec = WavSpec {
        channels: 1,
        sample_rate,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = WavWriter::create(filename, spec)?;

    for &sample in samples {
        // Convert f32 [-1.0, 1.0] to i16
        let scaled = (sample * (i16::MAX as f32)) as i16;
        writer.write_sample(scaled)?;
    }

    writer.finalize()?;
    println!("Saved: {}", filename);
    Ok(())
}

pub fn measure_fractional_delay(freq_hz: f32, sample_rate: f32, output: &[f32]) -> f32 {
    // Calculate offset of a sinewave in samples, using quadrature cross-correlation.
    // Assumes original input starts at a phase of zero!

    let total_samples = output.len();

    // Ignore transient response
    let skip_samples = total_samples / 5;

    let omega = 2.0 * std::f32::consts::PI * freq_hz / sample_rate;

    let mut sum_cos = 0.0;
    let mut sum_sin = 0.0;

    // Calculate cross-correlation with quadrature
    for (i, y) in output.into_iter().enumerate() {
        if i >= skip_samples {
            let ref_cos = (omega * i as f32).cos();
            let ref_sin = (omega * i as f32).sin();
            sum_cos += y * ref_cos;
            sum_sin += y * ref_sin;
        }
    }

    // Phase shift in radians
    let phase_shift = f32::atan2(sum_cos, sum_sin);

    // Convert phase shift to samples:
    // Delay = - phase / omega
    let mut delay_samples = -phase_shift / omega;

    let period = sample_rate / freq_hz;
    while delay_samples < -0.5 {
        delay_samples += period;
    }

    delay_samples
}
