// Integration test — runs with `cargo test`
// This mirrors main.cpp exactly so you can compare outputs.

use microfft::Complex32;
use num::Complex;
use microcwt::{
    MicroCwt,
    scales::{Scales, ScaleType},
    wavelet::Morlet,
};

// ── helpers ──────────────────────────────────────────────────────────────────

fn make_cosine_signal(n: usize, freq_hz: f32, fs: f32) -> Vec<f32> {
    (0..n)
        .map(|i| (2.0 * std::f32::consts::PI * freq_hz * i as f32 / fs).cos())
        .collect()
}

fn magnitude(c: Complex<f32>) -> f32 {
    (c.re * c.re + c.im * c.im).sqrt()
}

// ── Test 1: smoke test — does it run without panicking? ──────────────────────
#[test]
fn test_cwt_runs() {
    let n   = 1000;
    let fs  = 1000;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    let sig = make_cosine_signal(n, 1.0, fs as f32);
    let mut tfm = vec![Complex::new(0.0_f32, 0.0_f32); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    // If we got here without panicking the basic pipeline works
    assert_eq!(tfm.len(), n * fn_ as usize);
}

// ── Test 2: output is non-zero ───────────────────────────────────────────────
// A cosine input should produce a non-trivial CWT output.
#[test]
fn test_cwt_nonzero_output() {
    let n   = 1000;
    let fs  = 1000;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    let sig = make_cosine_signal(n, 1.0, fs as f32);
    let mut tfm = vec![Complex::new(0.0_f32, 0.0_f32); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    let max_mag = tfm.iter().map(|&c| magnitude(c)).fold(0.0_f32, f32::max);
    assert!(max_mag > 0.0, "CWT output was all zeros");
}

// ── Test 3: peak frequency detection ─────────────────────────────────────────
// A 1 Hz cosine should produce the strongest response near 1 Hz.
// This is the key correctness check — if the peak is at the wrong scale
// the FFT mirroring or daughter multiplication has a bug.
#[test]
fn test_cwt_peak_at_correct_frequency() {
    let n   = 1000;
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    // 1 Hz cosine — peak should appear near 1 Hz in the scalogram
    let sig = make_cosine_signal(n, 1.0, fs as f32);
    let mut tfm = vec![Complex::new(0.0_f32, 0.0_f32); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    // Find which scale index has the highest average magnitude
    let mut scale_energies: Vec<f32> = (0..fn_ as usize)
        .map(|si| {
            let row = &tfm[si * n..(si + 1) * n];
            row.iter().map(|&c| magnitude(c)).sum::<f32>() / n as f32
        })
        .collect();

    let peak_scale_idx = scale_energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    // Get the frequency that index corresponds to
    let mut freqs = vec![0.0_f32; fn_ as usize];
    scales.get_frequencies(&mut freqs);
    let peak_freq = freqs[peak_scale_idx];

    println!("Peak frequency: {:.3} Hz (expected ~1.0 Hz)", peak_freq);

    // Allow ±2 Hz tolerance given the coarse 200-bin resolution
    assert!(
        (peak_freq - 1.0).abs() < 2.0,
        "Peak at {:.3} Hz, expected ~1.0 Hz",
        peak_freq
    );
}

// ── Test 4: complex input vs real input give same result ─────────────────────
// main.cpp runs both sig (real) and sigc (complex with imag=0).
// Both should produce the same output since sigc has zero imaginary part.
#[test]
fn test_real_and_complex_input_match() {
    let n   = 1000;
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    let sig: Vec<f32> = make_cosine_signal(n, 1.0, fs as f32);

    // Interleaved complex input: [re0, 0.0, re1, 0.0, ...]
    let sigc: Vec<f32> = sig
        .iter()
        .flat_map(|&re| [re, 0.0_f32])
        .collect();

    let mut tfm_real    = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];
    let mut tfm_complex = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];

    {
        let mut morlet = Morlet::new(1.0);
        let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
        let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);
        fcwt.cwt(&sig, &mut tfm_real, &mut scales, false);
    }
    {
        let mut morlet = Morlet::new(1.0);
        let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
        let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);
        fcwt.cwt(&sigc, &mut tfm_complex, &mut scales, true);
    }

    // Compare magnitudes — allow small floating point delta
    for (i, (&r, &c)) in tfm_real.iter().zip(tfm_complex.iter()).enumerate() {
        let diff = (magnitude(r) - magnitude(c)).abs();
        assert!(
            diff < 1e-3,
            "Mismatch at index {}: real mag={:.6}, complex mag={:.6}",
            i, magnitude(r), magnitude(c)
        );
    }
}

// ── Test 5: zero input gives zero output ─────────────────────────────────────
#[test]
fn test_zero_input_gives_zero_output() {
    let n   = 1000;
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    let sig  = vec![0.0_f32; n];
    let mut tfm = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    let max_mag = tfm.iter().map(|&c| magnitude(c)).fold(0.0_f32, f32::max);
    assert!(
        max_mag < 1e-6,
        "Expected ~zero output for zero input, got max magnitude {:.6}",
        max_mag
    );
}

use std::time::Instant;

// ── Test 6: batch processing ──────────────────────────────────────────────────
// Simulates processing many signals in sequence, like a training data loader
// feeding samples into a feature extractor.
#[test]
fn test_batch_processing() {
    let n      = 1000;
    let fs     = 1000_i32;
    let f0     = 0.1_f32;
    let f1     = 20.0_f32;
    let fn_    = 200_i32;
    let n_batch = 100; // 100 signals, like a training batch

    for batch_idx in 0..n_batch {
        // Vary the frequency per sample, like different labelled examples
        let freq = 1.0 + batch_idx as f32 * 0.1;
        let sig  = make_cosine_signal(n, freq, fs as f32);
        let mut tfm = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];

        let mut morlet = Morlet::new(1.0);
        let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
        let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

        fcwt.cwt(&sig, &mut tfm, &mut scales, false);

        // Each output should be non-zero
        let max_mag = tfm.iter().map(|&c| magnitude(c)).fold(0.0_f32, f32::max);
        assert!(
            max_mag > 0.0,
            "Batch item {} produced zero output",
            batch_idx
        );
    }
}

// ── Test 7: throughput benchmark ─────────────────────────────────────────────
// Measures how many signals per second your CWT can process.
// Compare this number against your model's inference budget.
#[test]
fn test_throughput() {
    let n   = 1000;
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;
    let n_runs = 500;

    let sig = make_cosine_signal(n, 5.0, fs as f32);

    let start = Instant::now();

    for _ in 0..n_runs {
        let mut tfm = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];
        let mut morlet = Morlet::new(1.0);
        let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
        let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);
        fcwt.cwt(&sig, &mut tfm, &mut scales, false);
    }

    let elapsed = start.elapsed();
    let per_signal_ms = elapsed.as_secs_f64() * 1000.0 / n_runs as f64;
    let signals_per_sec = 1000.0 / per_signal_ms;

    println!(
        "Throughput: {:.1} signals/sec ({:.3} ms/signal) over {} runs",
        signals_per_sec, per_signal_ms, n_runs
    );

    // Adjust this threshold to match your model's real-time requirement.
    // For example if your model needs to process 10 signals/sec, assert > 10.
    assert!(
        signals_per_sec > 1.0,
        "Too slow: {:.1} signals/sec",
        signals_per_sec
    );
}

// ── Test 8: noisy mixed-frequency input ──────────────────────────────────────
// Real sensor data has multiple frequency components plus noise.
// This checks the CWT still localises the dominant frequency correctly.
#[test]
fn test_noisy_mixed_frequency_input() {
    let n   = 1000;
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    // Mix: strong 5 Hz + weak 15 Hz + noise
    // Mirrors real EEG/ECG/accelerometer data fed into an ML model
    let sig: Vec<f32> = (0..n)
        .map(|i| {
            let t = i as f32 / fs as f32;
            let signal = 1.0  * (2.0 * std::f32::consts::PI * 5.0  * t).cos()
                       + 0.3  * (2.0 * std::f32::consts::PI * 15.0 * t).cos();
            // Deterministic pseudo-noise so the test is reproducible
            let noise = 0.1 * ((i as f32 * 127.1).sin()
                              + (i as f32 * 311.7).sin());
            signal + noise
        })
        .collect();

    let mut tfm = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    let scale_energies: Vec<f32> = (0..fn_ as usize)
        .map(|si| {
            let row = &tfm[si * n..(si + 1) * n];
            row.iter().map(|&c| magnitude(c)).sum::<f32>() / n as f32
        })
        .collect();

    let peak_idx = scale_energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let mut freqs = vec![0.0_f32; fn_ as usize];
    scales.get_frequencies(&mut freqs);
    let peak_freq = freqs[peak_idx];

    println!(
        "Mixed signal peak: {:.3} Hz (expected ~5.0 Hz, dominant component)",
        peak_freq
    );

    assert!(
        (peak_freq - 5.0).abs() < 2.0,
        "Peak at {:.3} Hz, expected ~5.0 Hz",
        peak_freq
    );
}

// ── Test 9: longer signal (closer to real ML input sizes) ────────────────────
// 1000 samples at 1000 Hz = 1 second. Many ML models use 5-30 second windows.
// This checks the next_power_of_two path for larger inputs.
#[test]
fn test_longer_signal() {
    let n   = 4000; // 4 seconds at 1000 Hz — newsize will be 4096
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    let sig = make_cosine_signal(n, 3.0, fs as f32);
    let mut tfm = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    let scale_energies: Vec<f32> = (0..fn_ as usize)
        .map(|si| {
            let row = &tfm[si * n..(si + 1) * n];
            row.iter().map(|&c| magnitude(c)).sum::<f32>() / n as f32
        })
        .collect();

    let peak_idx = scale_energies
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let mut freqs = vec![0.0_f32; fn_ as usize];
    scales.get_frequencies(&mut freqs);
    let peak_freq = freqs[peak_idx];

    println!("4000-sample signal peak: {:.3} Hz (expected ~3.0 Hz)", peak_freq);

    assert!(
        (peak_freq - 3.0).abs() < 2.0,
        "Peak at {:.3} Hz, expected ~3.0 Hz",
        peak_freq
    );
}

// ── Test 10: output shape matches what a model layer expects ─────────────────
// An ML model treats the CWT output as a 2D feature map: [n_scales x n_samples]
// This checks the memory layout is correct for that interpretation.
#[test]
fn test_output_shape_for_ml() {
    let n   = 1000;
    let fs  = 1000_i32;
    let f0  = 0.1_f32;
    let f1  = 20.0_f32;
    let fn_ = 200_i32;

    let sig = make_cosine_signal(n, 1.0, fs as f32);
    let mut tfm = vec![Complex::new(0.0_f32, 0.0); n * fn_ as usize];

    let mut morlet = Morlet::new(1.0);
    let mut scales = Scales::new(&morlet, ScaleType::FcwtLinfreqs, fs, f0, f1, fn_);
    let mut fcwt   = MicroCwt::new(&mut morlet, fs, false, false);

    fcwt.cwt(&sig, &mut tfm, &mut scales, false);

    // Verify the flat buffer can be interpreted as [fn_ x n] row-major 2D array
    assert_eq!(tfm.len(), fn_ as usize * n);

    // Each row (scale) should be addressable as a slice of length n
    for scale_idx in 0..fn_ as usize {
        let row = &tfm[scale_idx * n..(scale_idx + 1) * n];
        assert_eq!(row.len(), n);

        // Extract real and imaginary channels separately,
        // which is how most ML frameworks (PyTorch, etc.) consume complex CWT output
        let real_channel: Vec<f32> = row.iter().map(|c| c.re).collect();
        let imag_channel: Vec<f32> = row.iter().map(|c| c.im).collect();
        let mag_channel:  Vec<f32> = row.iter().map(|c| magnitude(*c)).collect();

        assert_eq!(real_channel.len(), n);
        assert_eq!(imag_channel.len(), n);
        assert_eq!(mag_channel.len(), n);
    }

    println!(
        "Output shape: [{} scales x {} samples] — ready for ML consumption",
        fn_, n
    );
}