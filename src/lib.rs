#![no_std]

use alloc::vec;
use alloc::vec::Vec;
use microfft::Complex32;
use num::Complex;

use crate::{scales::Scales, wavelet::Morlet, wavelet::Wavelet};
extern crate alloc;

pub mod scales;
pub mod wavelet;

pub struct MicroCwt<'a> {
    pub wavelet: &'a mut Morlet,
    _threads: i32,
    _size: i32,
    _fs: f32,
    _f0: f32,
    _f1: f32,
    _afn: f32,
    _use_optimalization_schemes: bool,
    pub use_normalization: bool,
}

impl<'a> MicroCwt<'a> {
    pub fn new(
        pwav: &'a mut Morlet,
        pthreads: i32,
        puse_optimalization_schemes: bool,
        puse_normalization: bool,
    ) -> Self {
        Self {
        wavelet: pwav,
        _threads: pthreads,
        _size: 0,
        _fs: 0.0,
        _f0: 0.0,
        _f1: 0.0,
        _afn: 0.0,
        _use_optimalization_schemes: puse_optimalization_schemes,
        use_normalization: puse_normalization,
        }
    }

    pub fn cwt(
        &mut self,
        input: &[f32],
        output: &mut [Complex<f32>],
        scales: &mut Scales,
        complex_input: bool,
    ) {
        // psize = number of *samples* (not bytes)
        let psize = if complex_input {
            input.len() / 2
        } else {
            input.len()
        };

        // next_power_of_two() on the sample count, NOT on input.len()
        let newsize = psize.next_power_of_two();

        // Ihat: frequency-domain buffer after forward FFT
        // O1:   per-scale scratch buffer for daughter multiplication + IFFT
        // Both zero-initialised to newsize (zero-padding)
        let mut ihat: Vec<Complex32> = vec![Complex32::new(0.0, 0.0); newsize];
        let mut o1:   Vec<Complex32> = vec![Complex32::new(0.0, 0.0); newsize];

        match complex_input {
            true => {
                // input is interleaved [re0, im0, re1, im1, ...]
                // Copy psize complex samples into ihat, rest stays zero-padded
                for i in 0..psize {
                    ihat[i] = Complex32::new(input[2 * i], input[2 * i + 1]);
                }
                fft_inplace(&mut ihat, false);
            }
            false => {
                // Real input: embed into complex buffer with imag = 0
                for i in 0..psize {
                    ihat[i] = Complex32::new(input[i], 0.0);
                }
                fft_inplace(&mut ihat, false);
            }
        }

        for i in 1..(newsize >> 1) {
            ihat[newsize - i] = Complex32::new(ihat[i].re, -ihat[i].im);
        }

        // Fills self.wavelet.mother (Box<[f32]>) via generate_freq()
        self.wavelet.generate_freq(newsize);

        let nscales = scales.nscales as usize;
        for i in 0..nscales {
            let scale   = scales.scales[i];
            let is_last = i == nscales - 1;
            let out_slice = &mut output[i * psize..(i + 1) * psize];

            self.convolve(&ihat, &mut o1, out_slice, scale, psize, newsize, is_last);
        }
    }

    fn convolve(
        &self,
        ihat:    &[Complex32],
        o1:      &mut Vec<Complex32>,
        out:     &mut [Complex<f32>],  // length == psize
        scale:   f32,
        psize:   usize,
        newsize: usize,
        _last:   bool,
    ) {
        // Zero scratch buffer before each scale
        for v in o1.iter_mut() {
            *v = Complex32::new(0.0, 0.0);
        }

        // Daughter wavelet × Ihat → o1
        self.daughter_wavelet_multiplication(ihat, o1, scale, newsize);

        // fft_inplace(inverse=true) applies the 1/N scale internally,
        // matching FFTW_BACKWARD which is unnormalised in C++ — the C++ then
        // calls fft_normalize separately when use_normalization=true.
        fft_inplace(o1, true);

        // Optional extra normalisation pass
        if self.use_normalization {
            let norm = 1.0 / newsize as f32;
            for v in o1.iter_mut() {
                v.re *= norm;
                v.im *= norm;
            }
        }

        // Copy first psize elements from o1 into the caller's output slice
        for (dst, src) in out.iter_mut().zip(o1.iter()) {
            *dst = Complex::new(src.re, src.im);
        }
    }

    // For standard Morlet: imag_frequency=false, doublesided=false, so only
    // the simple forward multiply executes — matching main.cpp exactly.
    fn daughter_wavelet_multiplication(
        &self,
        input:  &[Complex32],
        output: &mut [Complex32],
        scale:  f32,
        isize:  usize,
    ) {
        // Access mother, imag_frequency, doublesided directly as pub fields
        // from your wavelet.rs Morlet struct
        let mother      = &self.wavelet.mother;
        let imaginary   = self.wavelet.imag_frequency;
        let doublesided = self.wavelet.doublesided;

        let isizef    = isize as f32;
        let endpointf = (isizef / 2.0_f32).min((isizef * 2.0_f32) / scale);
        let step      = scale / 2.0_f32;
        let endpoint  = endpointf as usize;
        let maximum   = isizef - 1.0_f32;

        // Forward (positive-frequency) half
        for q in 0..endpoint {
            let tmp = (step * q as f32).min(maximum) as usize;
            let m   = mother[tmp];

            // (1 - 2*imaginary): 1.0 when false, -1.0 when true
            // When imaginary=true the AVX path also shuffles re/im
            // (shuffle mask 177 swaps adjacent pairs), so we replicate that:
            //   normal:    re *= m,       im *= m * sign
            //   imaginary: re = im*m,     im = -re*m   (swap + negate)
            output[q] = if imaginary {
                Complex32::new( input[q].im * m,
                               -input[q].re * m)
            } else {
                Complex32::new(input[q].re * m,
                               input[q].im * m)
            };
        }

        // Doublesided (negative-frequency mirror)
        // Not active for standard Morlet but translated faithfully
        if doublesided {
            let s1 = isize - 1;
            for q in 0..endpoint {
                let tmp = (step * q as f32).min(maximum) as usize;
                let m   = mother[tmp];
                let idx = s1 - q;

                // Sign flip is on re here (swapped vs forward half) —
                // matches C++ non-AVX doublesided path exactly
                output[idx] = if imaginary {
                    Complex32::new( input[idx].im * m,
                                   -input[idx].re * m)
                } else {
                    Complex32::new(-input[idx].re * m,
                                   input[idx].im * m)
                };
            }
        }
    }
}

// microfft requires const-generic fixed-size arrays so we dispatch on runtime
// length. Extend the match arms for longer signals if needed.
//
// IFFT identity used: IFFT(x) = conj(FFT(conj(x))) / N
// This avoids needing a separate inverse plan.
fn fft_inplace(buf: &mut [Complex32], inverse: bool) {
    // conjugate input for inverse transform
    if inverse {
        for v in buf.iter_mut() {
            v.im = -v.im;
        }
    }

    // dispatch to fixed-size microfft implementation
    match buf.len() {
        512   => { let a: &mut [_; 512]   = buf.try_into().unwrap(); microfft::complex::cfft_512(a); }
        1024  => { let a: &mut [_; 1024]  = buf.try_into().unwrap(); microfft::complex::cfft_1024(a); }
        2048  => { let a: &mut [_; 2048]  = buf.try_into().unwrap(); microfft::complex::cfft_2048(a); }
        4096  => { let a: &mut [_; 4096]  = buf.try_into().unwrap(); microfft::complex::cfft_4096(a); }
        8192  => { let a: &mut [_; 8192]  = buf.try_into().unwrap(); microfft::complex::cfft_8192(a); }
        16384 => { let a: &mut [_; 16384] = buf.try_into().unwrap(); microfft::complex::cfft_16384(a); }
        _ => panic!("Unsupported FFT size: {}. Add arm to fft_inplace.", buf.len()),
    }

    // conjugate output and scale by 1/N for inverse
    if inverse {
        let norm = 1.0 / buf.len() as f32;
        for v in buf.iter_mut() {
            v.im = -v.im;
            v.re *= norm;
            v.im *= norm;
        }
    }
}