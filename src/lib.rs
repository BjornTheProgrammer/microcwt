#![no_std]

use alloc::vec;
use alloc::vec::Vec;
use microfft::Complex32;
use num::Complex;

use crate::{scales::Scales, wavelet::Wavelet};
extern crate alloc;

pub mod scales;
pub mod wavelet;

pub struct MicroCwt<'a> {
    pub wavelet: &'a dyn Wavelet,
    threads: i32,
    size: i32,
    fs: f32,
    f0: f32,
    f1: f32,
    afn: f32,
    use_optimalization_schemes: bool,
    use_normalization: bool,
}

impl<'a> MicroCwt<'a> {
    pub fn new(
        pwav: &'a dyn Wavelet,
        pthreads: i32,
        puse_optimalization_schemes: bool,
        puse_normalization: bool,
    ) -> Self {
        Self {
            wavelet: pwav,
            threads: pthreads,
            use_optimalization_schemes: puse_optimalization_schemes,
            use_normalization: puse_normalization,
            size: 0,
            fs: 0.0,
            f0: 0.0,
            f1: 0.0,
            afn: 0.0,
        }
    }

    pub fn cwt(input: &[f32], output: Vec<Complex<f32>>, scales: &mut Scales, complex_input: bool) {
        let nt: usize = input.len().next_power_of_two();
        let newsize = 1 << nt;

        let mut o1 = vec![Complex32::new(0.0, 0.0); newsize];

        match complex_input {
            true => todo!(),
            false => todo!(),
        }
    }
}
