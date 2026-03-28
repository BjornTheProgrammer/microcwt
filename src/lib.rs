#![no_std]

use alloc::vec::Vec;
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

    pub fn cwt(
        input: Vec<f32>,
        size: usize,
        output: Vec<Complex<f32>>,
        scales: &mut Scales,
        complexinput: bool,
    ) {
        let nt: usize = size.next_power_of_two();
    }
}
