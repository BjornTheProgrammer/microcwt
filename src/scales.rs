use alloc::vec::Vec;
#[allow(unused_imports)]
use num::Float;
use num::traits::Pow;

use crate::wavelet::Wavelet;

pub enum ScaleType {
    FcwtLinscales,
    FcwtLogscales,
    FcwtLinfreqs,
}

pub struct Scales {
    pub scales: Vec<f32>,
    pub fs: i32,
    pub fourwavl: f32,
    pub nscales: i32,
}

impl Scales {
    pub fn new(
        wavelet: &dyn Wavelet,
        scale_type: ScaleType,
        afs: i32,
        af0: f32,
        af1: f32,
        afn: i32,
    ) -> Self {
        let mut scales = Self {
            fs: afs,
            scales: alloc::vec![0f32; afn as usize],
            fourwavl: wavelet.get_four_wavelen(),
            nscales: afn,
        };

        match scale_type {
            ScaleType::FcwtLogscales => scales.calculate_logscale_array(2.0f32, afs, af0, af1, afn),
            ScaleType::FcwtLinscales => scales.calculate_linscale_array(afs, af0, af1, afn),
            ScaleType::FcwtLinfreqs => scales.calculate_linfreq_array(afs, af0, af1, afn),
        }

        scales
    }

    fn calculate_logscale_array(&mut self, base: f32, fs: i32, nf0: f32, nf1: f32, afn: i32) {
        // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;

        let s0: f32 = fs as f32 / nf1;
        let s1: f32 = fs as f32 / nf0;

        //Cannot pass the nyquist frequency
        assert!(
            nf1 <= (fs / 2) as f32,
            "Max frequency cannot be higher than the Nyquist frequency (fs/2)"
        );

        let power0: f32 = s0.log(10.0) / base.log(10.0);
        let power1: f32 = s1.log(10.0) / base.log(10.0);
        let dpower: f32 = power1 - power0;

        for i in 0..afn {
            let power = power0 + (dpower / (afn - 1) as f32) * i as f32;
            self.scales[i as usize] = base.pow(power);
        }
    }

    fn calculate_linscale_array(&mut self, fs: i32, nf0: f32, nf1: f32, afn: i32) {
        //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
        let s0: f32 = fs as f32 / nf1;
        let s1: f32 = fs as f32 / nf0;

        //Cannot pass the nyquist frequency
        assert!(
            nf1 <= (fs / 2) as f32,
            "Max frequency cannot be higher than the Nyquist frequency (fs/2)"
        );
        let ds: f32 = s1 - s0;

        for i in 0..afn {
            self.scales[i as usize] = s0 + (ds / afn as f32) * i as f32;
        }
    }

    fn calculate_linfreq_array(&mut self, fs: i32, nf0: f32, nf1: f32, afn: i32) {
        //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;

        //Cannot pass the nyquist frequency
        assert!(
            nf1 <= (fs / 2) as f32,
            "Max frequency cannot be higher than the Nyquist frequency (fs/2)"
        );
        let df: f32 = nf1 - nf0;

        for i in 0..afn {
            self.scales[(afn - i - 1) as usize] =
                (fs as f32) / (nf0 + (df / afn as f32) * i as f32);
        }
    }

    pub fn get_scales(&self, pfreqs: &mut [f32]) {
        pfreqs.copy_from_slice(&self.scales[..pfreqs.len()]);
    }

    pub fn get_frequencies(&self, pfreqs: &mut [f32]) {
        for (dst, &s) in pfreqs.iter_mut().zip(self.scales.iter()) {
            *dst = self.fs as f32 / s;
        }
    }
}
