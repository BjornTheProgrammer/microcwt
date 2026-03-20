use alloc::{boxed::Box, vec};
use core::f32::consts::PI;
#[allow(unused_imports)]
use num::{Complex, Float as _};

pub trait Wavelet {
    fn generate_time(&self, real: &mut [f32], imag: &mut [f32], size: i32, scale: f32);
    fn generate_freq(&mut self, size: usize);
    fn get_support(&self, scale: f32) -> i32;
    fn get_wavelet(&self, scale: f32, pwav: &mut [Complex<f32>]);
    fn get_four_wavelen(&self) -> f32;
}

pub struct Morlet {
    pub four_wavelen: f32,
    pub fb: f32,
    fb2: f32,
    ifb: f32,
    pub imag_frequency: bool,
    pub doublesided: bool,
    pub mother: Box<[f32]>,
}

impl Morlet {
    pub fn new(bandwidth: f32) -> Self {
        let fb = bandwidth;
        let fb2 = 2.0f32 * fb * fb;
        let ifb = 1.0f32 / fb;

        Self {
            four_wavelen: 0.9876,
            fb,
            fb2,
            ifb,
            imag_frequency: false,
            doublesided: false,
            mother: Box::new([]),
        }
    }
}

// This is just PI.powf(-0.25)
pub const IPI4: f32 = 0.75112554446;

impl Wavelet for Morlet {
    fn generate_time(&self, real: &mut [f32], imag: &mut [f32], size: i32, scale: f32) {
        // Time domain because we know size from scale
        let width = self.get_support(scale);
        let norm = size as f32 * self.ifb * IPI4;

        for t in 0..(width * 2 + 1) {
            let tmp1 = (t - width) as f32 / scale;
            let tmp2 = (-(tmp1 * tmp1) / self.fb2).exp();

            real[t as usize] = norm * tmp2 * (tmp1 * 2.0f32 * PI).cos() / scale;
            imag[t as usize] = norm * tmp2 * (tmp1 * 2.0f32 * PI).sin() / scale;
        }
    }

    fn generate_freq(&mut self, size: usize) {
        // Frequency domain, because we only need size. Default scale is always 2;
        let width = size;

        let toradians: f32 = (2.0 * PI) / (size as f32);
        let norm: f32 = f32::sqrt(2.0 * PI) * IPI4;

        self.mother = (0..width)
            .map(|w| {
                let mut tmp1 = 2.0f32 * ((w as f32) * toradians) * self.fb - 2.0f32 * PI * self.fb;
                tmp1 = -(tmp1 * tmp1) / 2.0;
                norm * tmp1.exp()
            })
            .collect();
    }

    fn get_support(&self, scale: f32) -> i32 {
        return (self.fb * scale * 3.0f32) as i32;
    }

    fn get_wavelet(&self, scale: f32, pwav: &mut [Complex<f32>]) {
        let w = self.get_support(scale) as usize;
        let pn = pwav.len();
        let size = w * 2 + 1;
        let alloc_size = size.max(pn);

        let mut real = vec![0.0f32; alloc_size];
        let mut imag = vec![0.0f32; alloc_size];

        self.generate_time(&mut real, &mut imag, pn as i32, scale);

        for t in 0..pn {
            pwav[t] = Complex::new(real[t], imag[t]);
        }
    }

    fn get_four_wavelen(&self) -> f32 {
        self.four_wavelen
    }
}
