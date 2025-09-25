pub struct QuasiRandom {
    alpha: Vec<f64>,
    counter: u32,
}

impl QuasiRandom {
    pub fn new(dimensions: usize) -> Self {
        let phi = Self::compute_phi(dimensions);
        let mut alpha = Vec::with_capacity(dimensions);

        for i in 0..dimensions {
            alpha.push((1.0 / phi).powi(i as i32 + 1) % 1.0);
        }

        Self { alpha, counter: 0 }
    }

    fn compute_phi(d: usize) -> f64 {
        let mut x = 2.0f64;
        for _ in 0..30 {
            x = (1.0 + x).powf(1.0 / (d as f64 + 1.0));
        }
        x
    }

    pub fn next(&mut self) -> Vec<f64> {
        let offset = 0.5;
        let mut result = Vec::with_capacity(self.alpha.len());

        for &alpha_i in &self.alpha {
            result.push((offset + alpha_i * (self.counter as f64 + 1.0)) % 1.0);
        }

        self.counter += 1;
        result
    }
}