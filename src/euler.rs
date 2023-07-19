#[derive(Debug)]
pub struct Coef {
    pub coef: Vec<f64>,
}

impl Coef {
    pub fn new(dim: usize) -> Vec<f64> {
        vec![0.0; dim]
    }

    pub fn euler(self: &mut Coef, vec: &Vec<cgmath::Vector2<f64>>) {
        let dt = 1.0e-4;
        loop {
            let pre = self.err_square(vec);
            //
            for i in 0..self.coef.len() {
                let mut grad = 0.0;
                for j in 0..vec.len() {
                    grad += self.grad_U(i, vec[j].x, vec[j].y);
                }
                self.coef[i] += -dt * grad;
            }
            //
            let post = self.err_square(vec);
            //println!("{:?}, {:?}", pre, post);
            if (pre - post).powf(2.0).sqrt() < 1.0e-12 {
                break;
            }
        }
    }

    fn err_square(self: &Coef, vec: &Vec<cgmath::Vector2<f64>>) -> f64 {
        let mut diff = 0.0;
        for i in 0..self.coef.len() {
            diff = vec[i].y;
            for j in 0..vec.len() {
                diff += self.coef[i] * vec[j].x.powf(j as f64);
            }
            diff -= vec[i].y;
        }
        diff.powf(2.0) / vec.len() as f64
    }

    #[allow(non_snake_case)]
    fn grad_U(self: &Coef, dim: usize, x: f64, y: f64) -> f64 {
        let mut diff = y;
        for i in 0..self.coef.len() {
            diff -= self.coef[i] * x.powf(i as f64);
        }
        -2.0 * x.powf(dim as f64) * diff
    }
}
