#[derive(Debug, Clone)]
pub struct Coef {
    pub coef: Vec<f64>,
    pub dt: f64,
}

impl Coef {
    pub fn new(dim: usize) -> Self {
        Coef {
            coef: vec![0.0; dim],
            dt: 1.0e-3,
        }
    }

    #[allow(dead_code)]
    pub fn euler(self: &mut Coef, vec: &Vec<cgmath::Vector2<f64>>) {
        loop {
            let pre = self.err_abs(vec);
            //let tmp = self.coef.clone();
            //
            for i in 0..self.coef.len() {
                let mut grad = 0.0;
                for j in 0..vec.len() {
                    grad += self.grad_U(i, vec[j].x, vec[j].y);
                }
                self.coef[i] += -self.dt * grad;
            }
            //
            //println!("{:?}", self.coef);
            let post = self.err_abs(vec);
            //println!("{:?}, {:?}, {:?}, {:?}", pre, post, pre - post, self.dt);
            //if post > pre {
            //println!("{:?}, {:?}, {:?}, {:?}", pre, post, pre - post, self.dt);
            //self.coef = tmp.clone();
            //self.dt = 0.9999 * self.dt;
            //} else {
            //self.dt = 1.0001 * self.dt;
            //}
            if (pre > 0.0) && (post > 0.0) && (pre - post).abs() < 1.0e-10 {
                break;
            }
        }
    }

    #[allow(dead_code)]
    fn err_abs(self: &Coef, vec: &Vec<cgmath::Vector2<f64>>) -> f64 {
        let mut diff = 0.0;
        for i in 0..self.coef.len() {
            for j in 0..vec.len() {
                diff += self.coef[i] * vec[j].x.powf(i as f64);
            }
            diff -= vec[i].y;
        }
        diff.abs() / vec.len() as f64
    }

    #[allow(non_snake_case)]
    #[allow(dead_code)]
    fn grad_U(self: &Coef, dim: usize, x: f64, y: f64) -> f64 {
        let mut diff = y;
        for i in 0..self.coef.len() {
            diff -= self.coef[i] * x.powf(i as f64);
        }
        -2.0 * x.powf(dim as f64) * diff
    }

    pub fn eval(self: &mut Coef, x: f64) -> f64 {
        let mut ret = 0.0;
        for i in 0..self.coef.len() {
            ret += self.coef[i] * x.powf(i as f64);
        }
        ret
    }
}
