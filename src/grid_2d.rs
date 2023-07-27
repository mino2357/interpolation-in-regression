//use super::embedded_runge_kutta;
use super::vector2;

#[derive(Debug)]
pub struct Grid2D {
    pub point_2d: Vec<vector2::Vector2>,
}

impl Grid2D {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Grid2D {
            point_2d: vec![],
        }
    }

    #[allow(dead_code)]
    pub fn poly_fitting_by_euler_with_tol(&mut self, poly: &mut Vec<f64>, tol: f64) -> Vec<f64> {
        let mut dt = 1.0e-3;
        loop {
            let pre: f64 = self.potential(&poly);
            let tmp: Vec<f64> = self.euler_step(poly, dt);
            let post: f64 = self.potential(&tmp);
            // 以下経験的なパラメータあり
            //println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
            if pre > post {
                dt = (1.01 * dt).min(1.0e0);
                *poly = tmp;
            } else if post > pre {
                dt = (0.9 * dt).max(1.0e-5);
            }
            if post < tol || (pre - post).abs() < 1.0e-12 {
                //println!("dt: {:?}, pre: {:?}, post: {:?}", dt, pre, post);
                break;
            }
        }
        poly.to_vec()
    }

    #[allow(dead_code)]
    pub fn poly_fitting_by_classical_rk4_with_tol(&mut self, poly: &mut Vec<f64>, tol: f64) -> Vec<f64> {
        let mut dt = 1.0e-4;
        loop {
            let pre: f64 = self.potential(&poly);
            let tmp: Vec<f64> = self.classical_rk4_step(poly, dt);
            let post: f64 = self.potential(&tmp);
            //println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
            if pre > post {
                dt = (1.01 * dt).min(1.0e0);
                *poly = tmp;
            } else if post > pre {
                dt = (0.9 * dt).max(1.0e-5);
            }
            if post < tol || (pre - post).abs() < 1.0e-12 {
                //println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
                break;
            }
        }
        poly.to_vec()
    }

    #[allow(dead_code)]
    pub fn poly_fitting_by_classical_rk4(&mut self, poly: &mut Vec<f64>) -> Vec<f64> {
        let mut dt = 1.0e-4;
        loop {
            let pre: f64 = self.potential(&poly);
            let tmp: Vec<f64> = self.classical_rk4_step(poly, dt);
            let post: f64 = self.potential(&tmp);
            println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
            if pre > post {
                dt = (1.01 * dt).min(1.0e0);
                *poly = tmp;
            } else if post > pre {
                dt = (0.9 * dt).max(1.0e-5);
            }
            if post < 1.0e-3 || (pre - post).abs() < 1.0e-12 {
                println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
                break;
            }
        }
        poly.to_vec()
    }

    #[allow(dead_code)]
    pub fn poly_fitting_by_modified_euler(&mut self, poly: &mut Vec<f64>) -> Vec<f64> {
        let mut dt = 1.0e-3;
        loop {
            let pre: f64 = self.potential(&poly);
            let tmp: Vec<f64> = self.modified_euler_step(poly, dt);
            let post: f64 = self.potential(&tmp);
            println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
            if pre > post {
                dt = (1.01 * dt).min(1.0e0);
                *poly = tmp;
            } else if post > pre {
                dt = (0.9 * dt).max(1.0e-5);
            }
            if post < 1.0e-3 || (pre - post).abs() < 1.0e-12 {
                println!("dt: {:?}, pre: {:?}, post: {:?}", dt, pre, post);
                break;
            }
        }
        poly.to_vec()
    }

    #[allow(dead_code)]
    pub fn poly_fitting_by_euler(&mut self, poly: &mut Vec<f64>) -> Vec<f64> {
        let mut dt = 1.0e-3;
        loop {
            let pre: f64 = self.potential(&poly);
            let tmp: Vec<f64> = self.euler_step(poly, dt);
            let post: f64 = self.potential(&tmp);
            // 以下経験的なパラメータあり
            println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
            if pre > post {
                dt = (1.01 * dt).min(1.0e0);
                *poly = tmp;
            } else if post > pre {
                dt = (0.9 * dt).max(1.0e-5);
            }
            if post < 1.0e-3 || (pre - post).abs() < 1.0e-12 {
                println!("dt: {:?}, pre: {:?}, post: {:?}", dt, pre, post);
                break;
            }
        }
        poly.to_vec()
    }

    #[allow(dead_code)]
    pub fn classical_rk4_step(&mut self, poly: &mut Vec<f64>, dt: f64) -> Vec<f64> { // Heun's method
        let dim = poly.len();
        // k1
        let mut vec1 = poly.clone();
        for i in 0..dim {
            vec1[i] = poly[i];
        }
        let mut k1 = vec![0.0; dim];
        for i in 0..dim {
            k1[i] = self.potential_deriv(&vec1)[i];
        };
        // k2
        let mut vec2 = vec![0.0; dim];
        for i in 0..dim {
            vec2[i] = poly[i] + dt * 0.5 * k1[i];
        }
        let mut k2 = vec![0.0; dim];
        for i in 0..dim {
            k2[i] = self.potential_deriv(&vec2)[i];
        }
        // k3
        let mut vec3 = vec![0.0; dim];
        for i in 0..dim {
            vec3[i] = poly[i] + dt * 0.5 * k2[i];
        }
        let mut k3 = vec![0.0; dim];
        for i in 0..dim {
            k3[i] = self.potential_deriv(&vec3)[i];
        }
        // k4
        let mut vec4 = vec![0.0; dim];
        for i in 0..dim {
            vec4[i] = poly[i] + dt * k3[i];
        }
        let mut k4 = vec![0.0; dim];
        for i in 0..dim {
            k4[i] = self.potential_deriv(&vec4)[i];
        }
        // sum
        for i in 0..dim {
            poly[i] = poly[i] - 1.0 / 6.0 * dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        poly.clone()
    }

    #[allow(dead_code)]
    pub fn modified_euler_step(&mut self, poly: &mut Vec<f64>, dt: f64) -> Vec<f64> { // Heun's method
        let mut tmp = vec![0.0; poly.len()];
        for i in 0..poly.len() {
            tmp[i] = poly[i] - dt * self.potential_deriv(&poly)[i];
        }
        for i in 0..poly.len() {
            poly[i] = poly[i] - 0.5 * dt * (self.potential_deriv(&poly)[i] + self.potential_deriv(&tmp)[i]);
        }
        poly.clone()
    }

    #[allow(dead_code)]
    pub fn euler_step(&mut self, poly: &mut Vec<f64>, dt: f64) -> Vec<f64> {
        for i in 0..poly.len() {
            poly[i] = poly[i] - dt * self.potential_deriv(&poly)[i];
        }
        poly.clone()
    }

    #[allow(dead_code)]
    pub fn potential_deriv(&mut self, poly: &Vec<f64>) -> Vec<f64> {
        let mut du = vec![0.0; poly.len()];
        for i in 0..poly.len() {
            for j in 0..self.point_2d.len() {
                du[i] -= self.point_2d[j].x.powf(i as f64) * (self.point_2d[j].y - self.eval(poly, self.point_2d[j].x));
            }
        }
        du
    }

    #[allow(dead_code)]
    pub fn potential(&mut self, poly: &Vec<f64>) -> f64 {
        let mut u = 0.0;
        for i in 0..self.point_2d.len() {
            u += (self.point_2d[i].y - self.eval(poly, self.point_2d[i].x)).powf(2.0);
        }
        u
    }

    // FMA. ref: https://zenn.dev/herumi/articles/poly-evaluation-by-fma
    #[allow(dead_code)]
    pub fn eval(&mut self, poly: &Vec<f64>, x: f64) -> f64 {
        let degree = poly.len() - 1;
        let mut t = poly[degree];
        for i in 1..poly.len() {
            t = poly[degree - i] + x * t;
        }
        t
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn polynomial_eval_test() {
        let poly = vec![1.0, 2.0, 3.0];
        let mut test = Grid2D::new();
        assert_eq!(test.eval(&poly, 1.0), 6.0);
        assert_eq!(test.eval(&poly, 2.0), 17.0);
    }

    #[test]
    fn potential_eval_1() {
        let poly = vec![3.0, 2.0, 1.0];
        let mut test = Grid2D::new();
        test.point_2d.push(vector2::Vector2 {
            x: 0.0,
            y: 3.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: 1.0,
            y: 6.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: 2.0,
            y: 11.0,
        });
        assert_eq!(test.potential(&poly), 0.0);
    }

    #[test]
    fn potential_eval_2() {
        let poly = vec![3.0, 2.0, 1.0];
        let mut test = Grid2D::new();
        test.point_2d.push(vector2::Vector2 {
            x: 0.0,
            y: 4.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: 1.0,
            y: 7.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: 2.0,
            y: 12.0,
        });
        assert_eq!(test.potential(&poly), 3.0);
    }

    #[test]
    fn potential_deriv_1() {
        let poly = vec![3.0, 2.0, 1.0];
        let mut test = Grid2D::new();
        test.point_2d.push(vector2::Vector2 {
            x: 0.0,
            y: 3.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: 1.0,
            y: 6.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: 2.0,
            y: 11.0,
        });
        assert_eq!(test.potential_deriv(&poly)[0], 0.0);
        assert_eq!(test.potential_deriv(&poly)[1], 0.0);
        assert_eq!(test.potential_deriv(&poly)[2], 0.0);
    }

    #[test]
    fn potential_deriv_2() {
        let poly = vec![0.0, 0.0];
        let mut test = Grid2D::new();
        test.point_2d.push(vector2::Vector2 {
            x: 1.0,
            y: 1.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: -1.0,
            y: -1.0,
        });
        assert_eq!(test.potential_deriv(&poly)[0], 0.0);
        assert_eq!(test.potential_deriv(&poly)[1], -2.0);
    }

    #[test]
    fn poly_fitting_euler_1_step() {
        let mut poly = vec![0.0, 0.0];
        let mut test = Grid2D::new();
        test.point_2d.push(vector2::Vector2 {
            x: 1.0,
            y: 1.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: -1.0,
            y: -1.0,
        });
        let dt = 1.0e-3;
        assert_eq!(test.euler_step(&mut poly, dt)[0], 0.0);
        assert_eq!(test.euler_step(&mut poly, dt)[1], 0.0039959999999999996);
    }

    #[test]
    fn poly_fitting_euler_1() {
        let mut poly = vec![0.0, 0.0];
        let mut test = Grid2D::new();
        test.point_2d.push(vector2::Vector2 {
            x: 1.0,
            y: 1.0,
        });
        test.point_2d.push(vector2::Vector2 {
            x: -1.0,
            y: -1.0,
        });
        let a_0 = test.poly_fitting_by_euler(&mut poly)[0];
        let a_1 = test.poly_fitting_by_euler(&mut poly)[1];
        assert_eq!(a_0, 0.0);
        assert_eq!(a_1, 0.978498617059409);
    }
}
