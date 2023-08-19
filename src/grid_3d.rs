use super::point;
use super::two_variable_polynomial;

#[derive(Debug)]
pub struct Grid3D {
    pub points_3d: Vec<point::Point3>,
}

impl Grid3D {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Grid3D { points_3d: vec![] }
    }

    #[allow(dead_code)]
    pub fn push(&mut self, vec: point::Point3) {
        self.points_3d.push(point::Point3 {
            x: vec.x,
            y: vec.y,
            z: vec.z,
        });
    }

    #[allow(dead_code)]
    pub fn poly_fitting_by_euler_with_tol(
        &mut self,
        poly: &mut two_variable_polynomial::TwoPolynomial,
        tol: f64,
    ) -> Vec<f64> {
        let mut dt = 1.0e-3;
        loop {
            let pre: f64 = self.potential(&poly);
            let mut tmp = poly.clone();
            let tmp_a: two_variable_polynomial::TwoPolynomial = self.euler_step(&mut tmp, dt);
            let post: f64 = self.potential(&tmp_a);
            // 以下経験的なパラメータあり
            //println!("dt: {:?}, pre: {:?}, post: {:?}, coef: {:?}", dt, pre, post, poly);
            if pre > post {
                dt = (1.01 * dt).min(1.0e3);
                *poly = tmp;
            } else if post > pre {
                dt = (0.9 * dt).max(1.0e-5);
            }
            if post < tol {
                break;
            }
        }
        poly.two_poly.to_vec()
    }

    #[allow(dead_code)]
    pub fn euler_step(
        &mut self,
        poly: &mut two_variable_polynomial::TwoPolynomial,
        dt: f64,
    ) -> two_variable_polynomial::TwoPolynomial {
        let num_coef: usize = (poly.degree + 1) * (poly.degree + 1);
        for i in 0..num_coef {
            poly.two_poly[i] = poly.two_poly[i] - dt * self.potential_deriv(&poly)[i];
        }
        poly.clone()
    }

    #[allow(dead_code)]
    pub fn potential_deriv(&mut self, poly: &two_variable_polynomial::TwoPolynomial) -> Vec<f64> {
        let num_coef: usize = (poly.degree + 1) * (poly.degree + 1);
        let mut du = vec![0.0; num_coef];
        let n = poly.degree;
        for j in 0..self.points_3d.len() {
            for m in 0..(n + 1) {
                for i in ((m * (n + 1))..(m * (n + 1) + (n - m) + 1)).rev() {
                    //println!("{:?}, {:?}, {:?}, {:?}", i, n + 1, i / (n + 1), i % (n + 1));
                    let x_deg = i / (n + 1);
                    let y_deg = i % (n + 1);
                    du[i] -= self.points_3d[j].x.powf(x_deg as f64)
                        * self.points_3d[j].y.powf(y_deg as f64)
                        * (self.points_3d[j].z
                            - poly.eval_xy(self.points_3d[j].x, self.points_3d[j].y));
                }
            }
        }
        du
    }

    #[allow(dead_code)]
    pub fn potential(&mut self, poly: &two_variable_polynomial::TwoPolynomial) -> f64 {
        let mut u = 0.0;
        for i in 0..self.points_3d.len() {
            u += (self.points_3d[i].z
                - self.poly_eval(poly, self.points_3d[i].x, self.points_3d[i].y))
            .powf(2.0);
        }
        u / self.points_3d.len() as f64
    }

    #[allow(dead_code)]
    fn poly_eval(
        &mut self,
        two_poly: &two_variable_polynomial::TwoPolynomial,
        x: f64,
        y: f64,
    ) -> f64 {
        two_poly.eval_xy(x, y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_poly_deriv_1() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 1.0, 0.0, 0.0].to_vec(); // (x - 1)^2
        let mut data = Grid3D::new();
        data.points_3d.push(point::Point3 {
            x: 1.0,
            y: -100.0,
            z: 0.0,
        });
        data.points_3d.push(point::Point3 {
            x: 0.0,
            y: 100.0,
            z: 1.0,
        });
        data.points_3d.push(point::Point3 {
            x: 3.0,
            y: 100.0,
            z: 4.0,
        });
        assert_eq!(data.potential(&poly), 0.0);
    }

    #[test]
    fn two_poly_deriv_2() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 1.0, 0.0, 0.0].to_vec(); // (x - 1)^2
        let mut data = Grid3D::new();
        data.points_3d.push(point::Point3 {
            x: 1.0,
            y: -100.0,
            z: 0.0,
        });
        data.points_3d.push(point::Point3 {
            x: 0.0,
            y: 100.0,
            z: 1.0,
        });
        data.points_3d.push(point::Point3 {
            x: 3.0,
            y: 100.0,
            z: 4.0,
        });
        assert_eq!(data.potential_deriv(&poly)[0], 0.0);
        assert_eq!(data.potential_deriv(&poly)[1], 0.0);
        assert_eq!(data.potential_deriv(&poly)[2], 0.0);
        assert_eq!(data.potential_deriv(&poly)[3], 0.0);
        assert_eq!(data.potential_deriv(&poly)[4], 0.0);
        assert_eq!(data.potential_deriv(&poly)[5], 0.0);
    }

    #[test]
    fn potential_deriv_1() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [1.0, 0.0, 0.0, -2.0, 0.0, 0.0, 1.0, 0.0, 0.0].to_vec(); // (x - 1)^2
        let mut test = Grid3D::new();
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        });
        test.points_3d.push(point::Point3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        });
        test.points_3d.push(point::Point3 {
            x: 2.0,
            y: 0.0,
            z: 1.0,
        });
        assert_eq!(test.potential_deriv(&poly)[0], 0.0);
        assert_eq!(test.potential_deriv(&poly)[1], 0.0);
        assert_eq!(test.potential_deriv(&poly)[2], 0.0);
        assert_eq!(test.potential_deriv(&poly)[3], 0.0);
        assert_eq!(test.potential_deriv(&poly)[4], 0.0);
        assert_eq!(test.potential_deriv(&poly)[5], 0.0);
    }

    #[test]
    fn potential_deriv_2() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [1.0, -2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].to_vec(); // (y - 1)^2
        let mut test = Grid3D::new();
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        });
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        });
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 2.0,
            z: 1.0,
        });
        assert_eq!(test.potential_deriv(&poly)[0], 0.0);
        assert_eq!(test.potential_deriv(&poly)[1], 0.0);
        assert_eq!(test.potential_deriv(&poly)[2], 0.0);
        assert_eq!(test.potential_deriv(&poly)[3], 0.0);
        assert_eq!(test.potential_deriv(&poly)[4], 0.0);
        assert_eq!(test.potential_deriv(&poly)[5], 0.0);
    }

    #[test]
    fn potential_deriv_3() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].to_vec();
        let mut test = Grid3D::new();
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        });
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        });
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 2.0,
            z: 1.0,
        });
        assert_eq!(test.potential_deriv(&poly)[0], -2.0);
        assert_eq!(test.potential_deriv(&poly)[1], -2.0);
        assert_eq!(test.potential_deriv(&poly)[2], -4.0);
        assert_eq!(test.potential_deriv(&poly)[3], 0.0);
        assert_eq!(test.potential_deriv(&poly)[4], 0.0);
        assert_eq!(test.potential_deriv(&poly)[5], 0.0);
    }

    #[test]
    fn potential_deriv_4() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].to_vec();
        let mut test = Grid3D::new();
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        test.points_3d.push(point::Point3 {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        });
        test.points_3d.push(point::Point3 {
            x: 2.0,
            y: 0.0,
            z: 4.0,
        });
        assert_eq!(test.potential_deriv(&poly)[0], -5.0);
        assert_eq!(test.potential_deriv(&poly)[1], -1.0);
        assert_eq!(test.potential_deriv(&poly)[2], -1.0);
        assert_eq!(test.potential_deriv(&poly)[3], -9.0);
        assert_eq!(test.potential_deriv(&poly)[4], -1.0);
        assert_eq!(test.potential_deriv(&poly)[5], 0.0);
    }

    #[test]
    fn euler_step() {
        let mut poly = two_variable_polynomial::TwoPolynomial::new(2);
        poly.two_poly = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0].to_vec();
        let mut test = Grid3D::new();
        test.points_3d.push(point::Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        test.points_3d.push(point::Point3 {
            x: 1.0,
            y: 0.0,
            z: 1.0,
        });
        test.points_3d.push(point::Point3 {
            x: 2.0,
            y: 0.0,
            z: 4.0,
        });
        test.points_3d.push(point::Point3 {
            x: -1.0,
            y: 0.0,
            z: 1.0,
        });
        test.points_3d.push(point::Point3 {
            x: -2.0,
            y: 0.0,
            z: 4.0,
        });
        let tol = 1.0e-1;
        let coef = test.poly_fitting_by_euler_with_tol(&mut poly, tol);
        assert_eq!(coef[0], 0.2257267851632633);
        assert_eq!(coef[1], 0.0);
        assert_eq!(coef[2], 0.0);
        assert_eq!(coef[3], 0.0);
        assert_eq!(coef[4], 0.0);
        assert_eq!(coef[5], 0.0);
    }
}
