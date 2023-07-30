#[derive(Debug, Clone)]
pub struct TwoPolynomial {
    pub two_poly: Vec<f64>,
    pub degree: usize,
}

impl TwoPolynomial {
    #[allow(dead_code)]
    pub fn new(degree_: usize) -> Self {
        TwoPolynomial {
            two_poly: vec![0.0; (degree_ + 1) * (degree_ + 1)],
            degree: degree_,
        }
    }

    #[allow(dead_code)]
    pub fn eval_xy(&self, x: f64, y: f64) -> f64 {
        let n = self.degree;
        let mut t = self.eval_y(n, y);
        for i in 1..(n + 1) {
            t = self.eval_y(n - i, y) + x * t;
        }
        t
    }

    #[allow(dead_code)]
    pub fn eval_y(&self, order: usize, y: f64) -> f64 {
        let n = self.degree;
        let m = order;
        let mut t = self.two_poly[m * (n + 1) + (n - m)];
        for i in ((m * (n + 1))..(m * (n + 1) + (n - m))).rev() {
            t = self.two_poly[i] + y * t;
        }
        t
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn two_poly_degree_1() {
        let test = TwoPolynomial::new(3);
        assert_eq!(test.degree, 3);
    }

    #[test]
    fn two_poly_degree_2() {
        let test = TwoPolynomial::new(4);
        assert_eq!(test.two_poly.len(), 25);
    }

    #[test]
    fn two_poly_eval_1() {
        let mut test = TwoPolynomial::new(3);
        test.two_poly = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0].to_vec();
        assert_eq!(test.eval_y(0, 1.0), 10.0);
        assert_eq!(test.eval_y(1, 1.0), 18.0);
        assert_eq!(test.eval_y(2, 1.0), 19.0);
        assert_eq!(test.eval_y(3, 1.0), 13.0);
    }

    #[test]
    fn two_poly_eval_2() {
        let mut test = TwoPolynomial::new(3);
        test.two_poly = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0].to_vec();
        assert_eq!(test.eval_y(0, 2.0), 49.0);
        assert_eq!(test.eval_y(1, 2.0), 45.0);
        assert_eq!(test.eval_y(2, 2.0), 29.0);
        assert_eq!(test.eval_y(3, 2.0), 13.0);
    }

    #[test]
    fn two_poly_eval_3() {
        let mut test = TwoPolynomial::new(3);
        test.two_poly = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0].to_vec();
        assert_eq!(test.eval_xy(1.0, 1.0), 60.0);
        assert_eq!(test.eval_xy(2.0, 1.0), 226.0);
    }

    #[test]
    fn two_poly_eval_4() {
        let mut test = TwoPolynomial::new(2);
        test.two_poly = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0].to_vec();
        assert_eq!(test.eval_xy(1.0, 1.0), 22.0);
        assert_eq!(test.eval_xy(2.0, 1.0), 52.0);
        assert_eq!(test.eval_xy(2.0, 2.0), 73.0);
        assert_eq!(test.eval_xy(1.0, 2.0), 38.0);
        assert_eq!(test.eval_xy(0.0, 0.0), 1.0);
    }
}