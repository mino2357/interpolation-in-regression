use crate::grid_3d::Grid3D;
use crate::kd_tree;
use crate::kd_tree::Grid2D;
use crate::point;
use crate::two_variable_polynomial;

#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct WaveEq {
    pub interior: kd_tree::Points2D,
    pub boundary: kd_tree::Points2D, // Circle
    pub value: Vec<f64>,
    pub near_points_boundary: Vec<Vec<usize>>,
    pub near_points_interior: Vec<Vec<usize>>,
    pub near: Vec<Vec<usize>>,
    pub poly: Vec<two_variable_polynomial::TwoPolynomial>,
}

impl WaveEq {
    #[allow(dead_code)]
    pub fn new() -> Self {
        WaveEq {
            interior: kd_tree::Points2D {
                points: Vec::<kd_tree::Grid2D>::new(),
            },
            boundary: kd_tree::Points2D {
                points: Vec::<kd_tree::Grid2D>::new(),
            },
            value: vec![0.0; 0],
            near_points_boundary: vec![vec![0 as usize, 0]; 0],
            near_points_interior: vec![vec![0 as usize, 0]; 0],
            near: vec![vec![0 as usize, 0]; 0],
            poly: vec![two_variable_polynomial::TwoPolynomial::new(2); 0],
        }
    }

    #[allow(dead_code)]
    pub fn make(num_points: usize, dim: usize, tol: f64) -> Self {
        let mut wave = WaveEq::new();
        wave.create(num_points);
        wave.set_initial_condition();
        wave.set_boundary_near_points();
        wave.set_interior_near_points();
        wave.set_interior_near_points_plus_boundary();
        wave.init_poly(dim);
        wave.set_init_poly(tol);
        wave.clone()
    }

    #[allow(dead_code)]
    pub fn poly_eval(&mut self, x: f64, y: f64) -> f64 {
        let vec = Grid2D::new(x, y);
        let tree = kd_tree::KDTree::construct_kd_tree(&self.interior);
        let dx = 0.01 * std::f64::consts::PI / self.boundary.points.len() as f64;
        let mut radius = std::f64::consts::PI / self.boundary.points.len() as f64;
        let mut _index: usize = 0;
        loop {
            let near = tree.neighbor_search(&vec, radius);
            if near.len() > 0 {
                _index = near[0];
                break;
            }
            radius += dx;
        }
        self.poly[_index].eval_xy(x, y)
    }

    #[allow(dead_code)]
    pub fn set_init_poly(&mut self, tol: f64) {
        for i in 0..self.interior.points.len() {
            let mut neighbor_vec = Grid3D::new();
            for j in &self.near_points_interior[i] {
                //if i != *j {
                let v_x = self.interior.points[*j].x;
                let v_y = self.interior.points[*j].y;
                let v_z = self.value[*j];
                let vec = point::Point3 {
                    x: v_x,
                    y: v_y,
                    z: v_z,
                };
                neighbor_vec.push(vec);
                //}
                for k in &self.near_points_boundary[i] {
                    let v_x = self.boundary.points[*k].x;
                    let v_y = self.boundary.points[*k].y;
                    let v_z = 0.0;
                    let vec = point::Point3 {
                        x: v_x,
                        y: v_y,
                        z: v_z,
                    };
                    neighbor_vec.push(vec);
                }
            }
            println!("{} / {}", i, self.interior.points.len());
            neighbor_vec.poly_fitting_by_euler_with_tol(&mut self.poly[i], tol);
        }
    }

    #[allow(dead_code)]
    pub fn init_poly(&mut self, dim: usize) {
        for _ in 0..self.interior.points.len() {
            self.poly
                .push(two_variable_polynomial::TwoPolynomial::new(dim));
        }
    }

    #[allow(dead_code)]
    pub fn set_interior_near_points_plus_boundary(&mut self) {
        for i in 0..self.near_points_interior.len() {
            self.near_points_boundary.push(vec![0 as usize; 0]);
            for k in 0..self.near.len() {
                if self.near[k][0] == i {
                    self.near_points_boundary[i].push(k);
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn set_interior_near_points(&mut self) {
        let tree = kd_tree::KDTree::construct_kd_tree(&self.interior);
        let dx = 0.01 * std::f64::consts::PI / self.boundary.points.len() as f64;
        for i in 0..self.interior.points.len() {
            let mut radius = std::f64::consts::PI / self.boundary.points.len() as f64;
            loop {
                let near = tree.neighbor_search(&self.interior.points[i], radius);
                if near.len() > 3 {
                    self.near_points_interior.push(near);
                    break;
                }
                radius += dx;
            }
        }
    }

    #[allow(dead_code)]
    pub fn set_boundary_near_points(&mut self) {
        let tree = kd_tree::KDTree::construct_kd_tree(&self.interior);
        let dx = 0.01 * std::f64::consts::PI / self.boundary.points.len() as f64;
        for i in 0..self.boundary.points.len() {
            let mut radius = std::f64::consts::PI / self.boundary.points.len() as f64;
            loop {
                let near = tree.neighbor_search(&self.boundary.points[i], radius);
                if near.len() > 0 {
                    self.near.push(near);
                    break;
                }
                radius += dx;
            }
        }
    }

    #[allow(dead_code)]
    pub fn create(&mut self, num_points: usize) {
        for i in 0..num_points {
            let x = (2.0 * std::f64::consts::PI * i as f64 / num_points as f64).cos();
            let y = (2.0 * std::f64::consts::PI * i as f64 / num_points as f64).sin();
            self.boundary.push(x, y);
        }
        let num_interior = (((num_points as f64).powi(2)) / (4.0 * std::f64::consts::PI)) as usize;
        use rand::prelude::*;
        let seed: [u8; 32] = [1; 32];
        let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);
        let mut iter = 0;
        loop {
            let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
            let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
            if x_r * x_r + y_r * y_r < 1.0 {
                self.interior.push(x_r, y_r);
                iter += 1;
            }
            if iter == num_interior {
                break;
            }
        }
        for _ in 0..1000 {
            self.euler_step();
        }
    }

    #[allow(dead_code)]
    pub fn set_initial_condition(&mut self) {
        for i in 0..self.interior.points.len() {
            let x = self.interior.points[i].x;
            let y = self.interior.points[i].y;
            self.value.push((-10.0 * (x * x + y * y)).exp());
        }
    }

    #[allow(dead_code)]
    pub fn euler_step(&mut self) {
        let dt = 1.0e-3;
        for i in 0..self.interior.points.len() {
            let mut x = self.interior.points[i].x - dt * self.lennard_jones_potential_deriv(i).x;
            let mut y = self.interior.points[i].y - dt * self.lennard_jones_potential_deriv(i).y;
            if x * x + y * y < 1.0 {
                self.interior.points[i].x = x;
                self.interior.points[i].y = y;
            } else {
                loop {
                    let tmp = x;
                    x = -0.5 * y;
                    y = 0.5 * tmp;
                    if x * x + y * y < 1.0 {
                        self.interior.points[i].x = x;
                        self.interior.points[i].y = y;
                        break;
                    }
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn lennard_jones_potential_deriv(&mut self, index: usize) -> kd_tree::Grid2D {
        let epsilon: f64 = 1.0;
        let sigma: f64 =
            0.9 * (2.0_f64).powf(-1.0 / 6.0) / (self.interior.points.len() as f64).sqrt();
        let mut f_x = 0.0;
        let mut f_y = 0.0;
        let sigma2 = sigma * sigma;
        let sigma4 = sigma2 * sigma2;
        let sigma8 = sigma4 * sigma4;
        let sigma6 = sigma2 * sigma4;
        let sigma12 = sigma8 * sigma4;
        for i in 0..self.interior.points.len() {
            if i != index {
                let dx = self.interior.points[index].x - self.interior.points[i].x;
                let dy = self.interior.points[index].y - self.interior.points[i].y;
                let r = (dx * dx + dy * dy).sqrt();
                let r2 = r * r;
                let r4 = r2 * r2;
                let r8 = r4 * r4;
                let r13 = r8 * r4 * r;
                let r7 = r2 * r4 * r;
                f_x += 4.0 * epsilon * (12.0 * sigma12 / r13 - 6.0 * sigma6 / r7) * dx;
                f_y += 4.0 * epsilon * (12.0 * sigma12 / r13 - 6.0 * sigma6 / r7) * dy;
            }
        }
        for i in 0..self.boundary.points.len() {
            let dx = self.interior.points[index].x - self.boundary.points[i].x;
            let dy = self.interior.points[index].y - self.boundary.points[i].y;
            let r = (dx * dx + dy * dy).sqrt();
            let r2 = r * r;
            let r4 = r2 * r2;
            let r8 = r4 * r4;
            let r13 = r8 * r4 * r;
            let r7 = r2 * r4 * r;
            f_x += 4.0 * epsilon * (12.0 * sigma12 / r13 - 6.0 * sigma6 / r7) * dx;
            f_y += 4.0 * epsilon * (12.0 * sigma12 / r13 - 6.0 * sigma6 / r7) * dy;
        }
        kd_tree::Grid2D { x: f_x, y: f_y }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wave_init() {
        let wave = WaveEq::new();
        //println!("{:#?}", wave);
        assert_eq!(wave.boundary.points.len(), 0);
        assert_eq!(wave.interior.points.len(), 0);
    }

    #[test]
    fn wave_init_boundary_interior() {
        let mut wave = WaveEq::new();
        wave.create(25);
        for i in 0..wave.interior.points.len() {
            println!(
                "{} {}",
                wave.interior.points[i].x, wave.interior.points[i].y
            )
        }
        for i in 0..wave.boundary.points.len() {
            println!(
                "{} {}",
                wave.boundary.points[i].x, wave.boundary.points[i].y
            )
        }
        assert_eq!(wave.boundary.points.len(), 25);
        assert_eq!(wave.interior.points.len(), 49);
    }

    #[test]
    fn set_initial() {
        let mut wave = WaveEq::new();
        wave.create(25);
        wave.set_initial_condition();
        assert_eq!(wave.value.len(), 49);
        assert_eq!(wave.boundary.points.len(), 25);
        assert_eq!(wave.interior.points.len(), 49);
    }

    #[test]
    fn set_near_boundary_points() {
        let mut wave = WaveEq::new();
        wave.create(25);
        wave.set_initial_condition();
        wave.set_boundary_near_points();
        assert_eq!(wave.near.len(), 25);
    }

    #[test]
    fn set_near_interior_points() {
        let mut wave = WaveEq::new();
        wave.create(25);
        wave.set_initial_condition();
        wave.set_boundary_near_points();
        wave.set_interior_near_points();
        assert_eq!(wave.near_points_interior.len(), 49);
    }

    #[test]
    fn set_near_interior_points_plus_boundary() {
        let mut wave = WaveEq::new();
        wave.create(25);
        wave.set_initial_condition();
        wave.set_boundary_near_points();
        wave.set_interior_near_points();
        wave.set_interior_near_points_plus_boundary();
        let mut x = 0.0;
        for i in 0..wave.near_points_boundary[0].len() {
            let index = wave.near_points_boundary[0][i];
            //println!("{} {}", wave.boundary.points[index].x, wave.boundary.points[index].y);
            x = wave.boundary.points[index].x;
        }
        assert_eq!(x, -0.8090169943749473);
    }

    #[test]
    fn set_near_poly_fit() {
        let mut wave = WaveEq::new();
        wave.create(22);
        wave.set_initial_condition();
        wave.set_boundary_near_points();
        wave.set_interior_near_points();
        wave.set_interior_near_points_plus_boundary();
        wave.init_poly(3);
        wave.set_init_poly(1.0e-6);
        let x = wave.poly_eval(-1.0, 0.0);
        assert_eq!(x, 0.00000000769584173687417);
        let x = wave.poly_eval(-0.9, 0.0);
        assert_eq!(x, 7.181925602865518e-9);
        let x = wave.poly_eval(-0.8, 0.0);
        assert_eq!(x, 6.712246014828336e-9);
        let x = wave.poly_eval(-0.7, 0.0);
        assert_eq!(x, 6.284063050153314e-9);
        let x = wave.poly_eval(-0.6, 0.0);
        assert_eq!(x, -0.0017798942283290195);
        let x = wave.poly_eval(-0.5, 0.0);
        assert_eq!(x, 0.0017354278901027462);
        let x = wave.poly_eval(-0.4, 0.0);
        assert_eq!(x, 0.007959655646481115);
        let x = wave.poly_eval(-0.3, 0.0);
        assert_eq!(x, -0.07064301492441283);
        let x = wave.poly_eval(-0.2, 0.0);
        assert_eq!(x, 0.15512311432500947);
        let x = wave.poly_eval(-0.1, 0.0);
        assert_eq!(x, 0.8552289444067229);
        let x = wave.poly_eval(0.0, 0.0);
        assert_eq!(x, 0.9614778960644093);
    }
}
