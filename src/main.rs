use rand::prelude::*;

mod point;
mod grid_2d;
mod grid_3d;
mod two_variable_polynomial;
mod visualization;

fn main() {
    let mut poly = two_variable_polynomial::TwoPolynomial::new(6);

    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let data_size: usize = 21;
    let mut points = grid_3d::Grid3D::new();
    for _ in 0..data_size {
        let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
        let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
        let z_r = 2.0 * (rng.gen::<f64>() - 0.5);
        points.points_3d.push(point::Point3 {
            x: x_r,
            y: y_r,
            z: z_r,
        });
    }

    let ratio = 0.99;
    let mut tol = ratio * points.potential(&poly);
    let mut counter = 0;
    let max_counter = 10;

    loop {
        visualization::draw_3d_graph(&poly, &points, counter);
        points.poly_fitting_by_euler_with_tol(&mut poly, tol);
        println!("tol: {:?}", tol);
        println!("{}", format!("{:04}", counter).to_string() + ".png");
        
        tol = ratio * tol;
        counter += 1;
        if counter == max_counter {
            break;
        }
    }
    visualization::gen_apng(max_counter);
}
