mod grid_3d;
mod kd_tree;
mod point;
mod two_variable_polynomial;
mod visualization;
mod wave_eqation;

fn main() {
    let mut wave: wave_eqation::WaveEq = wave_eqation::WaveEq::make(40, 2, 1.0e-9, 4);
    let n = 15;
    let mut vec = grid_3d::Grid3D::new();
    for i in -n..n {
        for j in -n..n {
            let x = i as f64 / n as f64;
            let y = j as f64 / n as f64;
            if x * x + y * y < 1.0 {
                println!("{} {} {}", x, y, wave.poly_eval(x, y));
                let p = point::Point3::new(x, y, wave.poly_eval(x, y));
                vec.push(p);
            }
        }
        println!("");
    }

    visualization::draw_3d_points(&vec, 0);
}
