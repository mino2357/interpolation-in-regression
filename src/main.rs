mod grid_3d;
mod kd_tree;
mod point;
mod two_variable_polynomial;
mod visualization;
mod wave_eqation;

fn main() {
    let mut a = wave_eqation::WaveEq::make(50, 2, 1.0e-9);
    let n = 25;
    for i in -n..n {
        for j in -n..n {
            let x = i as f64 / n as f64;
            let y = j as f64 / n as f64;
            if x * x + y * y < 1.0 {
                println!("{} {} {}", x, y, a.poly_eval(x, y));
            }
        }
    }
}
