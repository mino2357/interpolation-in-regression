mod grid_3d;
mod kd_tree;
mod point;
mod two_variable_polynomial;
mod visualization;
mod wave_eqation;

fn main() {
    let dt = 1.0e-3;
    let tol = 1.0e-4;
    let degree = 2;
    let num_poiunt = 70;
    let num_neighbor = 5;
    let mut wave: wave_eqation::WaveEq = wave_eqation::WaveEq::make(num_poiunt, degree, tol, num_neighbor);

    for t in 0..1000 {
        let n = 25;
        let mut vec = grid_3d::Grid3D::new();
        for i in -n..n {
            for j in -n..n {
                let x = i as f64 / n as f64;
                let y = j as f64 / n as f64;
                if x * x + y * y < 1.0 {
                    let p = point::Point3::new(x, y, wave.poly_eval(x, y));
                    vec.push(p);
                }
            }
        }
        
        let mut p_in = grid_3d::Grid3D::new();
        for i in 0..wave.interior.points.len() {
            let p = point::Point3::new(wave.interior.points[i].x, wave.interior.points[i].y, wave.value[i]);
            p_in.push(p);
        }

        visualization::draw_3d_points(&vec, &p_in, t);
        for j in 0..10 {
            println!("{} {}", t, j);
            wave.step(tol, dt);
        }
    }
}
