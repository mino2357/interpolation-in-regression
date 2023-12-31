#[derive(Debug)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    #[allow(dead_code)]
    pub fn new(x_: f64, y_: f64, z_: f64) -> Self {
        Point3 {
            x: x_,
            y: y_,
            z: z_,
        }
    }
}
