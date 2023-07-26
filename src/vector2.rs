#[derive(Debug)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64,
}

impl Vector2 {
    #[allow(dead_code)]
    pub fn new(x_: f64, y_: f64) -> Self {
        Vector2 {
            x: x_,
            y: y_,
        }
    }
}