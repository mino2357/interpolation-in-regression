use super::point;

#[derive(Debug)]
pub struct Grid3D {
    pub points_3d: Vec<point::Point3>,
}

impl Grid3D {
    #[allow(dead_code)]
    pub fn new() -> Self {
        Grid3D {
            points_3d: vec![],
        }
    }

}