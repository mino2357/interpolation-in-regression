use rand::prelude::*;
use plotters::prelude::*;

use crate::euler::Coef;

mod euler;

const OUT_FILE_NAME: &str = "test.png";
const NUM_POINTS: usize = 1000;
const GRAPH_MARGIN: f64 = 0.1;
const COEF_A: f64 = -2.0;
const COEF_B: f64 = -0.5;
const COEF_C: f64 = -0.6;
const COEF_D: f64 = 12.0;
const DIM: usize = 5;

fn main() {
    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let mut vec = Vec::new();
    for _ in 0..NUM_POINTS {
        let tmp = 2.0 * (rng.gen::<f64>() - 0.5);
        // Vector2. ref. https://docs.rs/cgmath/latest/cgmath/struct.Point2.html
        vec.push(cgmath::Vector2::new(
            tmp, // + 0.2 * (rng.gen::<f64>() - 0.5),
            COEF_A * tmp * tmp * tmp
                + COEF_B * tmp * tmp
                + COEF_C * tmp
                + COEF_D
                + 0.1 * (rng.gen::<f64>() - 0.5),
        ));
    }

    let mut coef: Coef = euler::Coef {
        coef: Coef::new(DIM),
    };

    coef.euler(&vec);

    println!("{:?}", coef);

    // ref. https://qiita.com/lo48576/items/343ca40a03c3b86b67cb
    let x_max = vec.iter().fold(0.0 / 0.0, |m, v| v.x.max(m));
    let x_min = vec.iter().fold(0.0 / 0.0, |m, v| v.x.min(m));
    let y_max = vec.iter().fold(0.0 / 0.0, |m, v| v.y.max(m));
    let y_min = vec.iter().fold(0.0 / 0.0, |m, v| v.y.min(m));

    // draw
    let root = BitMapBackend::new(OUT_FILE_NAME, (1200, 800)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    // setting graph
    let mut chart = ChartBuilder::on(&root)
        .caption("test", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            ((x_min - GRAPH_MARGIN) as f32)..((x_max + GRAPH_MARGIN) as f32),
            ((y_min - GRAPH_MARGIN) as f32)..((y_max + GRAPH_MARGIN) as f32),
        )
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (-500..=500).map(|x| x as f64 / 500.0).map(|x| {
                (
                    x as f32,
                    (coef.coef[3] * x * x * x
                        + coef.coef[2] * x * x
                        + coef.coef[1] * x
                        + coef.coef[0]) as f32,
                )
            }),
            &BLUE,
        ))
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            (-500..=500).map(|x| x as f64 / 500.0).map(|x| {
                (
                    x as f32,
                    (COEF_A * x * x * x + COEF_B * x * x + COEF_C * x + COEF_D) as f32,
                )
            }),
            &GREEN,
        ))
        .unwrap();

    chart
        .draw_series(PointSeries::of_element(
            (0..NUM_POINTS).map(|i| (vec[i].x as f32, vec[i].y as f32)),
            1,
            ShapeStyle::from(&RED).filled(),
            &|coord, size, style| EmptyElement::at(coord) + Circle::new((0, 0), size, style),
        ))
        .unwrap();
}
