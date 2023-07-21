use rand::prelude::*;
use plotters::prelude::*;
use crate::polynomial_interpolation::Coef;

mod embedded_runge_kutta;
mod polynomial_interpolation;

const OUT_FILE_NAME: &str = "test.png";
const NUM_POINTS: usize = 1000;
const GRAPH_MARGIN: f64 = 0.1;

fn twice(x: f64) -> f64 {
    2.0 * x
}

fn func(f: Box<dyn Fn(f64) -> f64>) -> f64 {
    f(2.0)
}

fn main() {
    println!("{:?}", func(Box::new(twice)));

    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let mut vec = Vec::new();
    for _ in 0..NUM_POINTS {
        let tmp = 2.0 * (rng.gen::<f64>() - 0.5);
        // Vector2. ref. https://docs.rs/cgmath/latest/cgmath/struct.Point2.html
        vec.push(cgmath::Vector2::new(
            1.0 * tmp, // + 0.01 * (rng.gen::<f64>() - 0.5),
            (4.0 * tmp).sin() + 0.1 * (rng.gen::<f64>() - 0.5),
        ));
    }

    let mut coef: Coef = polynomial_interpolation::Coef::new(7);

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
            (0..=1000)
                .map(|x| x_min + x as f64 * (x_max - x_min) / 1000.0)
                .map(|x| (x as f32, coef.eval(x) as f32)),
            &BLUE,
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
