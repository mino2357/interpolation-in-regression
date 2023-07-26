use rand::prelude::*;
use plotters::prelude::*;

mod vector2;
mod grid_2d;

const OUT_FILE_NAME: &str = "test.png";
const NUM_POINTS: usize = 8;
const GRAPH_MARGIN: f64 = 0.1;

fn main() {
    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let mut points = grid_2d::Grid2D::new();
    for _ in 0..NUM_POINTS {
        let tmp = 2.0 * (rng.gen::<f64>() - 0.5);
        points.point_2d.push(vector2::Vector2 {
            x: tmp,
            y: tmp * tmp * tmp + 4.0 * tmp * tmp - 2.0 * tmp + 1.0, // 2.0 * (rng.gen::<f64>() - 0.5),
        });
    }

    let mut poly = vec![0.0; NUM_POINTS];

    points.poly_fitting_by_euler(&mut poly);
    //points.poly_fitting_by_modified_euler(&mut poly);
    //points.poly_fitting_by_classical_rk4(&mut poly);

    println!("end: {:?}", poly);

    // ref. https://qiita.com/lo48576/items/343ca40a03c3b86b67cb
    let x_max = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.x.max(m));
    let x_min = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.x.min(m));
    let y_max = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.y.max(m));
    let y_min = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.y.min(m));

    // draw
    let root = BitMapBackend::new(OUT_FILE_NAME, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    // setting graph
    let mut chart = ChartBuilder::on(&root)
        .caption("regression", ("sans-serif", 50).into_font())
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
                .map(|x| (x as f32, points.eval(&poly, x) as f32)),
            &BLUE,
        ))
        .unwrap();

    chart
        .draw_series(PointSeries::of_element(
            (0..NUM_POINTS).map(|i| (points.point_2d[i].x as f32, points.point_2d[i].y as f32)),
            2,
            ShapeStyle::from(&RED).filled(),
            &|coord, size, style| EmptyElement::at(coord) + Circle::new((0, 0), size, style),
        ))
        .unwrap();
}
