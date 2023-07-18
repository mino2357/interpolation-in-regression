use rand::prelude::*;
use plotters::prelude::*;

const OUT_FILE_NAME: &str = "test.png";
const NUM_POINTS: usize = 10000;
const MAX_ITER: usize = 10000000;
const GRAPH_MARGIN: f64 = 0.1;
const COEF_A: f64 = 1.0;
const COEF_B: f64 = 2.0;
const COEF_C: f64 = 20.0;

// ax^2+bx+c
#[derive(Debug)]
struct Coef2 {
    a: f64,
    b: f64,
    c: f64,
}

impl Coef2 {
    fn euler(self: &mut Coef2, vec: &Vec<cgmath::Vector2<f64>>) {
        let dt = 1.0e-4;
        for _i in 0..MAX_ITER {
            let pre = self.err_square(vec);
            self.a += vec.iter().map(|v| -dt * self.dUda(v.x, v.y)).sum::<f64>();
            self.b += vec.iter().map(|v| -dt * self.dUdb(v.x, v.y)).sum::<f64>();
            self.c += vec.iter().map(|v| -dt * self.dUdc(v.x, v.y)).sum::<f64>();
            let post = self.err_square(vec);
            //println!("{:?}, {:?}", pre, post);
            if (pre - post).powf(2.0) < 0.001 * dt * dt {
                break;
            }
        }
        println!("{:?}", self.err_square(vec));
    }

    fn err_square(self: &Coef2, vec: &Vec<cgmath::Vector2<f64>>) -> f64 {
        vec.iter()
            .map(|v| (self.a * v.x * v.x + self.b * v.x + self.c - v.y).powf(2.0))
            .sum::<f64>()
            / NUM_POINTS as f64
    }

    #[allow(non_snake_case)]
    fn dUda(self: &Coef2, x: f64, y: f64) -> f64 {
        -2.0 * x * x * (y - self.a * x * x - self.b * x - self.c)
    }

    #[allow(non_snake_case)]
    fn dUdb(self: &Coef2, x: f64, y: f64) -> f64 {
        -2.0 * x * (y - self.a * x * x - self.b * x - self.c)
    }

    #[allow(non_snake_case)]
    fn dUdc(self: &Coef2, x: f64, y: f64) -> f64 {
        -2.0 * (y - self.a * x * x - self.b * x - self.c)
    }
}

fn main() {
    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let mut vec = Vec::new();
    for _ in 0..NUM_POINTS {
        let tmp = 2.0 * (rng.gen::<f64>() - 0.5);
        // Vector2. ref. https://docs.rs/cgmath/latest/cgmath/struct.Point2.html
        vec.push(cgmath::Vector2::new(
            tmp + 0.2 * (rng.gen::<f64>() - 0.5),
            COEF_A * tmp * tmp + COEF_B * tmp + COEF_C + 0.6 * (rng.gen::<f64>() - 0.5),
        ));
    }

    // ref. https://qiita.com/lo48576/items/343ca40a03c3b86b67cb
    let x_max = vec.iter().fold(0.0 / 0.0, |m, v| v.x.max(m));
    let x_min = vec.iter().fold(0.0 / 0.0, |m, v| v.x.min(m));
    let y_max = vec.iter().fold(0.0 / 0.0, |m, v| v.y.max(m));
    let y_min = vec.iter().fold(0.0 / 0.0, |m, v| v.y.min(m));

    let mut coef = Coef2 {
        a: 0.0,
        b: 0.0,
        c: 0.0,
    };

    coef.euler(&vec);

    println!("{:?}", coef);

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
            (-500..=500)
                .map(|x| x as f64 / 500.0)
                .map(|x| (x as f32, (coef.a * x * x + coef.b * x + coef.c) as f32)),
            &BLUE,
        ))
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            (-500..=500)
                .map(|x| x as f64 / 500.0)
                .map(|x| (x as f32, (COEF_A * x * x + COEF_B * x + COEF_C) as f32)),
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
