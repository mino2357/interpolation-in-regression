use rand::prelude::*;
use plotters::prelude::*;
use apng::{load_dynamic_image, Encoder, Frame, PNGImage};

use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;

mod vector2;
mod grid_2d;

const NUM_POINTS: usize = 100;

pub fn draw_graph(points: &mut grid_2d::Grid2D, poly: &Vec<f64>, counter: i32) {

    let out_file_name = format!("{:04}", counter).to_string() + ".png";

    let x_max = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.x.max(m));
    let x_min = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.x.min(m));
    let y_max = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.y.max(m));
    let y_min = points.point_2d.iter().fold(0.0 / 0.0, |m, v| v.y.min(m));

    // draw
    let root = BitMapBackend::new(&out_file_name, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let graph_margin_x = 0.1 * (x_max - x_min);
    let graph_margin_y = 0.1 * (y_max - y_min);

    // setting graph
    let mut chart = ChartBuilder::on(&root)
        .caption("regression", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            ((x_min - graph_margin_x) as f32)..((x_max + graph_margin_x) as f32),
            ((y_min - graph_margin_y) as f32)..((y_max + graph_margin_y) as f32),
        )
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    let margin = 0.02 * (x_max - x_min);

    chart
        .draw_series(LineSeries::new(
            (0..=1000)
                .map(|x| (x_min - margin) + x as f64 * (x_max - x_min + 2.0 * margin) / 1000.0)
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

fn main() {
    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let mut points = grid_2d::Grid2D::new();
    for _ in 0..NUM_POINTS {
        let tmp = 1.9 * (rng.gen::<f64>() - 0.5);
        points.point_2d.push(vector2::Vector2 {
            x: tmp,
            y: (8.0 * tmp).sin(), // tmp * tmp * tmp + 4.0 * tmp * tmp - 2.0 * tmp + 1.0, // 2.0 * (rng.gen::<f64>() - 0.5),
        });
    }

    let mut poly = vec![0.0; 13];

    let mut tol = 0.95 * points.potential(&poly);

    let mut counter = 0;
    let max_counter = 130;
    
    loop {
        points.poly_fitting_by_euler_with_tol(&mut poly, tol);
        //points.poly_fitting_by_classical_rk4_with_tol(&mut poly, tol);
    
        println!("coef: {:?}", poly);

        draw_graph(&mut points, &poly, counter);
        tol = 0.95 * tol;
        println!("{}", format!("{:04}", counter).to_string() + ".png");
        counter += 1;
        if tol < 1.0e-4 || counter == max_counter {
            break;
        }
    }

    let mut files = vec![];

    for i in 0..max_counter {
        files.push(format!("{:04}", i).to_string() + ".png");
    }

    let mut png_images: Vec<PNGImage> = Vec::new();

    for f in files.iter() {
        let mut file = File::open(f).unwrap();
        let mut buffer = vec![];
        file.read_to_end(&mut buffer).unwrap();
        let img = image::load_from_memory(&buffer).unwrap();
        png_images.push(load_dynamic_image(img).unwrap());
    }

    let path = Path::new(r"regression.png");
    let mut out = BufWriter::new(File::create(path).unwrap());

    let config = apng::create_config(&png_images, None).unwrap();
    let mut encoder = Encoder::new(&mut out, config).unwrap();

    for image in png_images.iter() {
        let frame = Frame {
            delay_num: Some(1),
            delay_den: Some(8), // 2, 3, 4, 5, 6, 7
            ..Default::default()
        };
        encoder.write_frame(image, frame).unwrap();
    }

    match encoder.finish_encode() {
        Ok(_n) => println!("success"),
        Err(err) => eprintln!("{}", err),
    }
}
