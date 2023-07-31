use rand::prelude::*;
use plotters::prelude::*;
use apng::{load_dynamic_image, Encoder, Frame, PNGImage};
use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;

mod point;
mod grid_2d;
mod grid_3d;
mod two_variable_polynomial;

fn draw_3d_graph(poly: &two_variable_polynomial::TwoPolynomial, points: &grid_3d::Grid3D, counter: i32) {
    let out_file_name = format!("{:04}", counter).to_string() + ".png";

    let root = BitMapBackend::new(&out_file_name, (2560, 1440)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let chart_res: i32 = 100;

    let z_max = points.points_3d.iter().fold(0.0 / 0.0, |m, v| v.z.max(m));
    let z_min = points.points_3d.iter().fold(0.0 / 0.0, |m, v| v.z.min(m));

    let z_margin = 10.0 * (z_max - z_min);

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption("Regression", ("sans-serif", 40))
        .build_cartesian_3d(-1.2..1.2, (z_min - z_margin)..(z_max + z_margin), -1.2..1.2)
        .unwrap();

    chart.with_projection(|mut pb| {
        pb.yaw = 0.1 + 0.005 * counter as f64;
        pb.into_matrix()
    });

    chart.configure_axes().draw().unwrap();

    let mut data = vec![];

    for x in (-chart_res..chart_res).map(|v| v as f64 / chart_res as f64) {
        let mut row = vec![];
        for y in (-chart_res..chart_res).map(|v| v as f64 / chart_res as f64) {
            row.push((x, poly.eval_xy(x, y), y));
        }
        data.push(row);
    }

    chart.draw_series(
        (0..(2*chart_res-1))
            .map(|x| std::iter::repeat(x).zip(0..(2*chart_res-1)))
            .flatten()
            .map(|(x,z)| {
                Polygon::new(vec![
                    data[x as usize][z as usize],
                    data[(x+1) as usize][z as usize],
                    data[(x+1) as usize][(z+1) as usize],
                    data[x as usize][(z+1) as usize],
                ], &BLUE.mix(0.3))
            })
    ).unwrap();

    chart
    .draw_series(PointSeries::of_element(
        (0..points.points_3d.len()).map(|i| (points.points_3d[i].x, points.points_3d[i].z, points.points_3d[i].y)),
        4,
        ShapeStyle::from(&RED).filled(),
        &|coord, size, style| EmptyElement::at(coord) + Circle::new((0, 0), size, style),
    ))
    .unwrap();
}

fn main() {
    let mut poly = two_variable_polynomial::TwoPolynomial::new(6);

    let seed: [u8; 32] = [1; 32];
    let mut rng: rand::rngs::StdRng = rand::SeedableRng::from_seed(seed);

    let data_size: usize = 21;
    let mut points = grid_3d::Grid3D::new();
    for _ in 0..data_size {
        let x_r = 2.0 * (rng.gen::<f64>() - 0.5);
        let y_r = 2.0 * (rng.gen::<f64>() - 0.5);
        let z_r = 2.0 * (rng.gen::<f64>() - 0.5);
        points.points_3d.push(point::Point3 {
            x: x_r,
            y: y_r,
            z: z_r, //(4.0 * x_r).sin() * (5.0 * y_r).cos() // 0.0 + 0.0 * x_r + 0.0 * y_r + 1.0 * x_r * x_r + 0.0 * x_r * y_r - 1.0 * y_r * y_r,
        });
    }

    let ratio = 0.99;

    let mut tol = ratio * points.potential(&poly);

    let mut counter = 0;
    let max_counter = 1200;

    loop {
        points.poly_fitting_by_euler_with_tol(&mut poly, tol);
    
        //println!("tol: {:?}, coef: {:?}", tol, poly);
        println!("tol: {:?}", tol);

        draw_3d_graph(&poly, &points, counter);
        tol = ratio * tol;
        
        println!("{}", format!("{:04}", counter).to_string() + ".png");
        
        counter += 1;
        if counter == max_counter {
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

    let path = Path::new(r"regression2.png");
    let mut out = BufWriter::new(File::create(path).unwrap());

    let config = apng::create_config(&png_images, None).unwrap();
    let mut encoder = Encoder::new(&mut out, config).unwrap();

    for image in png_images.iter() {
        let frame = Frame {
            delay_num: Some(1),
            delay_den: Some(30), // 2, 3, 4, 5, 6, 7
            ..Default::default()
        };
        encoder.write_frame(image, frame).unwrap();
    }

    match encoder.finish_encode() {
        Ok(_n) => println!("success"),
        Err(err) => eprintln!("{}", err),
    }
}
