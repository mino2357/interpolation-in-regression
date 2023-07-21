#[allow(dead_code)]
fn twice(x: f64) -> f64 {
    2.0 * x
}

// ref. https://bablovia.hatenablog.com/entry/2020/11/17/165344
#[allow(dead_code)]
fn func(f: Box<dyn Fn(f64) -> f64>) -> f64 {
    f(2.0)
}

#[allow(dead_code)]
fn id(x: f64) -> f64 {
    x
}

#[allow(dead_code)]
fn test_euler(mut x: f64, f: Box<dyn Fn(f64) -> f64>) -> f64 {
    let dt = 1.0e-3;
    for _i in 0..1000 {
        x = x + dt * f(x);
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(func(Box::new(twice)), 4.0);
    }

    #[test]
    fn euler() {
        assert_eq!(test_euler(1.0, Box::new(id)), 2.716923932235896);
    }
}
