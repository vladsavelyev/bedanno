use std::io;
use std::env;
use std::fs;
use bedanno;

fn main() {
    env::args().next();

    let input: Box<dyn io::Read> = match env::args().next() {
        Some(path) => { Box::new(fs::File::open(path).unwrap()) }
        None => { Box::new(io::stdin()) }
    };

    let reference: Box<dyn io::Read> = match env::args().next() {
        Some(path) => { Box::new(fs::File::open(path).unwrap()) }
        None => { Box::new(io::stdin()) }
    };

    let output = io::stdout();

    bedanno::run(
        input,
        reference,
        output,
    ).expect("Should be able to annotate BED file")
}
