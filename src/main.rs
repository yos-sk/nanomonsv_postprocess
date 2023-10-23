use clap::Parser;
use std::process;

mod realignment;


#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.1.0", about = "Post process of nanomonsv", long_about = None)]
struct Arguments {
    #[arg(short = 'b', long)]
    input_bed: String
}

fn main() {
    let arguments = Arguments::parse();
    if let Err(error) = realignment::run(
        &arguments.input_bed,
    ) {
        eprintln!("{}", error);
        process::exit(1);
    }
}
