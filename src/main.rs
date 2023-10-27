use clap::Parser;
use std::process;

mod realignment;


#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.1.0", about = "Post process of nanomonsv", long_about = None)]
struct Arguments {
    #[arg(short = 'i', long)]
    input_bed: String,

    #[arg(short = 's', long)]
    support_read_file: String,

    #[arg(short = 'b', long)]
    nanomonsv_bp_file: String,
}

fn main() {
    let arguments = Arguments::parse();
    if let Err(error) = realignment::run(
        &arguments.input_bed,
        &arguments.support_read_file,
        &arguments.nanomonsv_bp_file,
    ) {
        eprintln!("{}", error);
        process::exit(1);
    }
}
