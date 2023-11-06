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
    bam_file: String,

    #[arg(short = 'd', long, default_value = "99.0")]
    min_identity: String,
}

fn main() {
    let arguments = Arguments::parse();
    let min_identity: f64 = arguments.min_identity.parse().unwrap_or(99.0);
    if let Err(error) = realignment::run(
        &arguments.input_bed,
        &arguments.support_read_file,
        &arguments.bam_file,
        min_identity,
    ) {
        eprintln!("{}", error);
        process::exit(1);
    }
}
