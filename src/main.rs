use clap::{Parser, Subcommand};
use std::process;

mod realignment;
mod filt;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.1.0", about = "Post process of nanomonsv", long_about = None)]

struct Arguments {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Realignment {
        #[arg(short = 'i', long)]
        input_bed: String,

        #[arg(short = 's', long)]
        support_read_file: String,

        #[arg(short = 'b', long)]
        bam_file: String,

        #[arg(short = 'd', long, default_value = "99.0")]
        min_identity: String,
    },

    Filt {
        #[arg(short = 'i', long)]
        identical_file: String,

        #[arg(short = 'n', long)]
        nanomonsv_result: String,

        #[arg(short = 's', long)]
        support_read_file: String,

        #[arg(short = 'b', long)]
        bam_file: String,
    },
}

fn main() {
    let arguments = Arguments::parse();
    match &arguments.command {
        Commands::Realignment {
            input_bed,
            support_read_file,
            bam_file,
            min_identity,
        } => {
            let min_identity: f64 = min_identity.parse().unwrap_or(99.0);
            if let Err(error) = realignment::run(
                input_bed,
                support_read_file,
                bam_file,
                min_identity,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }  
        },

        Commands::Filt {
            identical_file,
            nanomonsv_result,
            support_read_file,
            bam_file,
        } => {
            if let Err(error) = filt::run(
                identical_file,
                nanomonsv_result,
                support_read_file,
                bam_file,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }  
        },
    }
}
