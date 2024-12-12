use clap::{Parser, Subcommand};
use std::process;

mod extract_seq;
mod filt;
mod merge;
mod realignment;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.2.5", about = "Post process of nanomonsv", long_about = None)]

struct Arguments {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    ExtractSeq {
        #[arg(short = 'b', long)]
        input_bed: String,

        #[arg(short = 'f', long)]
        reference_fasta: String,
    },

    Realignment {
        #[arg(short = 'i', long)]
        input_bed: String,

        #[arg(short = 's', long)]
        support_read_file: String,

        #[arg(short = 'b', long)]
        bam_file: String,

        #[arg(short = 'd', long, default_value = "99.0")]
        min_identity: String,

        #[arg(short = 'l', long, default_value = "180")]
        min_length: String,
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

    Merge {
        #[arg(short = 'i', long)]
        input_file_1: String,

        #[arg(short = 'j', long)]
        input_file_2: String,

        #[arg(short = 'o', long)]
        output_file: String,

        #[arg(short = 'm', long, default_value = "50")]
        coord_margin: String,

        #[arg(short = 's', long, default_value = "98.0")]
        inserted_seq_identity: String,
    },
}

fn main() {
    let arguments = Arguments::parse();
    match &arguments.command {
        Commands::ExtractSeq {
            input_bed,
            reference_fasta,
        } => {
            if let Err(error) = extract_seq::run(input_bed, reference_fasta) {
                eprintln!("{}", error);
                process::exit(1);
            }
        }

        Commands::Realignment {
            input_bed,
            support_read_file,
            bam_file,
            min_identity,
            min_length,
        } => {
            let min_identity: f64 = min_identity.parse().unwrap_or(99.0);
            let min_length: usize = min_length.parse().unwrap_or(180);
            if let Err(error) = realignment::run(
                input_bed,
                support_read_file,
                bam_file,
                min_identity,
                min_length,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        }

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
        }

        Commands::Merge {
            input_file_1,
            input_file_2,
            output_file,
            coord_margin,
            inserted_seq_identity,
        } => {
            let coord_margin: isize = coord_margin.parse().unwrap_or(50);
            let inserted_seq_identity: f64 = inserted_seq_identity.parse().unwrap_or(98.0);
            if let Err(error) = merge::run(
                input_file_1,
                input_file_2,
                output_file,
                coord_margin,
                inserted_seq_identity,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        }
    }
}
