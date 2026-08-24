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

        /// Contig names of haplotype 1 (bam_refiner output). Without both
        /// lists the haplotype-combination constraint on grouping is skipped.
        #[arg(long, default_value = "")]
        hap1_list: String,

        /// Contig names of haplotype 2.
        #[arg(long, default_value = "")]
        hap2_list: String,

        /// When two calls share a haplotype on one side, how far apart that
        /// breakpoint may be before they are treated as separate loci.
        #[arg(long, default_value = "100")]
        max_position_diff: String,
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
            hap1_list,
            hap2_list,
            max_position_diff,
        } => {
            let min_identity: f64 = min_identity.parse().unwrap_or(99.0);
            let min_length: usize = min_length.parse().unwrap_or(180);
            let max_position_diff: usize = max_position_diff.parse().unwrap_or(100);
            if let Err(error) = realignment::run(
                input_bed,
                support_read_file,
                bam_file,
                min_identity,
                min_length,
                hap1_list,
                hap2_list,
                max_position_diff,
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
