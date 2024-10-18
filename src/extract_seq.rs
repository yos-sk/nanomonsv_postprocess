use rust_htslib::faidx;
use std::cmp::{max, min};
use std::error::Error;
use std::io::BufRead;

use nanomonsv_postprocess::open_file;

pub fn run(input_bed: &str, reference_fasta: &str) -> Result<(), Box<dyn Error>> {
    let bed_reader = open_file(input_bed).expect(&format!("Could not open file {}", input_bed));
    let fasta_reader = faidx::Reader::from_path(reference_fasta)
        .expect(&format!("Could not open file {}", reference_fasta));

    for line in bed_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let contig = split_line[0];
        let pos = split_line[1].parse::<isize>().unwrap();
        let start = max(0, pos - 100) as usize;
        let length = fasta_reader.fetch_seq_len(contig) as isize;
        let end = min(pos + 101, length) as usize;

        let seq = fasta_reader
            .fetch_seq_string(contig, start, end)
            .expect(&format!("Could not fetch {}:{}-{}", contig, start, end));

        println!("{}\t{}", line, seq);
    }
    Ok(())
}
