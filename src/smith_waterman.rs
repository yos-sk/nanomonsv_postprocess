use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;

use std::error::Error;

pub fn run(seq1: &Vec<u8>, seq2: &Vec<u8>) -> Result<f64, Box<dyn Error>> {

    let score = |a: u8, b: u8| if a == b { 1i32 } else { -2i32 };
    // Gap open score: -2, gap extension score: -1 
    let mut aligner = Aligner::with_capacity(seq1.len(), seq2.len(), -2, -1, &score);
    let alignment = aligner.global(seq1, seq2);

    let mut m = 0;
    let mut d = 0;
    for stat in alignment.operations.iter() {
        if *stat == Match {
            m += 1;
        } else {
            d += 1;
        }
    }

    let identity: f64 = m as f64 / (m as f64 + d as f64) * 100.0;
    //println!("Smith-Waterman result: {:?}", alignment);
    //println!("Identity: {}", identity);
    Ok(identity)
}

// TODO: byte列バージョンに書き換える必要あり。
pub fn reverse_complement(sequence: &str) -> String {
    // complement
    let complement = sequence
        .chars()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => c,
        })
        .collect::<String>();

    // reverse
    let rev_comp = complement.chars().rev().collect::<String>();

    rev_comp
}
