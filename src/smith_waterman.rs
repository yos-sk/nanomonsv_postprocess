use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;

use std::error::Error;

/// Local-align `seq1` against `seq2` and its reverse complement.
///
/// Returns the identity (%) of the better orientation among those clearing both
/// thresholds, or `None` when neither does. The caller needs the value, not just
/// a yes/no, so that competing breakpoint pairs can be ranked against each other.
pub fn run(
    seq1: &Vec<u8>,
    seq2: &Vec<u8>,
    identity_th: f64,
    length_th: usize,
) -> Result<Option<f64>, Box<dyn Error>> {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -2i32 };
    // Gap open score: -2, gap extension score: -1
    let mut aligner1 = Aligner::with_capacity(seq1.len(), seq2.len(), -10, -1, &score);
    let alignment1 = aligner1.local(seq1, seq2);

    let mut m = 0;
    let mut d = 0;
    let mut length = 0;
    for stat in alignment1.operations.iter() {
        length += 1;
        if *stat == Match {
            m += 1;
        } else {
            d += 1;
        }
    }
    let identity1: f64 = m as f64 / (m as f64 + d as f64) * 100.0;
    let length1 = length;

    // reverse complement
    let seq2_c = reverse_complement(seq2);
    let mut aligner2 = Aligner::with_capacity(seq1.len(), seq2_c.len(), -10, -1, &score);
    let alignment2 = aligner2.local(seq1, &seq2_c);

    m = 0;
    d = 0;
    length = 0;
    for stat in alignment2.operations.iter() {
        length += 1;
        if *stat == Match {
            m += 1;
        } else {
            d += 1;
        }
    }
    let identity2: f64 = m as f64 / (m as f64 + d as f64) * 100.0;
    let length2 = length;

    let forward = identity1 >= identity_th && length1 >= length_th;
    let revcomp = identity2 >= identity_th && length2 >= length_th;

    Ok(match (forward, revcomp) {
        (true, true) => Some(identity1.max(identity2)),
        (true, false) => Some(identity1),
        (false, true) => Some(identity2),
        (false, false) => None,
    })
}

pub fn reverse_complement(sequence: &Vec<u8>) -> Vec<u8> {
    let mut complement: Vec<u8> = Vec::with_capacity(sequence.len());
    // complement
    for &base in sequence.iter().rev() {
        match base {
            b'A' => complement.push(b'T'),
            b'C' => complement.push(b'G'),
            b'G' => complement.push(b'C'),
            b'T' => complement.push(b'A'),
            _ => complement.push(base),
        }
    }

    complement
}
