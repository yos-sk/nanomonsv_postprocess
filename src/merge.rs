use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use std::collections::HashSet;
use std::error::Error;
use std::io::BufRead;
use std::path::Path;
use std::fs::File;
use std::io::Write;

use nanomonsv_postprocess::open_file;

#[derive(Eq, Hash, PartialEq, Clone)]
struct NanomonsvInfo {
    chrom1: String,
    pos1: usize,
    dir1: String,
    chrom2: String,
    pos2: usize,
    dir2: String,
    inserted_seq: Vec<u8>,
    sv_id: String,
    checked_read_num_tumor: usize,
    supporting_read_num_tumor: usize,
    checked_read_num_control: usize,
    supporting_read_num_control: usize,
    is_filter: String,
    identical_svs: String,
    supporting_read_num_total: usize,
    haplotype_bp1: String,
    haplotype_bp2: String,
    sv_type: String,
}

pub fn run(input_file_1: &str, input_file_2: &str, output_file: &str, coord_margin: isize, inserted_seq_identity: f64) -> Result<(), Box<dyn Error>> {
    let file1_reader = open_file(input_file_1).expect(&format!("Could not open file {}", input_file_1));
    let mut file1_info_db: HashSet<NanomonsvInfo> = HashSet::new();
    let mut header = String::new();

    for (i, line) in file1_reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            header = line.to_string();
            continue;
        }
        let split_line: Vec<&str> = line.split('\t').collect();

        let info = NanomonsvInfo{
            chrom1: split_line[0].to_string(),
            pos1: split_line[1].parse::<usize>().unwrap(),
            dir1: split_line[2].to_string(),
            chrom2: split_line[3].to_string(),
            pos2: split_line[4].parse::<usize>().unwrap(),
            dir2: split_line[5].to_string(),
            inserted_seq: split_line[6].as_bytes().to_vec(),
            sv_id: split_line[7].to_string(),
            checked_read_num_tumor: split_line[8].parse::<usize>().unwrap(),
            supporting_read_num_tumor: split_line[9].parse::<usize>().unwrap(),
            checked_read_num_control: split_line[10].parse::<usize>().unwrap(),
            supporting_read_num_control: split_line[11].parse::<usize>().unwrap(),
            is_filter: split_line[12].to_string(),
            identical_svs: split_line[13].to_string(),
            supporting_read_num_total: split_line[14].parse::<usize>().unwrap(),
            haplotype_bp1: split_line[15].to_string(),
            haplotype_bp2: split_line[16].to_string(),
            sv_type: split_line[17].to_string(),
        };

        file1_info_db.insert(info);
    }

    let file2_reader = open_file(input_file_2).expect(&format!("Could not open file {}", input_file_2));
    let mut file2_info_db: HashSet<NanomonsvInfo> = HashSet::new();

    for (i, line) in file2_reader.lines().enumerate() {
        let line = line?;
        if i == 0 {
            continue;
        }
        let split_line: Vec<&str> = line.split('\t').collect();

        let chrom1 = split_line[0].to_string();
        let pos1 = split_line[1].parse::<usize>().unwrap();
        let chrom2 =  split_line[3].to_string();
        let pos2 = split_line[4].parse::<usize>().unwrap();

        if chrom1 > chrom2 || (chrom1 == chrom2 && pos1 > pos2) {
            let info = NanomonsvInfo{
                chrom1: split_line[3].to_string(),
                pos1: split_line[4].parse::<usize>().unwrap(),
                dir1: split_line[5].to_string(),
                chrom2: split_line[0].to_string(),
                pos2: split_line[1].parse::<usize>().unwrap(),
                dir2: split_line[2].to_string(),
                inserted_seq: split_line[6].as_bytes().to_vec(),
                sv_id: split_line[7].to_string(),
                checked_read_num_tumor: split_line[8].parse::<usize>().unwrap(),
                supporting_read_num_tumor: split_line[9].parse::<usize>().unwrap(),
                checked_read_num_control: split_line[10].parse::<usize>().unwrap(),
                supporting_read_num_control: split_line[11].parse::<usize>().unwrap(),
                is_filter: split_line[12].to_string(),
                identical_svs: split_line[13].to_string(),
                supporting_read_num_total: split_line[14].parse::<usize>().unwrap(),
                haplotype_bp1: split_line[16].to_string(),
                haplotype_bp2: split_line[15].to_string(),
                sv_type: split_line[17].to_string(),
            };
            file2_info_db.insert(info);
        } else {
            let info = NanomonsvInfo{
                chrom1: split_line[0].to_string(),
                pos1: split_line[1].parse::<usize>().unwrap(),
                dir1: split_line[2].to_string(),
                chrom2: split_line[3].to_string(),
                pos2: split_line[4].parse::<usize>().unwrap(),
                dir2: split_line[5].to_string(),
                inserted_seq: split_line[6].as_bytes().to_vec(),
                sv_id: split_line[7].to_string(),
                checked_read_num_tumor: split_line[8].parse::<usize>().unwrap(),
                supporting_read_num_tumor: split_line[9].parse::<usize>().unwrap(),
                checked_read_num_control: split_line[10].parse::<usize>().unwrap(),
                supporting_read_num_control: split_line[11].parse::<usize>().unwrap(),
                is_filter: split_line[12].to_string(),
                identical_svs: split_line[13].to_string(),
                supporting_read_num_total: split_line[14].parse::<usize>().unwrap(),
                haplotype_bp1: split_line[15].to_string(),
                haplotype_bp2: split_line[16].to_string(),
                sv_type: split_line[17].to_string(),
            };
            file2_info_db.insert(info);
        }

    }
    let _ = merge(&file1_info_db, &file2_info_db, output_file, coord_margin, inserted_seq_identity, &header);

    Ok(())
}

fn merge(db1: &HashSet<NanomonsvInfo>, db2: &HashSet<NanomonsvInfo>, output_file: &str, coord_margin: isize, inserted_seq_identity: f64, header: &str) -> Result<(), Box<dyn Error>> {
    let mut merged_info_db: HashSet<NanomonsvInfo> = HashSet::new();
    let mut merged_sv_id_1: HashSet<String> = HashSet::new();
    let mut merged_sv_id_2: HashSet<String> = HashSet::new();
    for info1 in db1 {
        for info2 in db2 {
            if info1.chrom1 != info2.chrom1 {
                continue;
            }

            if (info1.pos1 as isize - info2.pos1 as isize).abs() > coord_margin {
                continue;
            }

            if info1.dir1 != info2.dir1 {
                continue;
            }

            if info1.chrom2 != info2.chrom2 {
                continue;
            }

            if (info1.pos2 as isize - info2.pos2 as isize).abs() > coord_margin {
                continue;
            }

            if info1.dir2 != info2.dir2 {
                continue;
            }

            let ins_seq_aln = match smith_waterman(&info1.inserted_seq, &info2.inserted_seq, inserted_seq_identity) {
                Ok(flag) => flag,
                _ => false,
            };

            if !ins_seq_aln {
                continue;
            }
            
            // merge info1 and info2
            let merged_info = NanomonsvInfo {
                chrom1: info1.chrom1.clone(),
                pos1: info1.pos1,
                dir1: info1.dir1.clone(),
                chrom2: info1.chrom2.clone(),
                pos2: info1.pos2,
                dir2: info1.dir2.clone(),
                inserted_seq: info1.inserted_seq.clone(),
                sv_id: "1_".to_string() + &info1.sv_id + ",2_" + &info2.sv_id,
                checked_read_num_tumor: info1.checked_read_num_tumor + info2.checked_read_num_tumor,
                supporting_read_num_tumor: info1.supporting_read_num_tumor + info2.supporting_read_num_tumor,
                checked_read_num_control: info1.checked_read_num_control + info2.checked_read_num_control,
                supporting_read_num_control: info1.supporting_read_num_control + info2.supporting_read_num_control,
                is_filter: "1_".to_string() + &info1.is_filter + ",2_" + &info2.is_filter,
                identical_svs: info1.identical_svs.split(",").map(|s| format!("1_{}", s)).collect::<Vec<String>>().join(",") + "," + &(&info2.identical_svs.split(",").map(|s| format!("2_{}", s)).collect::<Vec<String>>().join(",")),
                supporting_read_num_total: info1.supporting_read_num_total + info2.supporting_read_num_total,
                haplotype_bp1: "1_".to_string() + &info1.haplotype_bp1 + ",2_" + &info2.haplotype_bp1,
                haplotype_bp2: "1_".to_string() + &info1.haplotype_bp2 + ",2_" + &info2.haplotype_bp2,
                sv_type: info1.sv_type.clone(),
            };

            merged_info_db.insert(merged_info);
            merged_sv_id_1.insert(info1.sv_id.clone());
            merged_sv_id_2.insert(info2.sv_id.clone());
        }
    }

    // output merge results
    let path = Path::new(output_file);
    let mut wfile = File::create(&path)?;
    writeln!(wfile, "{}", header)?;

    for info in merged_info_db {
        writeln!(wfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                info.chrom1,
                info.pos1,
                info.dir1,
                info.chrom2,
                info.pos2,
                info.dir2,
                String::from_utf8_lossy(&info.inserted_seq).to_string(),
                info.sv_id,
                info.checked_read_num_tumor,
                info.supporting_read_num_tumor,
                info.checked_read_num_control,
                info.supporting_read_num_control,
                info.is_filter,
                info.identical_svs,
                info.supporting_read_num_total,
                info.haplotype_bp1,
                info.haplotype_bp2,
                info.sv_type,
        )?;
    }

    for info in db1 {
        if merged_sv_id_1.contains(&info.sv_id) {
            continue;
        }
        writeln!(wfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                info.chrom1,
                info.pos1,
                info.dir1,
                info.chrom2,
                info.pos2,
                info.dir2,
                String::from_utf8_lossy(&info.inserted_seq).to_string(),
                "1_".to_string() + &info.sv_id,
                info.checked_read_num_tumor,
                info.supporting_read_num_tumor,
                info.checked_read_num_control,
                info.supporting_read_num_control,
                "1_".to_string() + &info.is_filter,
                info.identical_svs.split(",").map(|s| format!("1_{}", s)).collect::<Vec<String>>().join(","),
                info.supporting_read_num_total,
                "1_".to_string() + &info.haplotype_bp1,
                "1_".to_string() + &info.haplotype_bp2,
                info.sv_type,
        )?;
    }

    for info in db2 {
        if merged_sv_id_2.contains(&info.sv_id) {
            continue;
        }
        writeln!(wfile, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                info.chrom1,
                info.pos1,
                info.dir1,
                info.chrom2,
                info.pos2,
                info.dir2,
                String::from_utf8_lossy(&info.inserted_seq).to_string(),
                "2_".to_string() + &info.sv_id,
                info.checked_read_num_tumor,
                info.supporting_read_num_tumor,
                info.checked_read_num_control,
                info.supporting_read_num_control,
                "2_".to_string() + &info.is_filter,
                info.identical_svs.split(",").map(|s| format!("2_{}", s)).collect::<Vec<String>>().join(","),
                info.supporting_read_num_total,
                "2_".to_string() + &info.haplotype_bp1,
                "2_".to_string() + &info.haplotype_bp2,
                info.sv_type,
        )?;
    }

    Ok(())    
}

fn smith_waterman(seq1: &Vec<u8>, seq2: &Vec<u8>, identity_th: f64) -> Result<bool, Box<dyn Error>> {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -2i32 };
    // Gap open score: -5, gap extension score: -1
    let mut aligner = Aligner::new(-5, -1, &score);
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

    Ok(identity > identity_th)
}