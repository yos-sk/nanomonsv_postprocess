use std::error::Error;
use std::io::BufRead;
use std::collections::HashMap;
use rust_htslib::bam::{IndexedReader, Read, Header};

use nanomonsv_postprocess::open_file;

#[derive(Eq, Hash, PartialEq, Clone)]
struct NanomonsvInfo {
    bp1_contig: String,
    bp1_pos: usize,
    bp1_strand: String,
    bp2_contig: String,
    bp2_pos: usize,
    bp2_strand: String,
    insert_seq: String,
    sv_id: String,
    total_read: usize,
    support_read: usize,
    control_total_read: usize,
    control_support_read: usize,
    is_filter: String,
}
impl NanomonsvInfo {
    fn new(
        t_bp1_c: String, 
        t_bp1_p: usize, 
        t_bp1_st: String, 
        t_bp2_c: String, 
        t_bp2_p: usize, 
        t_bp2_st: String, 
        t_insert_seq: String, 
        t_sv_id: String,
        t_total_read: usize, 
        t_support_read: usize, 
        t_ctl_total_read: usize, 
        t_ctl_support_read: usize, 
        t_is_filter: String
    ) -> Self {
        NanomonsvInfo {
            bp1_contig: t_bp1_c,
            bp1_pos: t_bp1_p,
            bp1_strand: t_bp1_st,
            bp2_contig: t_bp2_c,
            bp2_pos: t_bp2_p,
            bp2_strand: t_bp2_st,
            insert_seq: t_insert_seq,
            sv_id: t_sv_id,
            total_read: t_total_read,
            support_read: t_support_read,
            control_total_read: t_ctl_total_read,
            control_support_read: t_ctl_support_read,
            is_filter: t_is_filter,
        }
    }
}
 
pub fn run(identical_list: &str, nanomonsv_result: &str, support_read_file: &str, bam_file: &str) -> Result<(), Box<dyn Error>> {
    let nanomonsv_reader = open_file(nanomonsv_result).expect(&format!("Could not open file {}", nanomonsv_result));
    let mut nanomonsv_db: HashMap<String, NanomonsvInfo> = HashMap::new();

    for (i, line) in nanomonsv_reader.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let bp1_contig = split_line[0].to_string();
        let bp1_position = split_line[1].parse::<usize>().unwrap();
        let bp1_strand = split_line[2].to_string();
        let bp2_contig = split_line[3].to_string();
        let bp2_position = split_line[4].parse::<usize>().unwrap();
        let bp2_strand = split_line[5].to_string();
        let insert_seq = split_line[6].to_string();
        let sv_id = split_line[7].to_string();
        let total_read = split_line[8].parse::<usize>().unwrap();
        let support_read = split_line[9].parse::<usize>().unwrap();
        let control_total_read = split_line[10].parse::<usize>().unwrap();
        let control_support_read = split_line[11].parse::<usize>().unwrap();
        let is_filter = split_line[12].to_string();
        
        let sv_info = NanomonsvInfo::new(bp1_contig, 
                                         bp1_position, 
                                         bp1_strand, 
                                         bp2_contig, 
                                         bp2_position, 
                                         bp2_strand, 
                                         insert_seq, 
                                         sv_id.clone(),
                                         total_read, 
                                         support_read, 
                                         control_total_read, 
                                         control_support_read, 
                                         is_filter);
        nanomonsv_db.insert(sv_id.clone(), sv_info);
    }

    eprintln!("Reading nanomonsv file finished.");

    let mut new_nanomonsv_result: Vec<NanomonsvInfo> = Vec::new();

    // for identical SV
    let identical_reader = open_file(identical_list).expect(&format!("Could not open file {}", identical_list));
    let mut identical_list: Vec<String> = Vec::new();

    for line in identical_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let key: Vec<String> = split_line[0].split(',').map(String::from).collect();
        let bp1_cat = split_line[3].to_string();
        let bp1_info: Vec<&str> = split_line[4].split(',').collect();
        let bp2_cat = split_line[5].to_string();
        let bp2_info: Vec<&str> = split_line[6].split(',').collect();
        
        for sv_id in key.iter() {
            identical_list.push(sv_id.to_string());
        }

        if bp1_cat != "Unassigned" && bp1_cat != "Ambiguous" && bp2_cat != "Unassigned" && bp2_cat != "Ambiguous" {
            let bp1_sv_id = bp1_info[0].to_string();
            let bp1_bp_num = bp1_info[1].parse::<usize>().unwrap();
            let bp2_sv_id = bp2_info[0].to_string();
            let bp2_bp_num = bp2_info[1].parse::<usize>().unwrap();

            let bp1_info = match nanomonsv_db.get(&bp1_sv_id) {
                Some(value) => value.clone(),
                None => {
                    panic!("Error while reading nanomonsv result file");
                },
            };
            let bp2_info = match nanomonsv_db.get(&bp2_sv_id) {
                Some(value) => value.clone(),
                None => {
                    panic!("Error while reading nanomonsv result file");
                },
            };

            if bp1_bp_num == 1 && bp2_bp_num == 1 { 
                let new_sv_info =  NanomonsvInfo::new(bp1_info.bp1_contig.clone(), 
                                                      bp1_info.bp1_pos, 
                                                      bp1_info.bp1_strand.clone(), 
                                                      bp2_info.bp1_contig.clone(), 
                                                      bp2_info.bp1_pos, 
                                                      bp2_info.bp1_strand.clone(), 
                                                      bp1_info.insert_seq.clone(), 
                                                      bp1_info.sv_id.clone() + "_id",
                                                      bp1_info.total_read, 
                                                      bp1_info.support_read, 
                                                      bp1_info.control_total_read, 
                                                      bp1_info.control_support_read, 
                                                      "PASS".to_string());
                new_nanomonsv_result.push(new_sv_info);
            } else if bp1_bp_num == 1 && bp2_bp_num == 2 {
                let new_sv_info =  NanomonsvInfo::new(bp1_info.bp1_contig.clone(), 
                                                      bp1_info.bp1_pos, 
                                                      bp1_info.bp1_strand.clone(), 
                                                      bp2_info.bp2_contig.clone(), 
                                                      bp2_info.bp2_pos, 
                                                      bp2_info.bp2_strand.clone(), 
                                                      bp1_info.insert_seq.clone(), 
                                                      bp1_info.sv_id.clone() + "_id",
                                                      bp1_info.total_read, 
                                                      bp1_info.support_read, 
                                                      bp1_info.control_total_read, 
                                                      bp1_info.control_support_read, 
                                                      "PASS".to_string());
                new_nanomonsv_result.push(new_sv_info);
            } else if bp1_bp_num == 2 && bp2_bp_num == 1 {
                let new_sv_info =  NanomonsvInfo::new(bp1_info.bp2_contig.clone(), 
                                                      bp1_info.bp2_pos, 
                                                      bp1_info.bp2_strand.clone(), 
                                                      bp2_info.bp1_contig.clone(), 
                                                      bp2_info.bp1_pos, 
                                                      bp2_info.bp1_strand.clone(), 
                                                      bp1_info.insert_seq.clone(), 
                                                      bp1_info.sv_id.clone() + "_id",
                                                      bp1_info.total_read, 
                                                      bp1_info.support_read, 
                                                      bp1_info.control_total_read, 
                                                      bp1_info.control_support_read, 
                                                      "PASS".to_string());
                new_nanomonsv_result.push(new_sv_info);
            } else {
                let new_sv_info =  NanomonsvInfo::new(bp1_info.bp2_contig.clone(), 
                                                      bp1_info.bp2_pos, 
                                                      bp1_info.bp2_strand.clone(), 
                                                      bp2_info.bp2_contig.clone(), 
                                                      bp2_info.bp2_pos, 
                                                      bp2_info.bp2_strand.clone(), 
                                                      bp1_info.insert_seq.clone(), 
                                                      bp1_info.sv_id.clone() + "_id",
                                                      bp1_info.total_read, 
                                                      bp1_info.support_read, 
                                                      bp1_info.control_total_read, 
                                                      bp1_info.control_support_read, 
                                                      "PASS".to_string());
                new_nanomonsv_result.push(new_sv_info);
            }
        } else {
            let mut max_sv_id = String::new();
            let mut n_read = 0;
            for sv_id in key.iter() {
                let info = match nanomonsv_db.get(sv_id) {
                    Some(value) => value.clone(),
                    None => {
                        panic!("Error while reading nanomonsv result file");
                    },
                };

                if info.support_read > n_read {
                    max_sv_id = sv_id.to_string();
                    n_read = info.support_read;
                }
            }

            let sv_info = match nanomonsv_db.get(&max_sv_id) {
                Some(value) => value.clone(),
                None => {
                    panic!("Error while reading nanomonsv result file");
                },
            };
            new_nanomonsv_result.push(sv_info);
        }    
    }
    eprintln!("Reading identical SV file finished.");

    let reader = open_file(support_read_file).expect(&format!("Could not open {}", support_read_file));
    let mut read_db: HashMap<String, Vec<String>> = HashMap::new();
    // collect support reads of SVs
    for line in reader.lines() {
        // eprintln!("line: {}", i);
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let sv_id = split_line[7].to_string();
        let read_id = split_line[8].to_string();
        if let Some(value) = read_db.get_mut(&sv_id) {
            value.push(read_id.clone());
        } else {
            read_db.insert(sv_id.clone(), vec![read_id.clone()]);
        }
    }

    eprintln!("Reading support read file finished.");

    let mut bam = IndexedReader::from_path(&bam_file).expect(&format!("Could not open {}", bam_file));
    // let mut hap_db: HashMap<String, HashMap<usize, Vec<usize>>> = HashMap::new();
    let header = Header::from_template(bam.header());
    let mut reference_sizes: HashMap<String, usize> = HashMap::new();

    // Extract contig size
    for (key, records) in header.to_hashmap() {
        for record in records {
            if key != "SQ" {
                continue;
            } 
            let size = record["LN"].parse::<usize>().unwrap();
            reference_sizes.insert(record["SN"].clone(), size);
        }
    }

    // count hap1, hap2, unassigned reads
    for (sv_id, info) in nanomonsv_db.iter() {
        if !identical_list.contains(sv_id) {
            new_nanomonsv_result.push(info.clone());        
        }
        
        let support_reads: Vec<String> = if let Some(value) = read_db.get(sv_id) {
            value.to_vec()
        } else {
            Vec::new()
        };

        if support_reads.is_empty() {
            continue;
        }

        let start = if info.bp1_pos >  100 {
            info.bp1_pos - 100
        } else {
            0
        };
        let max_size = if let Some(length) = reference_sizes.get(&info.bp1_contig) {
            *length
        } else {
            0
        };
        let end = if info.bp1_pos + 100 > max_size {
            max_size
        } else {
            info.bp1_pos + 100
        };
        let _ = bam.fetch((&info.bp1_contig, start as u64, end as u64));
        let mut hap1 = 0;
        let mut hap2 = 0;
        let mut unassign = 0;
        for record in bam.records() {
            let record = record?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if support_reads.contains(&qname) {
                //let mut hap = 0;
                for aux in record.aux_iter() {
                    let (tag, value) = aux?;
                    if tag != b"HP" {
                        continue;
                    }
                    if value == rust_htslib::bam::record::Aux::String("HP1") {
                    // if value == rust_htslib::bam::record::Aux::U8(1) {
                        hap1 += 1;
                    } else if value == rust_htslib::bam::record::Aux::String("HP2") {
                    // else if value == rust_htslib::bam::record::Aux::U8(2) {
                        hap2 += 2;
                    } else {
                        unassign += 1;
                    }
                }
            } else {
                continue;
            }
        }
        eprint!("{}\t{},{},{}", sv_id, unassign, hap1, hap2);

        let start = if info.bp2_pos >  100 {
            info.bp2_pos - 100
        } else {
            0
        };
        let max_size = if let Some(length) = reference_sizes.get(&info.bp2_contig) {
            *length
        } else {
            0
        };
        let end = if info.bp2_pos + 100 > max_size {
            max_size
        } else {
            info.bp2_pos + 100
        };
        let _ = bam.fetch((&info.bp2_contig, start as u64, end as u64));
        let mut hap1 = 0;
        let mut hap2 = 0;
        let mut unassign = 0;
        for record in bam.records() {
            let record = record?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if support_reads.contains(&qname) {
                //let mut hap = 0;
                for aux in record.aux_iter() {
                    let (tag, value) = aux?;
                    if tag != b"HP" {
                        continue;
                    }
                    if value == rust_htslib::bam::record::Aux::String("HP1") {
                        hap1 += 1;
                    } else if value == rust_htslib::bam::record::Aux::String("HP2") {
                        hap2 += 1;
                    } else {
                        unassign += 1;
                    }
                }
            }
        }
        if !identical_list.contains(sv_id) {
            eprintln!("\t{},{},{}\tsingle", unassign, hap1, hap2);
        } else {
            eprintln!("\t{},{},{}\tidentical", unassign, hap1, hap2);
        }
    }
    eprintln!("Counting SV haplotype finished.");
    
    // write new nanomonsv results
    for result in new_nanomonsv_result.iter() {
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", result.bp1_contig, result.bp1_pos, result.bp1_strand, result.bp2_contig, result.bp2_pos, result.bp2_strand, result.insert_seq, result.sv_id, result.total_read, result.support_read, result.control_total_read, result.control_support_read, result.is_filter);
    }

    eprintln!("Recording new SV results finished.");
    Ok(())
}