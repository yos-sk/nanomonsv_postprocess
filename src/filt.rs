use rust_htslib::bam::{Header, IndexedReader, Read};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::io::BufRead;

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
    identical: String,
    total_support_read: usize,
    bp1_haplotype: String,
    bp2_haplotype: String,
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
        t_is_filter: String,
        t_identical: String,
        t_total_support_read: usize,
        t_bp1_haplotype: String,
        t_bp2_haplotype: String,
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
            identical: t_identical,
            total_support_read: t_total_support_read,
            bp1_haplotype: t_bp1_haplotype,
            bp2_haplotype: t_bp2_haplotype,
        }
    }
}

pub fn run(
    identical_list: &str,
    nanomonsv_result: &str,
    support_read_file: &str,
    bam_file: &str,
) -> Result<(), Box<dyn Error>> {
    let nanomonsv_reader =
        open_file(nanomonsv_result).expect(&format!("Could not open file {}", nanomonsv_result));
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

        let sv_info = NanomonsvInfo::new(
            bp1_contig,
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
            is_filter,
            "-".to_string(),
            support_read,
            "-".to_string(),
            "-".to_string(),
        );
        nanomonsv_db.insert(sv_id.clone(), sv_info);
    }

    eprintln!("Reading nanomonsv file finished.");

    let mut new_nanomonsv_result: Vec<NanomonsvInfo> = Vec::new();

    // for identical SV
    let identical_reader =
        open_file(identical_list).expect(&format!("Could not open file {}", identical_list));
    let mut identical_list: Vec<String> = Vec::new();

    for line in identical_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let key: Vec<String> = split_line[0].split(',').map(String::from).collect();
        let bp1_read: Vec<usize> = split_line[1]
            .split(',')
            .map(|s| s.parse::<usize>().unwrap())
            .collect();
        let bp1_read_sum = bp1_read.iter().sum();
        let bp2_read: Vec<usize> = split_line[2]
            .split(',')
            .map(|s| s.parse::<usize>().unwrap())
            .collect();
        let bp2_read_sum = bp2_read.iter().sum();
        if bp1_read_sum != bp2_read_sum {
            eprintln!(
                "{:?} has inconsistent support reads: {} {}",
                key, bp1_read_sum, bp2_read_sum
            );
        }
        let bp1_cat = split_line[3].to_string();
        let t_bp1_info: Vec<&str> = split_line[4].split(';').collect();
        let mut bp1_info: Vec<(String, usize, usize)> = Vec::new();
        for info in &t_bp1_info {
            if *info == "-" {
                bp1_info.push(("-".to_string(), 0, 0));
                break;
            }
            let split_info: Vec<&str> = info.split(",").collect();
            let t_sv_id: String = split_info[0].to_string();
            let bp: usize = split_info[1].parse::<usize>().unwrap();
            let n_read: usize = split_info[2].parse::<usize>().unwrap();
            bp1_info.push((t_sv_id, bp, n_read));
        }
        let bp2_cat = split_line[5].to_string();
        let t_bp2_info: Vec<&str> = split_line[6].split(';').collect();
        let mut bp2_info: Vec<(String, usize, usize)> = Vec::new();
        for info in &t_bp2_info {
            if *info == "-" {
                bp2_info.push(("-".to_string(), 0, 0));
                break;
            }
            let split_info: Vec<&str> = info.split(",").collect();
            let t_sv_id: String = split_info[0].to_string();
            let bp: usize = split_info[1].parse::<usize>().unwrap();
            let n_read: usize = split_info[2].parse::<usize>().unwrap();
            bp2_info.push((t_sv_id, bp, n_read));
        }
        //let mut identical: String = String::new();

        for sv_id in &key {
            identical_list.push(sv_id.to_string());
        }

        let mut c_sv_id: String = String::new();
        let mut max_n_read = 0;
        for info_1 in &bp1_info {
            if info_1.0 == "-" {
                continue;
            }
            for info_2 in &bp2_info {
                if info_1.0 == info_2.0 {
                    if c_sv_id.is_empty() {
                        c_sv_id = info_1.0.to_string(); 
                        max_n_read = info_1.2 + info_2.2;
                    } else {
                        if info_1.2 + info_2.2 > max_n_read {
                            c_sv_id = info_1.0.to_string();
                            max_n_read = info_1.2 + info_2.2;
                        }
                    }
                } else {
                    continue;
                }
            }
        }

        if c_sv_id.is_empty() {
            for info in &bp1_info {
                if info.2 > max_n_read {
                    c_sv_id = info.0.to_string();
                    max_n_read = info.2;
                } 
            }

            for info in &bp2_info {
                if info.2 > max_n_read {
                    c_sv_id = info.0.to_string();
                    max_n_read = info.2;
                }
            }
        }

        if !c_sv_id.is_empty() {
            let sv_info = match nanomonsv_db.get(&c_sv_id) {
                Some(value) => value.clone(),
                None => {
                    panic!("Error while reading nanomonsv result file: {}", c_sv_id);
                }
            };
            let new_sv_info = NanomonsvInfo::new(
                sv_info.bp1_contig.clone(),
                sv_info.bp1_pos,
                sv_info.bp1_strand.clone(),
                sv_info.bp2_contig.clone(),
                sv_info.bp2_pos,
                sv_info.bp2_strand.clone(),
                sv_info.insert_seq.clone(),
                sv_info.sv_id.clone(),
                sv_info.total_read,
                sv_info.support_read,
                sv_info.control_total_read,
                sv_info.control_support_read,
                "PASS".to_string(),
                key.join(","),
                bp1_read_sum,
                bp1_cat.clone(),
                bp2_cat.clone(),
            );
            new_nanomonsv_result.push(new_sv_info);
        } else {
            let mut max_sv_id = String::new();
            let mut s_read = 0;
            for sv_id in key.iter() {
                let info = match nanomonsv_db.get(sv_id) {
                    Some(value) => value.clone(),
                    None => {
                        panic!("Error while reading nanomonsv result file");
                    }
                };

                if info.support_read > s_read {
                    max_sv_id = sv_id.to_string();
                    s_read = info.support_read;
                }
            }

            let sv_info = match nanomonsv_db.get(&max_sv_id) {
                Some(value) => value.clone(),
                None => {
                    panic!("Error while reading nanomonsv result file");
                }
            };
            let new_sv_info = NanomonsvInfo::new(
                sv_info.bp1_contig.clone(),
                sv_info.bp1_pos,
                sv_info.bp1_strand.clone(),
                sv_info.bp2_contig.clone(),
                sv_info.bp2_pos,
                sv_info.bp2_strand.clone(),
                sv_info.insert_seq.clone(),
                sv_info.sv_id.clone()/,
                sv_info.total_read,
                sv_info.support_read,
                sv_info.control_total_read,
                sv_info.control_support_read,
                "PASS".to_string(),
                key.join(","),
                bp1_read_sum,
                bp1_cat.clone(),
                bp2_cat.clone(),
            );
            new_nanomonsv_result.push(new_sv_info);
        }
    }
    eprintln!("Reading identical SV file finished.");

    let reader =
        open_file(support_read_file).expect(&format!("Could not open {}", support_read_file));
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

    let mut bam =
        IndexedReader::from_path(&bam_file).expect(&format!("Could not open {}", bam_file));
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
        /*if !identical_list.contains(sv_id) {
            new_nanomonsv_result.push(info.clone());
        }*/

        let support_reads: Vec<String> = if let Some(value) = read_db.get(sv_id) {
            value.to_vec()
        } else {
            Vec::new()
        };

        if support_reads.is_empty() {
            continue;
        }

        let start = if info.bp1_pos > 100 {
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
        let mut counted_qname: HashSet<String> = HashSet::new();
        let mut hap1 = 0;
        let mut hap2 = 0;
        let mut unassign = 0;
        for record in bam.records() {
            let record = record?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if support_reads.contains(&qname) {
                if counted_qname.contains(&qname) {
                    continue;
                } else {
                    counted_qname.insert(qname.clone());
                }
                //let mut hap = 0;
                for aux in record.aux_iter() {
                    let (tag, value) = aux?;
                    if tag != b"HP" {
                        continue;
                    }
                    //if value == rust_htslib::bam::record::Aux::String("HP1") {
                    if value == rust_htslib::bam::record::Aux::U8(1) {
                        hap1 += 1;
                    // } else if value == rust_htslib::bam::record::Aux::String("HP2") {
                    } else if value == rust_htslib::bam::record::Aux::U8(2) {
                        hap2 += 1;
                    } else {
                        unassign += 1;
                    }
                }
            } else {
                continue;
            }
        }
        eprint!("{}\t{},{},{}", sv_id, unassign, hap1, hap2);

        let bp1_cat = if hap1 > hap2 {
            "haplotype1"
        } else if hap2 > hap1 {
            "haplotype2"
        } else {
            "Unassigned"
        };
        let bp1_read_sum = unassign + hap2 + hap1;

        let start = if info.bp2_pos > 100 {
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
        let mut counted_qname: HashSet<String> = HashSet::new();
        let mut hap1 = 0;
        let mut hap2 = 0;
        let mut unassign = 0;
        for record in bam.records() {
            let record = record?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if support_reads.contains(&qname) {
                if counted_qname.contains(&qname) {
                    continue;
                } else {
                    counted_qname.insert(qname.clone());
                }
                //let mut hap = 0;
                for aux in record.aux_iter() {
                    let (tag, value) = aux?;
                    if tag != b"HP" {
                        continue;
                    }
                    if value == rust_htslib::bam::record::Aux::U8(1) {
                        hap1 += 1;
                    } else if value == rust_htslib::bam::record::Aux::U8(2) {
                        hap2 += 1;
                    } else {
                        unassign += 1;
                    }
                }
            }
        }
        let bp2_cat = if hap1 > hap2 {
            "haplotype1"
        } else if hap2 > hap1 {
            "haplotype2"
        } else {
            "Unassigned"
        };
        let bp2_read_sum = unassign + hap2 + hap1;

        if bp1_read_sum != bp2_read_sum {
            eprintln!(
                "{} has inconsistent support reads: {} {}",
                sv_id, bp1_read_sum, bp2_read_sum
            );
        }
        if !identical_list.contains(sv_id) {
            eprintln!("\t{},{},{}\tsingle", unassign, hap1, hap2);
            let new_sv_info = NanomonsvInfo::new(
                info.bp1_contig.clone(),
                info.bp1_pos,
                info.bp1_strand.clone(),
                info.bp2_contig.clone(),
                info.bp2_pos,
                info.bp2_strand.clone(),
                info.insert_seq.clone(),
                info.sv_id.clone(),
                info.total_read,
                info.support_read,
                info.control_total_read,
                info.control_support_read,
                "PASS".to_string(),
                sv_id.clone(),
                hap1 + hap2 + unassign,
                bp1_cat.to_string(),
                bp2_cat.to_string(),
            );
            new_nanomonsv_result.push(new_sv_info.clone());
        } else {
            eprintln!("\t{},{},{}\tidentical", unassign, hap1, hap2);
        }
    }
    eprintln!("Counting SV haplotype finished.");

    // write new nanomonsv results
    println!("Chr_1\tPos_1\tDir_1\tChr_2\tPos_2\tDir_2\tInserted_Seq\tSV_ID\tChecked_Read_Num_Tumor\tSupporting_Read_Num_Tumor\tChecked_Read_Num_Control\tSupporting_Read_Num_Control\tIs_Filter\tIdentical_SVs\tSupporting_Read_Num_Total\tHaplotype_bp1\tHaplotype_bp2");
    for sv_info in new_nanomonsv_result.iter() {       
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            sv_info.bp1_contig,
            sv_info.bp1_pos,
            sv_info.bp1_strand,
            sv_info.bp2_contig,
            sv_info.bp2_pos,
            sv_info.bp2_strand,
            sv_info.insert_seq,
            sv_info.sv_id,
            sv_info.total_read,
            sv_info.support_read,
            sv_info.control_total_read,
            sv_info.control_support_read,
            sv_info.is_filter,
            sv_info.identical,
            sv_info.total_support_read,
            sv_info.bp1_haplotype,
            sv_info.bp2_haplotype,
        );
    }
    eprintln!("Recording new SV results finished.");
    Ok(())
}
