use rust_htslib::bam::{Header, IndexedReader, Read};
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;
use std::io::BufRead;

use nanomonsv_postprocess::open_file;
#[path = "./smith_waterman.rs"]
mod smith_waterman;

#[derive(Eq, Hash, PartialEq, Clone)]
struct SVInfo {
    sv_id: String,
    bp1_contig: String,
    bp1_pos: usize,
    bp1_seq: Vec<u8>,
    bp2_contig: String,
    bp2_pos: usize,
    bp2_seq: Vec<u8>,
}
impl SVInfo {
    fn new() -> Self {
        SVInfo {
            sv_id: String::new(),
            bp1_contig: String::new(),
            bp1_pos: 0,
            bp1_seq: Vec::new(),
            bp2_contig: String::new(),
            bp2_pos: 0,
            bp2_seq: Vec::new(),
        }
    }

    fn add_bp1_info(&mut self, sv_id: String, bp1_contig: String, bp1_pos: usize, bp1_seq: &[u8]) {
        self.sv_id = sv_id;
        self.bp1_contig = bp1_contig;
        self.bp1_pos = bp1_pos;
        self.bp1_seq = bp1_seq.to_vec();
    }

    fn add_bp2_info(&mut self, sv_id: String, bp2_contig: String, bp2_pos: usize, bp2_seq: &[u8]) {
        if self.sv_id != sv_id {
            eprintln!("Inconsistent SV information: {} {}", self.sv_id, sv_id);
        }
        self.bp2_contig = bp2_contig;
        self.bp2_pos = bp2_pos;
        self.bp2_seq = bp2_seq.to_vec();
    }
}

#[derive(Eq, PartialEq, Clone)]
struct IdenticalInfo {
    identical_bp1: HashSet<(String, usize)>,
    identical_bp2: HashSet<(String, usize)>,
}
impl IdenticalInfo {
    fn new() -> Self {
        IdenticalInfo {
            identical_bp1: HashSet::new(),
            identical_bp2: HashSet::new(),
        }
    }

    fn add_identical_bp1(&mut self, sv_id: String, bp_num: usize) {
        self.identical_bp1.insert((sv_id, bp_num));
    }

    fn add_identical_bp2(&mut self, sv_id: String, bp_num: usize) {
        self.identical_bp2.insert((sv_id, bp_num));
    }
}

pub fn run(
    input_bed: &str,
    support_read_file: &str,
    bam_file: &str,
    min_identity: f64,
    min_length: usize,
) -> Result<(), Box<dyn Error>> {
    // open bed file including SV breakpoint information with ±100 bp sequence
    let bp_reader = open_file(input_bed).expect(&format!("Could not open file {}", input_bed));
    let mut bp_info_db: HashSet<SVInfo> = HashSet::new();

    // collect SV breakpoint information
    let mut line1_info = String::new();
    for (i, line) in bp_reader.lines().enumerate() {
        let line = line?;

        match i % 2 {
            0 => {
                line1_info = line;
            }
            1 => {
                let mut info = SVInfo::new();
                let split_line_1: Vec<&str> = line1_info.split('\t').collect();
                let contig_1 = split_line_1[0].to_string();
                let position_1 = split_line_1[1].parse::<usize>().unwrap();
                let sv_id_1 = split_line_1[3].to_string();
                let seq_1 = split_line_1[4].as_bytes();
                info.add_bp1_info(sv_id_1, contig_1, position_1, &seq_1);
                let split_line_2: Vec<&str> = line.split('\t').collect();
                let contig_2 = split_line_2[0].to_string();
                let position_2 = split_line_2[1].parse::<usize>().unwrap();
                let sv_id_2 = split_line_2[3].to_string();
                let seq_2 = split_line_2[4].as_bytes();
                info.add_bp2_info(sv_id_2, contig_2, position_2, &seq_2);
                bp_info_db.insert(info.clone());
            }
            _ => (),
        }
    }

    // grouping identical SVs by smith-waterman algorithm
    //let threshold = 99.0;
    let mut identical_pairs: Vec<(String, String, usize)> = Vec::new();
    // realignment of 2 * 2 breakpoint combination
    for sv_info_1 in &bp_info_db {
        let sv_type1 = &sv_info_1.sv_id[0..1];
        for sv_info_2 in &bp_info_db {
            if sv_info_1.sv_id >= sv_info_2.sv_id {
                continue;
            }
            let sv_type2 = &sv_info_2.sv_id[0..1];

            if sv_type1 == "d" && sv_type2 == "i" {
                continue;
            }
            // smith_waterman algorithm
            let result1 = match smith_waterman::run(
                &sv_info_1.bp1_seq,
                &sv_info_2.bp1_seq,
                min_identity,
                min_length,
            ) {
                Ok(flag) => flag,
                _ => false,
            };
            let result2 = match smith_waterman::run(
                &sv_info_1.bp1_seq,
                &sv_info_2.bp2_seq,
                min_identity,
                min_length,
            ) {
                Ok(flag) => flag,
                _ => false,
            };

            if result1 && !result2 {
                let pattern = 0;
                let result = match smith_waterman::run(
                    &sv_info_1.bp2_seq,
                    &sv_info_2.bp2_seq,
                    min_identity,
                    min_length,
                ) {
                    Ok(flag) => flag,
                    _ => false,
                };
                if result {
                    eprintln!("{}\t{}\t{}", sv_info_1.sv_id, sv_info_2.sv_id, pattern);
                    // update breakpoint information
                    identical_pairs.push((
                        sv_info_1.sv_id.clone(),
                        sv_info_2.sv_id.clone(),
                        pattern,
                    ));
                }
            } else {
                let pattern = 1;
                let result = match smith_waterman::run(
                    &sv_info_1.bp2_seq,
                    &sv_info_2.bp1_seq,
                    min_identity,
                    min_length,
                ) {
                    Ok(flag) => flag,
                    _ => false,
                };
                if result {
                    eprintln!("{}\t{}\t{}", sv_info_1.sv_id, sv_info_2.sv_id, pattern);
                    // update breakpoint information
                    identical_pairs.push((
                        sv_info_1.sv_id.clone(),
                        sv_info_2.sv_id.clone(),
                        pattern,
                    ));
                }
            }
        }
    }

    // grouping
    let mut groups: Vec<Vec<(String, String, usize)>> = Vec::new();
    for pair in identical_pairs {
        let mut found = false;
        for group in &mut groups {
            if group
                .iter()
                .any(|x| x.0 == pair.0 || x.1 == pair.0 || x.0 == pair.1 || x.1 == pair.1)
            {
                group.push(pair.clone());
                found = true;
                break;
            }
        }
        if !found {
            groups.push(vec![pair.clone()]);
        }
    }

    let mut group_db: HashMap<Vec<String>, IdenticalInfo> = HashMap::new();

    for pairs in groups.iter_mut() {
        let mut info = IdenticalInfo::new();
        let mut group_db_key: Vec<String> = Vec::new();
        pairs.sort_by_key(|pair| pair.0.clone());
        // eprintln!("{:?}", pairs);
        for pair in pairs {
            match pair.2 {
                0 => {
                    if info.identical_bp2.contains(&(pair.0.clone(), 1))
                        || info.identical_bp2.contains(&(pair.1.clone(), 1))
                    {
                        info.add_identical_bp1(pair.0.clone(), 2);
                        info.add_identical_bp2(pair.0.clone(), 1);
                        info.add_identical_bp1(pair.1.clone(), 2);
                        info.add_identical_bp2(pair.1.clone(), 1);
                    } else {
                        info.add_identical_bp1(pair.0.clone(), 1);
                        info.add_identical_bp2(pair.0.clone(), 2);
                        info.add_identical_bp1(pair.1.clone(), 1);
                        info.add_identical_bp2(pair.1.clone(), 2);
                    }
                }
                1 => {
                    if info.identical_bp2.contains(&(pair.0.clone(), 1))
                        || info.identical_bp2.contains(&(pair.1.clone(), 2))
                    {
                        info.add_identical_bp1(pair.0.clone(), 2);
                        info.add_identical_bp2(pair.0.clone(), 1);
                        info.add_identical_bp1(pair.1.clone(), 1);
                        info.add_identical_bp2(pair.1.clone(), 2);
                    } else {
                        info.add_identical_bp1(pair.0.clone(), 1);
                        info.add_identical_bp2(pair.0.clone(), 2);
                        info.add_identical_bp1(pair.1.clone(), 2);
                        info.add_identical_bp2(pair.1.clone(), 1);
                    }
                }
                _ => (),
            }
            if !group_db_key.contains(&pair.0) {
                group_db_key.push(pair.0.clone());
            }
            if !group_db_key.contains(&pair.1) {
                group_db_key.push(pair.1.clone());
            }
        }
        group_db.insert(group_db_key.clone(), info.clone());
    }

    let _ = classify_haplotype_fetch_bam(&group_db, &bp_info_db, support_read_file, bam_file);
    Ok(())
}

fn classify_haplotype_fetch_bam(
    group_db: &HashMap<Vec<String>, IdenticalInfo>,
    bp_info_db: &HashSet<SVInfo>,
    support_read_file: &str,
    bam_file: &str,
) -> Result<(), Box<dyn Error>> {
    let reader =
        open_file(support_read_file).expect(&format!("Could not open {}", support_read_file));
    let mut read_db: HashMap<String, Vec<String>> = HashMap::new();
    // collect support reads of SVs
    for line in reader.lines() {
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

    let mut bam =
        IndexedReader::from_path(&bam_file).expect(&format!("Could not open {}", bam_file));
    let mut hap_db: HashMap<String, HashMap<usize, Vec<usize>>> = HashMap::new();
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

    for info in bp_info_db.iter() {
        let sv_id = &info.sv_id;
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
        for record in bam.records() {
            let record = record?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if support_reads.contains(&qname) {
                // inhibit duplicate count
                if counted_qname.contains(&qname) {
                    continue;
                } else {
                    counted_qname.insert(qname.clone());
                }
                let mut hap = 0;
                for aux in record.aux_iter() {
                    let (tag, value) = aux?;
                    if tag != b"HP" {
                        continue;
                    }
                    // if value == rust_htslib::bam::record::Aux::String("HP1") {
                    if value == rust_htslib::bam::record::Aux::U8(1) {
                        hap = 1;
                    //v} else if value == rust_htslib::bam::record::Aux::String("HP2") {
                    } else if value == rust_htslib::bam::record::Aux::U8(2) {
                        hap = 2;
                    } else {
                        hap = 0;
                    }
                }
                if let Some(bp_map) = hap_db.get_mut(sv_id) {
                    if let Some(hap_cnt_vec) = bp_map.get_mut(&1) {
                        hap_cnt_vec[hap] += 1;
                    } else {
                        let new_hap_cnt_vec = if hap == 0 {
                            vec![1, 0, 0]
                        } else if hap == 1 {
                            vec![0, 1, 0]
                        } else {
                            vec![0, 0, 1]
                        };
                        bp_map.insert(1, new_hap_cnt_vec);
                    }
                } else {
                    let new_hap_cnt_vec = if hap == 0 {
                        vec![1, 0, 0]
                    } else if hap == 1 {
                        vec![0, 1, 0]
                    } else {
                        vec![0, 0, 1]
                    };
                    let tmp_bp_map: HashMap<usize, Vec<usize>> =
                        HashMap::from([(1, new_hap_cnt_vec)]);
                    hap_db.insert(sv_id.clone(), tmp_bp_map);
                }
            } else {
                continue;
            }
        }
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
        for record in bam.records() {
            let record = record?;
            let qname = String::from_utf8(record.qname().to_vec())?;
            if support_reads.contains(&qname) {
                // inhibit duplicate count
                if counted_qname.contains(&qname) {
                    continue;
                } else {
                    counted_qname.insert(qname.clone());
                }
                let mut hap = 0;
                for aux in record.aux_iter() {
                    let (tag, value) = aux?;
                    if tag != b"HP" {
                        continue;
                    }
                    if value == rust_htslib::bam::record::Aux::U8(1) {
                        hap = 1;
                    } else if value == rust_htslib::bam::record::Aux::U8(2) {
                        hap = 2;
                    } else {
                        hap = 0;
                    }
                }
                if let Some(bp_map) = hap_db.get_mut(sv_id) {
                    if let Some(hap_cnt_vec) = bp_map.get_mut(&2) {
                        hap_cnt_vec[hap] += 1;
                    } else {
                        let new_hap_cnt_vec = if hap == 0 {
                            vec![1, 0, 0]
                        } else if hap == 1 {
                            vec![0, 1, 0]
                        } else {
                            vec![0, 0, 1]
                        };
                        bp_map.insert(2, new_hap_cnt_vec);
                    }
                } else {
                    let new_hap_cnt_vec = if hap == 0 {
                        vec![1, 0, 0]
                    } else if hap == 1 {
                        vec![0, 1, 0]
                    } else {
                        vec![0, 0, 1]
                    };
                    let tmp_bp_map: HashMap<usize, Vec<usize>> =
                        HashMap::from([(2, new_hap_cnt_vec)]);
                    hap_db.insert(sv_id.clone(), tmp_bp_map);
                }
            } else {
                continue;
            }
        }
    }
    // using group_db and hap_db, determine which breakpoint is true
    let mut assigned = 0;
    let mut ambiguous = 0;
    let mut unassigned = 0;
    for (key, value) in group_db.iter() {
        let mut sv_bp1_cnt = vec![0, 0, 0];
        let mut max_sv1_id = vec![(String::new(), 0); 3];
        let mut max_sv1_cnt = vec![0, 0, 0];
        // iterate by value, then count haplotype
        // breakpoint 1
        for item in value.identical_bp1.iter() {
            let sv_id = (&item.0).to_string();
            let bp_num = item.1;
            if let Some(cnt_map) = hap_db.get(&sv_id) {
                if let Some(cnt_vec) = cnt_map.get(&bp_num) {
                    sv_bp1_cnt[0] += cnt_vec[0];
                    sv_bp1_cnt[1] += cnt_vec[1];
                    sv_bp1_cnt[2] += cnt_vec[2];
                    if cnt_vec[0] > max_sv1_cnt[0] {
                        max_sv1_cnt[0] = cnt_vec[0];
                        max_sv1_id[0] = (sv_id.clone(), bp_num);
                    }
                    if cnt_vec[1] > max_sv1_cnt[1] {
                        max_sv1_cnt[1] = cnt_vec[1];
                        max_sv1_id[1] = (sv_id.clone(), bp_num);
                    }
                    if cnt_vec[2] > max_sv1_cnt[2] {
                        max_sv1_cnt[2] = cnt_vec[2];
                        max_sv1_id[2] = (sv_id.clone(), bp_num);
                    }
                }
            }
        }

        let mut sv_bp2_cnt = vec![0, 0, 0];
        let mut max_sv2_id = vec![(String::new(), 0); 3];
        let mut max_sv2_cnt = vec![0, 0, 0];
        // breakpoint 2
        for item in value.identical_bp2.iter() {
            let sv_id = (&item.0).to_string();
            let bp_num = item.1;
            if let Some(cnt_map) = hap_db.get(&sv_id) {
                if let Some(cnt_vec) = cnt_map.get(&bp_num) {
                    sv_bp2_cnt[0] += cnt_vec[0];
                    sv_bp2_cnt[1] += cnt_vec[1];
                    sv_bp2_cnt[2] += cnt_vec[2];
                    if cnt_vec[0] > max_sv2_cnt[0] {
                        max_sv2_cnt[0] = cnt_vec[0];
                        max_sv2_id[0] = (sv_id.clone(), bp_num);
                    }
                    if cnt_vec[1] > max_sv2_cnt[1] {
                        max_sv2_cnt[1] = cnt_vec[1];
                        max_sv2_id[1] = (sv_id.clone(), bp_num);
                    }
                    if cnt_vec[2] > max_sv2_cnt[2] {
                        max_sv2_cnt[2] = cnt_vec[2];
                        max_sv2_id[2] = (sv_id.clone(), bp_num);
                    }
                }
            }
        }

        eprintln!(
            "{:?}\t{:?}\t{:?}\t{:?}\t{:?}",
            key, value.identical_bp1, value.identical_bp2, sv_bp1_cnt, sv_bp2_cnt
        );
        print!(
            "{}\t{}\t{}",
            key.join(","),
            sv_bp1_cnt
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<String>>()
                .join(","),
            sv_bp2_cnt
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<String>>()
                .join(",")
        );

        if sv_bp1_cnt[1] == 0 && sv_bp1_cnt[2] == 0 {
            unassigned += 1;
            print!("\tUnassigned\t-");
        } else if sv_bp1_cnt[1] > 0 && sv_bp1_cnt[2] > 0 {
            ambiguous += 1;
            print!("\tAmbiguous\t-");
        } else {
            assigned += 1;
            if sv_bp1_cnt[1] > 0 {
                print!("\thaplotype1\t{},{}", max_sv1_id[1].0, max_sv1_id[1].1);
            } else {
                print!("\thaplotype2\t{},{}", max_sv1_id[2].0, max_sv1_id[2].1);
            }
        }

        if sv_bp2_cnt[1] == 0 && sv_bp2_cnt[2] == 0 {
            unassigned += 1;
            println!("\tUnassigned\t-");
        } else if sv_bp2_cnt[1] > 0 && sv_bp2_cnt[2] > 0 {
            ambiguous += 1;
            println!("\tAmbiguous\t-");
        } else {
            assigned += 1;
            if sv_bp2_cnt[1] > 0 {
                println!("\thaplotype1\t{},{}", max_sv2_id[1].0, max_sv2_id[1].1);
            } else {
                println!("\thaplotype2\t{},{}", max_sv2_id[2].0, max_sv2_id[2].1);
            }
        }
    }

    eprintln!("Assigned: {}", assigned);
    eprintln!("Unassigned: {}", unassigned);
    eprintln!("Ambiguous: {}", ambiguous);
    Ok(())
}
