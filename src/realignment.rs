use std::error::Error;
use std::io::BufRead;
use std::collections::HashMap;
use std::collections::HashSet;

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


pub fn run(bp_file: &str, support_read_file: &str, nanomonsv_bp_file: &str) -> Result<(), Box<dyn Error>> {
    // open bed file including SV breakpoint information with ±100 bp sequence
    let bp_reader = open_file(bp_file).expect(&format!("Could not open file {}", bp_file));
    let mut bp_info_db: HashSet<SVInfo> = HashSet::new();

    // collect SV breakpoint information
    let mut line1_info = String::new();
    for (i, line) in bp_reader.lines().enumerate() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let contig = split_line[0].to_string();
        let position = split_line[1].parse::<usize>().unwrap();
        let sv_id = split_line[3].to_string();
        let seq = split_line[4].as_bytes();

        match i % 2 {
            0 => {
                line1_info = line;
            },
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
            },
            _ => (),
        }    
    }

    // grouping identical SVs by smith-waterman algorithm
    let threshold = 99.0;
    let mut identical_pairs: Vec<(String, String, usize)> = Vec::new();
    // realignment of 2 * 2 breakpoint combination
    for sv_info_1 in &bp_info_db {
        let sv_type1 = &sv_info_1.sv_id[0..1];
        for sv_info_2 in &bp_info_db {
            if sv_info_1.sv_id >= sv_info_2.sv_id {
                continue;
            }
            let sv_type2 = &sv_info_2.sv_id[0..1];
            if sv_type1 != sv_type2 {
                continue;
            }
            // smith_waterman algorithm
            let result1 = match smith_waterman::run(&sv_info_1.bp1_seq, &sv_info_2.bp1_seq) {
                Ok(id) => id,
                _ => -1.0,
            };
            let result2 = match smith_waterman::run(&sv_info_1.bp1_seq, &sv_info_2.bp2_seq) {
                Ok(id) => id,
                _ => -1.0,
            };
            let mut bp1_max_id = 0.0;
            let mut bp2_max_id = 0.0;
            let mut pattern = 0;
            if result1 > result2 {
                bp1_max_id = result1;
                pattern = 0;
                match smith_waterman::run(&sv_info_1.bp2_seq, &sv_info_2.bp2_seq) {
                    Ok(id) => bp2_max_id = id,
                    _ => bp2_max_id = -1.0,
                }
            } else {
                bp1_max_id = result2;
                pattern = 1;
                match smith_waterman::run(&sv_info_1.bp2_seq, &sv_info_2.bp1_seq) {           
                    Ok(id) => bp2_max_id = id,
                    _ => bp2_max_id = -1.0,
                }
            }

            if bp1_max_id >= threshold && bp2_max_id >= threshold {
                eprintln!("{}\t{}", sv_info_1.sv_id, sv_info_2.sv_id);
                // update breakpoint information
                identical_pairs.push((sv_info_1.sv_id.clone(), sv_info_2.sv_id.clone(), pattern))

            }
        }
    }

    // grouping
    let mut groups: Vec<Vec<(String, String, usize)>> = Vec::new();
    for pair in identical_pairs {
        let mut found = false;
        for group in &mut groups {
            if group.iter().any(|x| x.0 == pair.0 || x.1 == pair.0 || x.0 == pair.1 || x.1 == pair.1) {
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

    for pairs in groups.iter() {
        let mut info = IdenticalInfo::new();
        let mut group_db_key: Vec<String> = Vec::new();
        // TODO: asserion of inconsistent case. e.g) (r3, r4, 0), (r3, r5, 0), (r4, r5, 1)
        for pair in pairs {
            match pair.2 {
                0 => {
                    info.add_identical_bp1(pair.0.clone(), 1);
                    info.add_identical_bp2(pair.0.clone(), 2);
                    info.add_identical_bp1(pair.1.clone(), 1);
                    info.add_identical_bp2(pair.1.clone(), 2);
                },
                1 => {
                    info.add_identical_bp1(pair.0.clone(), 1);
                    info.add_identical_bp2(pair.0.clone(), 2);
                    info.add_identical_bp1(pair.1.clone(), 2);
                    info.add_identical_bp2(pair.1.clone(), 1);
                },
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

    /*
    for (key, value) in group_db.iter() {
        println!("{:?} {:?} {:?}", key, value.identical_bp1, value.identical_bp2);
    }
    */
    
    
    let _ = classify_haplotype(&group_db, &bp_info_db, support_read_file, nanomonsv_bp_file);
    Ok(())
}


fn classify_haplotype(group_db: &HashMap<Vec<String>, IdenticalInfo>, bp_info_db: &HashSet<SVInfo>, support_read_file: &str, nanomonsv_bp_file: &str) -> Result<(), Box<dyn Error>> {
    let reader = open_file(support_read_file).expect(&format!("Could not open {}", support_read_file));
    let mut read_db: HashMap<String, String> = HashMap::new();
    // collect support reads of SVs
    for line in reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let sv_id = split_line[7].to_string();
        let read_id = split_line[8].to_string();
        read_db.insert(read_id, sv_id);
    }

    // count haplotypes of each break point 
    let bp_reader = open_file(nanomonsv_bp_file).expect(&format!("Could not open {}", nanomonsv_bp_file));
    let mut hap_db: HashMap<String, HashMap<usize, Vec<usize>>> = HashMap::new();

    let read_set: HashSet<String> = HashSet::new();

    for line in bp_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let contig = split_line[0].parse::<usize>().unwrap();
        let position = split_line[2].parse::<usize>().unwrap();
        let read_id = split_line[3].to_string();
        // search sv_id
        let sv_id = match read_db.get(&read_id) {
            Some(value) => value.to_string(),
            None => {
                continue;
            },
        };
        // search SVInfo
        let mut t_sv_info = SVInfo::new();
        for info in bp_info_db.iter() {
            if info.sv_id == sv_id {
                t_sv_info = info.clone();
                break;
            }
        }

        // determine breakpoint
        let bp_num = if t_sv_info.bp1_pos - 50 <= position && position <= t_sv_info.bp1_pos + 50 {
            1
        } else if t_sv_info.bp2_pos - 50 <= position && position <= t_sv_info.bp2_pos + 50 {
            2
        } else {
            0
        };

        if bp_num == 0 {
            continue;
        }

        // determine haplotype of breakpoint
        let split_info: Vec<&str> = split_line[6].split(',').collect();
        let mapq = split_info[5].parse::<usize>()?;

        let hap = if mapq == 30 {
            0
        } else if &read_id[0..2] == "h1" || &read_id[0..10] == "haplotype1" {
            1
        } else {
            2
        };
        // store haplotype of breakpoint into HashMap<sv_id, HashMap<bp_num. Vec::hap_count>>
        if let Some(bp_map) = hap_db.get_mut(&sv_id) {
            if let Some(hap_cnt_vec) = bp_map.get_mut(&bp_num) {
                hap_cnt_vec[hap] += 1;
            } else {
                let new_hap_cnt_vec = if hap == 0 {
                    vec![1, 0, 0]
                } else if hap == 1 {
                    vec![0, 1, 0]
                } else {
                    vec![0, 0, 1]
                };
                bp_map.insert(bp_num, new_hap_cnt_vec);
            }
        } else {
            let new_hap_cnt_vec = if hap == 0 {
                vec![1, 0, 0]
            } else if hap == 1 {
                vec![0, 1, 0]
                } else {
                vec![0, 0, 1]
            };
            let tmp_bp_map: HashMap<usize, Vec<usize>> = HashMap::from([(bp_num, new_hap_cnt_vec)]);
            hap_db.insert(sv_id, tmp_bp_map);
        }
    }

    // using group_db and hap_db, determine which breakpoint is true
    for (key, value) in group_db.iter() {
        let mut sv_bp1_cnt = vec![0, 0, 0];
        let mut max_sv1_id = vec![String::new(); 3];
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
                        max_sv1_id[0] = sv_id.clone();
                    }
                    if cnt_vec[1] > max_sv1_cnt[1] {
                        max_sv1_cnt[1] = cnt_vec[1];
                        max_sv1_id[1] = sv_id.clone();
                    }
                    if cnt_vec[2] > max_sv1_cnt[2] {
                        max_sv1_cnt[2] = cnt_vec[2];
                        max_sv1_id[2] = sv_id.clone();
                    }
                }
            }
        }

        let mut sv_bp2_cnt = vec![0, 0, 0];
        let mut max_sv2_id = vec![String::new(); 3];
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
                        max_sv2_id[0] = sv_id.clone();
                    }
                    if cnt_vec[1] > max_sv1_cnt[1] {
                        max_sv2_cnt[1] = cnt_vec[1];
                        max_sv2_id[1] = sv_id.clone();
                    }
                    if cnt_vec[2] > max_sv1_cnt[2] {
                        max_sv2_cnt[2] = cnt_vec[2];
                        max_sv2_id[2] = sv_id.clone();
                    }
                }
            }
        }
        // TODO: choose one SV breakpoint in each group 
        println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}", key, value.identical_bp1, value.identical_bp2, sv_bp1_cnt, sv_bp2_cnt);
        let max_index_bp1 = sv_bp1_cnt
        .iter()
        .enumerate()
        .max_by_key(|&(_, item)| item)
        .map(|(index, _)| index);
        
        match max_index_bp1 {
            Some(index) => {
                print!("{}\t", max_sv1_id[index])
            },
            None => eprintln!("No max value"),
        }

        let max_index_bp2 = sv_bp2_cnt
        .iter()
        .enumerate()
        .max_by_key(|&(_, item)| item)
        .map(|(index, _)| index);

        match max_index_bp2 {
            Some(index) => {
                println!("{}", max_sv2_id[index])
            },
            None => eprintln!("No max value"),
        }
    }
    Ok(())
}

