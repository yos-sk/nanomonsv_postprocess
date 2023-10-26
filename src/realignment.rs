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
#[derive(Eq, Hash, PartialEq, Clone)]
struct IdenticalInfo {
    identical_bp1: Vec<(String, usize)>,
    identical_bp2: Vec<(String, usize)>,
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


pub fn run(bp_file: &str) -> Result<(), Box<dyn Error>> {
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

    let mut group_db: HashMap<Vec<String>, IdenticalInfo> = Vec::new();
    let mut group_db_key: Vec<String> = Vec::new();

    for pairs in groups.iter() {
        let mut info = IdenticalInfo::new();
        // TODO: asserion of inconsistent case. e.g) (r3, r4, 0), (r3, r5, 0), (r4, r5, 1)
        for pair in pairs {
            match pair.2 {
                0 => {
                    IdenticalInfo.add_identical_bp1(pair.0.clone(), 1);
                    IdenticalInfo.add_identical_bp2(pair.0.clone(), 2);
                    IdenticalInfo.add_identical_bp1(pair.1.clone(), 1);
                    IdenticalInfo.add_identical_bp2(pair.1.clone(), 2);
                },
                1 => {
                    IdenticalInfo.add_identical_bp1(pair.0.clone(), 1);
                    IdenticalInfo.add_identical_bp2(pair.0.clone(), 2);
                    IdenticalInfo.add_identical_bp1(pair.1.clone(), 2);
                    IdenticalInfo.add_identical_bp2(pair.1.clone(), 1);
                },
                _ => ().
            }
            if !group_db_key.contains(&pair.0) {
                group_db_key.push(pair.0.clone());
            }
            if !gorup_db_key.contains(&pair.1) {
                group_db_key.push(pair.1.clone());
            }
        }
        group_db.insert(group_db_key, IdenticalInfo.clone());
    }

    /*
    for group in group_set {
        println!("{:?}", group);
    }
    */
    //classify_haplotype(&group_set, support_read_file);
    Ok(())
}

/*
fn classify_haplotype(group_set: &Vec<(HashSet<String>, HashSet<String>)>, support_read_file: &str, nanomonsv_bp_file: &str) -> Result<(), Box<dyn Error>> {
    let reader = open_file(support_read_file).expect(&format!("Could not open {}", support_read_file));
    let read_db: HashMap<String, String> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let sv_id = split_line[7].to_string();
        let read_id = split_line[8].to_string();
        read_db.insert(read_id, sv_id);
    }

    let bp_reader = open_file(nanomonsv_bp_file).expect(&format!("Could not open {}", nanomonsv_bp_file));
    let hap_db: HashMap<String, HashMap<Vec<usize>>> = HashMap::new();

    let read_set: HashSet<String> = HashSet::new();

    for line in bp_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let read_id = split_line[3].to_string();
        let sv_id = match read_db.get(&read_id) {
            Some(value) => value,
            None => {
                continue;
            },
        };
        // breakpoint numberを含んだidでgroup_setを検索しなければならない。
        let bp_num = if read_set.contains_key(read_id) {
            sv_id.push('_').push('2')
        } else {
            sv_id.push('_').push('1)
        };

        let sv_hap_num = match group_set.get(&bp_num) {
            Some(value) => {
                if value.0.contains(&bp_num) {
                    1
                } else if value.1.contains(&bp_num) {
                    2
                } else {
                    eprintln!("Failed to group breakpoints");
                    continue;
                }
            },
            None => {
                eprintln!("Failed to group SVs");
                continue;
            },
        };
                
        // TODO: break pointごとに，haplotyp1/2/unassignedの議論をしなければならなくて，それをやるには，breakpointの対応関係を明らかにしなければならない。
        let split_info: Vec<&str> = split_line[6].split(',').collect();
        let mapq = split_info[5].parse::<usize>()?;

        let hap = if mapq == 30 {
            0
        } else if read_id[0..2] == "h1" || read_id[0..10] == "haplotype1" {
            1
        } else {
            2
        };

        if let Some(bp_map) = hap_db.get_mut(&sv_id) {
            if let Some(hap_cnt_vec) = bp_map.get_mut(&sv_hap_num) {
                *hap_cnt_vec[hap] += 1;
            } else {
                let new_hap_cnt_vec = if hap == 0 {
                    vec![1, 0, 0]
                } else if hap == 1 {
                    vec![0, 1, 0]
                } else {
                    vec![0, 0, 1]
                };
                bp_map.insert(sv_hap_num, new_hap_cnt_vec);
            }
        } else {
            let new_hap_cnt_vec = if hap == 0 {
                vec![1, 0, 0]
            } else if hap == 1 {
                vec![0, 1, 0]
            } else {
                vec![0, 0, 1]
            };
            let tmp_bp_map: HashMap<usize, Vec<usize> = hashmap! {
                sv_hap_num => new_hap_cnt_vec,
            };
            hap_db.insert(sv_id, tmp_bp_map);
        }
    }

    

}
*/
