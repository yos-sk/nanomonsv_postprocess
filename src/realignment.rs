use std::error::Error;
use std::io::BufRead;
use std::collections::HashMap;
use std::collections::HashSet;

use nanomonsv_postprocess::open_file;
#[path = "./smith_waterman.rs"]
mod smith_waterman;

pub fn run(bp_file: &str, support_read_file: &str) -> Result<(Vec<HashSet<String>>), Box<dyn Error>> {
    let bp_reader = open_file(bp_file).expect(&format!("Could not open file {}", bp_file));

    let mut bp_seq: HashMap<String, Vec<Vec<u8>>> = HashMap::new();

    for line in bp_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        // let contig = split_line[0].to_string();
        // let position = split_line[1].parse::<usize>().unwrap();
        let sv_id = split_line[3].to_string();
        let seq = split_line[4];
        
        if let Some(entry) = bp_seq.get_mut(&sv_id) {
            entry.push(seq.as_bytes().to_vec());
        } else {
            bp_seq.insert(sv_id, vec![seq.as_bytes().to_vec()]);
        }
    }
    let mut sorted_pairs: Vec<(&String, &Vec<Vec<u8>>)> = bp_seq.iter().collect();
    sorted_pairs.sort_by_key(|&(k, _)| k);

    let threshold = 99.0;
    let mut identical_pairs: Vec<(String, String, usize)> = Vec::new();
    // realignment of 2 * 2 breakpoint combination
    for (&ref key1, &ref value1) in &sorted_pairs {
        let sv_type1 = &key1[0..1];
        for (&ref key2, &ref value2) in &sorted_pairs {
            if key1 >= key2 {
                continue;
            }
            let sv_type2 = &key2[0..1];
            if sv_type1 != sv_type2 {
                continue;
            }
            // smith_waterman algorithm
            let result1 = match smith_waterman::run(&value1[0], &value2[0]) {
                Ok(id) => id,
                _ => -1.0,
            };
            let result2 = match smith_waterman::run(&value1[0], &value2[1]) {
                Ok(id) => id,
                _ => -1.0,
            };
            let mut bp1_max_id = 0;
            let mut bp2_max_id = 0;
            let mut pattern = 0;
            if result1 > result2 {
                bp1_max_id = result1;
                match smith_waterman::run(&value1[1], &value2[1]) {
                    Ok(id) => bp2_max_id = id;
                    _ => bp2_max_id = -1.0;
                }
            } else {
                bp1_max_id = result2;
                pattern = 1;
                match smith_waterman::run(&value1[1], &value2[0]) {           
                    Ok(id) => bp2_max_id = id;
                    _ => bp_max_id = -1.0;
                }
            }

            if bp1_max_id >= threshold && bp2_max_id >= threshold {
                eprintln!("{}\t{}", key1, key2);
                identical_pairs.push((key1.to_string(), key2.to_string(), pattern));
            }
        }
    }

    // grouping
    let mut groups: Vec<Vec<(String, String)>> = Vec::new();
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

    let mut group_set: Vec<(HashSet<String>, HashSet<String>)> = Vec::new();

    for pairs in groups.iter() {
        let mut id_set_1: HashSet<String> = HashSet::new();
        let mut id_set_2: HashSet<String> = HashSet::new();
        // TODO: asserion of inconsistent case. e.g) (r3, r4, 0), (r3, r5, 0), (r4, r5, 1)
        for pair in pairs {
            match pair.2 {
                0 => {
                    let bp_1_1 = pair.0.clone().push('_').push('1');
                    let bp_2_1 = pair.1.clone().push('_').push('1');
                    let bp_1_2 = pair.0.clone().push('_').push('2');
                    let bp_2_2 = pair.1.clone().push('_').push('2');
                    id_set_1.insert(bp_1_1);
                    id_set_1.insert(bp_2_1);
                    id_set_2.insert(bp_1_2);
                    id_set_2.insert(bp_2_2);
                },
                1 => {
                    let bp_1_1 = pair.0.clone().push('_').push('1');
                    let bp_2_1 = pair.1.clone().push('_').push('2');
                    let bp_1_2 = pair.0.clone().push('_').push('2');
                    let bp_2_2 = pair.1.clone().push('_').push('1');
                    id_set_1.insert(bp_1_1);
                    id_set_1.insert(bp_2_2);
                    id_set_2.insert(bp_1_2);
                    id_set_2.insert(bp_2_1);
                },
                _ => ().
            }
        }
        group_set.push((id_set_1, id_set_2));
    }

    /*
    for group in group_set {
        println!("{:?}", group);
    }
    */
    classify_haplotype(&group_set, support_read_file);
    Ok(())
}

fn classify_haplotype(support_read_file: &str, nanomonsv_bp_file: &str) -> Result<(), Box<dyn Error>> {
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
    let hap_db: HashMap<String, Vec<usize>> = HashMap::new();

    for line in bp_reader.lines() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let read_id = split_line[3];
        let sv_id = match read_db.get(&read_id) {
            Some(value) => value,
            None => {
                continue;
            },
        };
        // TODO: break pointごとに，haplotyp1/2/unassignedの議論をしなければならなくて，それをやるには，breakpointの対応関係を明らかにしなければならない。
        let split_info: Vec<&str> = split_line[6].split(',').collect();
        let mapq = split_info[5].parse::<usize>()?;

        let hap = if read_id[0..2] == "h1" || read_id[0..10] == "haplotype1" {
            1
        } else {
            2
        };

        let f_hap = if mapq == 30 {
            0
        } else {
            hap
        };

        match hap_db.get_mut(&sv_id) {
            Some(value) => {
                *value[f_hap] += 1;
            },
            None => {
                if f_hap == 0 {
                    hap_db.insert(sv_id, vec![1, 0, 0]);
                } else if f_hap == 1 {
                    hap_db.insert(sv_id, vec![0, 1, 0]);
                } else {
                    hap_db.insert(sv_id, vec![0, 0, 1]);
                }
            }
        }
    }

}

