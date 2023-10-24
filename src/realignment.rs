use std::error::Error;
use std::io::BufRead;
use std::collections::HashMap;
use std::collections::HashSet;

use nanomonsv_postprocess::open_file;
#[path = "./smith_waterman.rs"]
mod smith_waterman;

pub fn run(bp_file: &str) -> Result<(), Box<dyn Error>> {
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
    let mut identical_pairs: Vec<(String, String)> = Vec::new();
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
            // eprintln!("{:?}", value1[0].);
            let result1 = match smith_waterman::run(&value1[0], &value2[0]) {
                Ok(id) => id,
                _ => -1.0,
            };
            let result2 = match smith_waterman::run(&value1[0], &value2[1]) {
                Ok(id) => id,
                _ => -1.0,
            };
            let result3 = match smith_waterman::run(&value1[1], &value2[0]) {
                Ok(id) => id,
                _ => -1.0,
            };
            let result4 = match smith_waterman::run(&value1[1], &value2[1]) {
                Ok(id) => id,
                _ => -1.0,
            };

            let bp1_max_id = if result1 > result2 {
                result1
            } else {
                result2
            };
            let bp2_max_id = if result3 > result4 {
                result3
            } else {
                result4
            };

            if bp1_max_id >= threshold && bp2_max_id >= threshold {
                eprintln!("{}\t{}", key1, key2);
                identical_pairs.push((key1.to_string(), key2.to_string()));
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

    let mut group_set: Vec<HashSet<String>> = Vec::new();

    for pairs in groups.iter() {
        let mut id_set: HashSet<String> = HashSet::new();
        for pair in pairs {
            id_set.insert(pair.0.clone());
            id_set.insert(pair.1.clone());
        }
        group_set.push(id_set);
    }

    for group in group_set {
        println!("{:?}", group);
    }
    Ok(())
}
