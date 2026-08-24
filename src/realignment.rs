use rust_htslib::bam::{Header, IndexedReader, Read};
use std::collections::BTreeMap;
use std::collections::BTreeSet;
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

/// Which breakpoint of which SV is the same physical junction as this pair's.
/// `BTreeSet` rather than `HashSet`: these sets are iterated when the record is
/// printed, and a randomised hash order would make the output irreproducible.
#[derive(Eq, PartialEq, Clone)]
struct IdenticalInfo {
    identical_bp1: BTreeSet<(String, usize)>,
    identical_bp2: BTreeSet<(String, usize)>,
}
impl IdenticalInfo {
    fn new() -> Self {
        IdenticalInfo {
            identical_bp1: BTreeSet::new(),
            identical_bp2: BTreeSet::new(),
        }
    }

    fn add_identical_bp1(&mut self, sv_id: String, bp_num: usize) {
        self.identical_bp1.insert((sv_id, bp_num));
    }

    fn add_identical_bp2(&mut self, sv_id: String, bp_num: usize) {
        self.identical_bp2.insert((sv_id, bp_num));
    }
}

/// Which haplotype a breakpoint sits on, for the two breakpoints of one call,
/// read in the orientation the group assigned to that call.
type Combo = (u8, u8);

/// contig -> 1 or 2, from the bam_refiner haplotype contig lists.
fn load_contig_haplotype(
    hap1_list: &str,
    hap2_list: &str,
) -> Result<HashMap<String, u8>, Box<dyn Error>> {
    let mut db: HashMap<String, u8> = HashMap::new();
    for (path, hap) in [(hap1_list, 1u8), (hap2_list, 2u8)] {
        let reader = open_file(path).expect(&format!("Could not open {}", path));
        for line in reader.lines() {
            let line = line?;
            if let Some(contig) = line.split_whitespace().next() {
                db.insert(contig.to_string(), hap);
            }
        }
    }
    Ok(db)
}

pub fn run(
    input_bed: &str,
    support_read_file: &str,
    bam_file: &str,
    min_identity: f64,
    min_length: usize,
    hap1_list: &str,
    hap2_list: &str,
    max_position_diff: usize,
) -> Result<(), Box<dyn Error>> {
    let contig_hap: HashMap<String, u8> = if hap1_list.is_empty() || hap2_list.is_empty() {
        eprintln!(
            "Warning: no haplotype contig lists given; calls will be grouped without \
             the haplotype-combination constraint"
        );
        HashMap::new()
    } else {
        load_contig_haplotype(hap1_list, hap2_list)?
    };
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

    // Compare every SV against every other one and keep, for each candidate
    // pair, how well its two breakpoints realign. Iteration goes over a sorted
    // vector rather than the HashSet directly: the set of pairs would be the
    // same either way, but the order decides ties later on, and a randomised
    // hash order would make the run irreproducible.
    let mut sv_infos: Vec<&SVInfo> = bp_info_db.iter().collect();
    sv_infos.sort_by(|a, b| a.sv_id.cmp(&b.sv_id));
    let sv_by_id: HashMap<&str, &SVInfo> =
        sv_infos.iter().map(|s| (s.sv_id.as_str(), *s)).collect();

    // (sv_id_1, sv_id_2, pattern, identity)
    let mut identical_pairs: Vec<(String, String, usize, f64)> = Vec::new();
    for sv_info_1 in &sv_infos {
        let sv_type1 = &sv_info_1.sv_id[0..1];
        for sv_info_2 in &sv_infos {
            if sv_info_1.sv_id >= sv_info_2.sv_id {
                continue;
            }
            let sv_type2 = &sv_info_2.sv_id[0..1];

            if sv_type1 == "d" && sv_type2 == "i" {
                continue;
            }
            // smith_waterman algorithm
            let result1 = smith_waterman::run(
                &sv_info_1.bp1_seq,
                &sv_info_2.bp1_seq,
                min_identity,
                min_length,
            )
            .unwrap_or(None);
            let result2 = smith_waterman::run(
                &sv_info_1.bp1_seq,
                &sv_info_2.bp2_seq,
                min_identity,
                min_length,
            )
            .unwrap_or(None);

            // A pair is scored by its weaker breakpoint: both junctions have to
            // look the same for the two calls to be the same event, so the
            // lower of the two identities is what should be ranked.
            let (pattern, first, second) = if result1.is_some() && result2.is_none() {
                (
                    0,
                    result1,
                    smith_waterman::run(
                        &sv_info_1.bp2_seq,
                        &sv_info_2.bp2_seq,
                        min_identity,
                        min_length,
                    )
                    .unwrap_or(None),
                )
            } else {
                (
                    1,
                    result2,
                    smith_waterman::run(
                        &sv_info_1.bp2_seq,
                        &sv_info_2.bp1_seq,
                        min_identity,
                        min_length,
                    )
                    .unwrap_or(None),
                )
            };

            let Some(second_identity) = second else {
                continue;
            };

            // Two calls are the same event seen once per haplotype, so they must
            // sit on different haplotype combinations. Same-combination calls --
            // (1,1) with (1,1), say -- are separate loci such as repeat copies,
            // and collapsing them is what over-filtered the centromeres.
            if !contig_hap.is_empty() {
                let combo_1 = oriented_combo(sv_info_1, &contig_hap, false);
                let combo_2 = oriented_combo(sv_info_2, &contig_hap, pattern == 1);
                match (combo_1, combo_2) {
                    (Some(c1), Some(c2)) if c1 != c2 => (),
                    // Unknown contig: too little information to declare the
                    // calls redundant, so leave them alone.
                    _ => continue,
                }
                // Sharing a haplotype on one side means sharing that breakpoint.
                if !shared_haplotype_sides_agree(
                    sv_info_1,
                    sv_info_2,
                    pattern == 1,
                    &contig_hap,
                    max_position_diff,
                ) {
                    continue;
                }
            }

            let identity = match first {
                Some(first_identity) => first_identity.min(second_identity),
                None => second_identity,
            };
            eprintln!(
                "{}\t{}\t{}\t{:.3}",
                sv_info_1.sv_id, sv_info_2.sv_id, pattern, identity
            );
            identical_pairs.push((
                sv_info_1.sv_id.clone(),
                sv_info_2.sv_id.clone(),
                pattern,
                identity,
            ));
        }
    }

    let groups = build_groups(&identical_pairs, &sv_by_id, &contig_hap);

    let mut group_db: BTreeMap<Vec<String>, IdenticalInfo> = BTreeMap::new();
    for members in groups {
        let mut info = IdenticalInfo::new();
        let mut key: Vec<String> = Vec::new();
        for (sv_id, swapped) in &members {
            // A swapped call contributes its bp2 to side 1 and its bp1 to side 2.
            let (side1_bp, side2_bp) = if *swapped { (2, 1) } else { (1, 2) };
            info.add_identical_bp1(sv_id.clone(), side1_bp);
            info.add_identical_bp2(sv_id.clone(), side2_bp);
            key.push(sv_id.clone());
        }
        group_db.insert(key, info);
    }

    let _ = classify_haplotype_fetch_bam(&group_db, &bp_info_db, support_read_file, bam_file);
    Ok(())
}

fn classify_haplotype_fetch_bam(
    group_db: &BTreeMap<Vec<String>, IdenticalInfo>,
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
        let mut sv_bp1_map: HashMap<usize, Vec<(String, usize, usize)>> = HashMap::new();
        // let mut max_sv1_id = vec![(String::new(), 0); 3];
        // let mut max_sv1_cnt = vec![0, 0, 0];
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
                    if cnt_vec[1] > 0 {
                        sv_bp1_map.entry(1)
                                    .or_insert_with(Vec::new)
                                    .push((sv_id, bp_num, cnt_vec[1]));
                    } else if cnt_vec[2] > 0 {
                        sv_bp1_map.entry(2)
                                    .or_insert_with(Vec::new)
                                    .push((sv_id, bp_num, cnt_vec[2]));
                    }
                }
            }
        }

        let mut sv_bp2_cnt = vec![0, 0, 0];
        // let mut max_sv2_id = vec![(String::new(), 0); 3];
        // let mut max_sv2_cnt = vec![0, 0, 0];
        let mut sv_bp2_map: HashMap<usize, Vec<(String, usize, usize)>> = HashMap::new();
        // breakpoint 2
        for item in value.identical_bp2.iter() {
            let sv_id = (&item.0).to_string();
            let bp_num = item.1;
            if let Some(cnt_map) = hap_db.get(&sv_id) {
                if let Some(cnt_vec) = cnt_map.get(&bp_num) {
                    sv_bp2_cnt[0] += cnt_vec[0];
                    sv_bp2_cnt[1] += cnt_vec[1];
                    sv_bp2_cnt[2] += cnt_vec[2];
                    if cnt_vec[1] > 0 {
                        sv_bp2_map.entry(1)
                                    .or_insert_with(Vec::new)
                                    .push((sv_id, bp_num, cnt_vec[1]));
                    } else if cnt_vec[2] > 0 {
                        sv_bp2_map.entry(2)
                                    .or_insert_with(Vec::new)
                                    .push((sv_id, bp_num, cnt_vec[2]));
                    }
                    /*
                    // Unassigned
                    if cnt_vec[0] > max_sv2_cnt[0] {
                        max_sv2_cnt[0] = cnt_vec[0];
                        max_sv2_id[0] = (sv_id.clone(), bp_num);
                    } else if cnt_cev[0] == max_sv2_cnt[0] {
                        if sv_id == max_sv1_id[0] {
                            max_sv2_id[0] = (sv_id.clone(), bp_num);
                        }
                    }
                    // haplotype 1
                    if cnt_vec[1] > max_sv2_cnt[1] {
                        max_sv2_cnt[1] = cnt_vec[1];
                        max_sv2_id[1] = (sv_id.clone(), bp_num);
                    }
                    // haplotype 2
                    if cnt_vec[2] > max_sv2_cnt[2] {
                        max_sv2_cnt[2] = cnt_vec[2];
                        max_sv2_id[2] = (sv_id.clone(), bp_num);
                    }
                    */
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
            print!("\tAmbiguous");
            if let Some(vec) = sv_bp1_map.get(&1) {
                for (i, (sv_id, bp, n_read)) in vec.into_iter().enumerate() {
                    if i == 0 {
                        print!("\t{},{},{}", sv_id, bp, n_read);
                    } else {
                        print!(";{},{},{}", sv_id, bp, n_read);
                    }
                }
            }
            if let Some(vec) = sv_bp1_map.get(&2) {
                for (sv_id, bp, n_read) in vec {
                    print!(";{},{},{}", sv_id, bp, n_read);
                }
            }
        } else {
            assigned += 1;
            if sv_bp1_cnt[1] > 0 {
                print!("\thaplotype1");
                if let Some(vec) = sv_bp1_map.get(&1) {
                    for (i, (sv_id, bp, n_read)) in vec.into_iter().enumerate() {
                        if i == 0 {
                            print!("\t{},{},{}", sv_id, bp, n_read);
                        } else {
                            print!(";{},{},{}", sv_id, bp, n_read);
                        }
                    }
                }
            } else {
                print!("\thaplotype2");
                if let Some(vec) = sv_bp1_map.get(&2) {
                    for (i, (sv_id, bp, n_read)) in vec.into_iter().enumerate() {
                        if i == 0 {
                            print!("\t{},{},{}", sv_id, bp, n_read);
                        } else {
                            print!(";{},{},{}", sv_id, bp, n_read);
                        }
                    }
                }
            }
        }

        if sv_bp2_cnt[1] == 0 && sv_bp2_cnt[2] == 0 {
            unassigned += 1;
            println!("\tUnassigned\t-");
        } else if sv_bp2_cnt[1] > 0 && sv_bp2_cnt[2] > 0 {
            ambiguous += 1;
            print!("\tAmbiguous");
            if let Some(vec) = sv_bp2_map.get(&1) {
                for (i, (sv_id, bp, n_read)) in vec.into_iter().enumerate() {
                    if i == 0 {
                        print!("\t{},{},{}", sv_id, bp, n_read);
                    } else {
                        print!(";{},{},{}", sv_id, bp, n_read);
                    }
                }
            }
            if let Some(vec) = sv_bp2_map.get(&2) {
                for (sv_id, bp, n_read) in vec {
                    print!(";{},{},{}", sv_id, bp, n_read);
                }
            }
            print!("\n");
        } else {
            assigned += 1;
            if sv_bp2_cnt[1] > 0 {
                print!("\thaplotype1");
                if let Some(vec) = sv_bp2_map.get(&1) {
                    for (i, (sv_id, bp, n_read)) in vec.into_iter().enumerate() {
                        if i == 0 {
                            print!("\t{},{},{}", sv_id, bp, n_read);
                        } else {
                            print!(";{},{},{}", sv_id, bp, n_read);
                        }
                    }
                }
            } else {
                print!("\thaplotype2");
                if let Some(vec) = sv_bp2_map.get(&2) {
                    for (i, (sv_id, bp, n_read)) in vec.into_iter().enumerate() {
                        if i == 0 {
                            print!("\t{},{},{}", sv_id, bp, n_read);
                        } else {
                            print!(";{},{},{}", sv_id, bp, n_read);
                        }
                    }
                }
            }
            print!("\n");
        }
    }

    eprintln!("Assigned: {}", assigned);
    eprintln!("Unassigned: {}", unassigned);
    eprintln!("Ambiguous: {}", ambiguous);
    Ok(())
}

/// The two breakpoints of a call, read in the given orientation.
fn oriented_breakpoints(sv: &SVInfo, swapped: bool) -> [(&str, usize); 2] {
    if swapped {
        [
            (sv.bp2_contig.as_str(), sv.bp2_pos),
            (sv.bp1_contig.as_str(), sv.bp1_pos),
        ]
    } else {
        [
            (sv.bp1_contig.as_str(), sv.bp1_pos),
            (sv.bp2_contig.as_str(), sv.bp2_pos),
        ]
    }
}

/// Whether the sides on which two calls share a haplotype point at the same locus.
///
/// If the calls really are one event seen once per haplotype, then wherever they
/// land on the *same* haplotype they must be the same breakpoint: same contig,
/// same position. On BL2009 HiFi 42 of 46 such sides sat within 100 bp of each
/// other (median 0 bp); the handful that did not were on a different contig
/// altogether or megabases away, i.e. separate loci that merely share sequence.
/// The opposite-haplotype side is not checked — different contig and a different
/// coordinate system are expected there.
fn shared_haplotype_sides_agree(
    sv1: &SVInfo,
    sv2: &SVInfo,
    swapped2: bool,
    hap: &HashMap<String, u8>,
    max_position_diff: usize,
) -> bool {
    let side_1 = oriented_breakpoints(sv1, false);
    let side_2 = oriented_breakpoints(sv2, swapped2);

    for side in 0..2 {
        let (contig_1, pos_1) = side_1[side];
        let (contig_2, pos_2) = side_2[side];
        let (hap_1, hap_2) = (hap.get(contig_1), hap.get(contig_2));
        // An unplaced contig leaves too little information to call the two
        // redundant, so refuse rather than guess.
        let (Some(hap_1), Some(hap_2)) = (hap_1, hap_2) else {
            return false;
        };
        if hap_1 != hap_2 {
            continue;
        }
        if contig_1 != contig_2 || pos_1.abs_diff(pos_2) > max_position_diff {
            return false;
        }
    }
    true
}

/// The haplotype combination of a call, read in the given orientation.
fn oriented_combo(sv: &SVInfo, hap: &HashMap<String, u8>, swapped: bool) -> Option<Combo> {
    let h1 = *hap.get(&sv.bp1_contig)?;
    let h2 = *hap.get(&sv.bp2_contig)?;
    Some(if swapped { (h2, h1) } else { (h1, h2) })
}

/// Group the calls that realignment says are the same event, highest identity
/// first, ties broken by SV id.
///
/// One somatic junction is called once per haplotype *per breakpoint*, so up to
/// four calls describe it: (1,1), (1,2), (2,1) and (2,2). A group therefore has
/// four slots and each holds at most one call — the first claimant, which by the
/// processing order is the one with the highest realignment identity.
///
/// The previous implementation instead merged every pair that shared a call into
/// one transitive group. In a centromere or another high-identity repeat the
/// breakpoint sequences all look alike, so a single group swallowed whole
/// families of calls and `filt` then kept only one of them. Two guards stop that
/// here: a call can only join a group if it realigns against *every* member
/// (no chaining through an intermediate), and only if its haplotype slot is
/// still free (repeat copies share a slot, so they cannot collapse).
///
/// `swapped` records the orientation a call was given inside its group: a
/// swapped call's bp2 is the junction side that the others call bp1.
fn build_groups(
    candidates: &[(String, String, usize, f64)],
    sv_by_id: &HashMap<&str, &SVInfo>,
    contig_hap: &HashMap<String, u8>,
) -> Vec<Vec<(String, bool)>> {
    let mut edges: Vec<&(String, String, usize, f64)> = candidates.iter().collect();
    edges.sort_by(|a, b| {
        b.3.partial_cmp(&a.3)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.0.cmp(&b.0))
            .then(a.1.cmp(&b.1))
    });

    let pattern_of: HashMap<(&str, &str), usize> = candidates
        .iter()
        .map(|(a, b, p, _)| ((a.as_str(), b.as_str()), *p))
        .collect();
    let pattern_between = |a: &str, b: &str| -> Option<usize> {
        if a <= b {
            pattern_of.get(&(a, b)).copied()
        } else {
            pattern_of.get(&(b, a)).copied()
        }
    };

    // group index -> members; a call belongs to exactly one group.
    let mut groups: Vec<Vec<(String, bool)>> = Vec::new();
    let mut group_of: HashMap<String, usize> = HashMap::new();

    for (id_a, id_b, pattern, _identity) in edges {
        let ga = *group_of
            .entry(id_a.clone())
            .or_insert_with(|| {
                groups.push(vec![(id_a.clone(), false)]);
                groups.len() - 1
            });
        let gb = *group_of
            .entry(id_b.clone())
            .or_insert_with(|| {
                groups.push(vec![(id_b.clone(), false)]);
                groups.len() - 1
            });
        if ga == gb {
            continue;
        }

        // Orient the second group so that this pair's correspondence holds:
        // differing orientations mean bp1 of one matches bp2 of the other.
        let orient_a = groups[ga]
            .iter()
            .find(|(id, _)| id == id_a)
            .map(|(_, s)| *s)
            .unwrap_or(false);
        let orient_b = groups[gb]
            .iter()
            .find(|(id, _)| id == id_b)
            .map(|(_, s)| *s)
            .unwrap_or(false);
        let wanted_b = orient_a ^ (*pattern == 1);
        let flip = wanted_b ^ orient_b;

        let merged_b: Vec<(String, bool)> = groups[gb]
            .iter()
            .map(|(id, s)| (id.clone(), s ^ flip))
            .collect();

        // Every cross pair must exist and agree with the orientations, so that
        // a group is a clique rather than a chain.
        let consistent = groups[ga].iter().all(|(m, om)| {
            merged_b.iter().all(|(n, on)| {
                pattern_between(m, n) == Some(if om == on { 0 } else { 1 })
            })
        });
        if !consistent {
            continue;
        }

        // One call per haplotype slot.
        if !contig_hap.is_empty() {
            let mut slots: BTreeSet<Combo> = BTreeSet::new();
            let mut all_known = true;
            for (id, swapped) in groups[ga].iter().chain(merged_b.iter()) {
                match sv_by_id
                    .get(id.as_str())
                    .and_then(|sv| oriented_combo(sv, contig_hap, *swapped))
                {
                    Some(combo) => {
                        if !slots.insert(combo) {
                            all_known = false;
                            break;
                        }
                    }
                    None => {
                        all_known = false;
                        break;
                    }
                }
            }
            if !all_known {
                continue;
            }
        } else if groups[ga].len() + merged_b.len() > 4 {
            continue;
        }

        for (id, _) in merged_b.iter() {
            group_of.insert(id.clone(), ga);
        }
        let mut merged = std::mem::take(&mut groups[ga]);
        merged.extend(merged_b);
        merged.sort();
        groups[ga] = merged;
        groups[gb].clear();
    }

    let mut out: Vec<Vec<(String, bool)>> =
        groups.into_iter().filter(|g| g.len() > 1).collect();
    out.sort();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sv(id: &str, c1: &str, c2: &str) -> SVInfo {
        sv_at(id, c1, 100, c2, 200)
    }

    fn sv_at(id: &str, c1: &str, p1: usize, c2: &str, p2: usize) -> SVInfo {
        let mut info = SVInfo::new();
        info.add_bp1_info(id.to_string(), c1.to_string(), p1, b"A");
        info.add_bp2_info(id.to_string(), c2.to_string(), p2, b"A");
        info
    }

    fn locus_hap() -> HashMap<String, u8> {
        let mut hap = HashMap::new();
        hap.insert("h1_a".to_string(), 1u8);
        hap.insert("h1_b".to_string(), 1u8);
        hap.insert("h2_a".to_string(), 2u8);
        hap
    }

    #[test]
    fn a_shared_haplotype_side_must_be_the_same_locus() {
        let hap = locus_hap();
        // Same hap1 breakpoint, hap2 partner differs: the real redundancy.
        let a = sv_at("a", "h1_a", 1_000, "h1_a", 5_000);
        let b = sv_at("b", "h1_a", 1_050, "h2_a", 5_000);
        assert!(shared_haplotype_sides_agree(&a, &b, false, &hap, 100));

        // Same haplotype but megabases apart: separate loci.
        let far = sv_at("c", "h1_a", 3_000_000, "h2_a", 5_000);
        assert!(!shared_haplotype_sides_agree(&a, &far, false, &hap, 100));

        // Same haplotype, different contig: separate loci.
        let other = sv_at("d", "h1_b", 1_000, "h2_a", 5_000);
        assert!(!shared_haplotype_sides_agree(&a, &other, false, &hap, 100));
    }

    #[test]
    fn the_opposite_haplotype_side_is_not_position_checked() {
        let hap = locus_hap();
        // hap1 side agrees; the hap2 partner sits at a quite different
        // coordinate, which is normal across haplotypes.
        let a = sv_at("a", "h1_a", 1_000, "h1_a", 5_000);
        let b = sv_at("b", "h1_a", 1_000, "h2_a", 4_000_000);
        assert!(shared_haplotype_sides_agree(&a, &b, false, &hap, 100));
    }

    struct Fixture {
        svs: Vec<SVInfo>,
        hap: HashMap<String, u8>,
    }

    impl Fixture {
        /// Calls named by their haplotype combination, e.g. "a11" sits on
        /// (hap1, hap1).
        fn new(names: &[(&str, u8, u8)]) -> Self {
            let mut hap = HashMap::new();
            hap.insert("h1_left".to_string(), 1u8);
            hap.insert("h2_left".to_string(), 2u8);
            hap.insert("h1_right".to_string(), 1u8);
            hap.insert("h2_right".to_string(), 2u8);
            let svs = names
                .iter()
                .map(|(id, a, b)| {
                    let c1 = if *a == 1 { "h1_left" } else { "h2_left" };
                    let c2 = if *b == 1 { "h1_right" } else { "h2_right" };
                    sv(id, c1, c2)
                })
                .collect();
            Fixture { svs, hap }
        }

        fn group(&self, candidates: &[(String, String, usize, f64)]) -> Vec<Vec<String>> {
            let by_id: HashMap<&str, &SVInfo> =
                self.svs.iter().map(|s| (s.sv_id.as_str(), s)).collect();
            build_groups(candidates, &by_id, &self.hap)
                .into_iter()
                .map(|g| g.into_iter().map(|(id, _)| id).collect())
                .collect()
        }
    }

    fn cand(a: &str, b: &str, identity: f64) -> (String, String, usize, f64) {
        (a.to_string(), b.to_string(), 0, identity)
    }

    #[test]
    fn all_four_haplotype_combinations_form_one_group() {
        let f = Fixture::new(&[("a", 1, 1), ("b", 1, 2), ("c", 2, 1), ("d", 2, 2)]);
        let ids = ["a", "b", "c", "d"];
        let mut candidates = Vec::new();
        for (i, x) in ids.iter().enumerate() {
            for y in ids.iter().skip(i + 1) {
                candidates.push(cand(x, y, 99.0));
            }
        }
        assert_eq!(
            f.group(&candidates),
            vec![vec![
                "a".to_string(),
                "b".to_string(),
                "c".to_string(),
                "d".to_string()
            ]]
        );
    }

    #[test]
    fn a_slot_taken_by_a_better_pair_is_not_taken_again() {
        // b and c both want the (1,2) slot next to a; b aligns better.
        let f = Fixture::new(&[("a", 1, 1), ("b", 1, 2), ("c", 1, 2)]);
        let candidates = vec![
            cand("a", "b", 99.5),
            cand("a", "c", 99.0),
            cand("b", "c", 99.0),
        ];
        assert_eq!(
            f.group(&candidates),
            vec![vec!["a".to_string(), "b".to_string()]]
        );
    }

    #[test]
    fn same_combination_calls_never_group() {
        // A repeat family: every call is (1,1), so none of them can be the same
        // event called on the other haplotype.
        let f = Fixture::new(&[("a", 1, 1), ("b", 1, 1), ("c", 1, 1)]);
        let candidates = vec![
            cand("a", "b", 99.0),
            cand("a", "c", 99.0),
            cand("b", "c", 99.0),
        ];
        assert!(f.group(&candidates).is_empty());
    }

    #[test]
    fn a_chain_without_a_closing_edge_does_not_merge() {
        // a~b and b~d realign, a~d does not: d must not join through b.
        let f = Fixture::new(&[("a", 1, 1), ("b", 1, 2), ("d", 2, 2)]);
        let candidates = vec![cand("a", "b", 99.5), cand("b", "d", 99.0)];
        assert_eq!(
            f.group(&candidates),
            vec![vec!["a".to_string(), "b".to_string()]]
        );
    }

    #[test]
    fn grouping_is_independent_of_input_order() {
        let f = Fixture::new(&[("a", 1, 1), ("b", 1, 2), ("c", 2, 1), ("d", 2, 2)]);
        let ids = ["a", "b", "c", "d"];
        let mut forward = Vec::new();
        for (i, x) in ids.iter().enumerate() {
            for y in ids.iter().skip(i + 1) {
                forward.push(cand(x, y, 99.0));
            }
        }
        let mut reverse = forward.clone();
        reverse.reverse();
        assert_eq!(f.group(&forward), f.group(&reverse));
    }

    /// A deletion has both breakpoints on one contig, so only (1,1) and (2,2)
    /// exist and a group can never exceed two calls.
    #[test]
    fn an_intra_contig_event_pairs_at_most_two_calls() {
        let mut hap = HashMap::new();
        hap.insert("h1".to_string(), 1u8);
        hap.insert("h2".to_string(), 2u8);
        let svs = vec![sv("d_1", "h1", "h1"), sv("d_2", "h2", "h2"), sv("d_3", "h1", "h1")];
        let by_id: HashMap<&str, &SVInfo> = svs.iter().map(|s| (s.sv_id.as_str(), s)).collect();
        let candidates = vec![
            cand("d_1", "d_2", 99.5),
            cand("d_2", "d_3", 99.0),
            cand("d_1", "d_3", 99.0),
        ];
        let groups = build_groups(&candidates, &by_id, &hap);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 2);
    }
}
