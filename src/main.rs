//
extern crate bio;
extern crate itertools;
use std::collections::HashMap;
use std::cmp;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use bio::io::fastq;
use bio::alignment::pairwise::*;

// fn check(rec: &fastq::Record, read: &str) -> (u16, Vec<(usize, char, char)>) {
//     let mut distance : u16 = 0;
//     let qual = rec.qual();
//     let mut dif : Vec<(usize, char, char)> = vec![];
//     let mut index : usize = 0;

//     for (i, j) in String::from_utf8_lossy(rec.seq()).chars().dropping(8).zip(read.chars()) {
//         if qual[index] > 63 {
//             if i != j {
//                 dif.push((index, i, j));
//                 distance += 1;
//             }
//         }
//         else {
//             distance += 1;
//         }
//         index += 1;
//     }
//     (distance, dif)
// }

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .map(|t| match t {
                 'A' => 'T',
                 'T' => 'A',
                 'G' => 'C',
                 'C' => 'G',
                 _ => 'N',
             }).rev().collect::<String>()
}

fn qual_check(a: &[u8], b: &[u8]) -> bool {
    for (i, j) in a.iter().zip(b.iter()) {
        if i < j {
            continue;
        }
        return false;
    }
    return true
}

fn data_stat(results: &HashMap<String, (String, Vec<u8>)>, output_file: &str) -> Result<String, Box<Error>> {

    // statistics on the datasets
    let wt_pac = "AGAGAAGATTTATCTGAAGTCGTTACGCGAG";
    let mut diff_counts : [usize; 31] = [0; 31];
    let mut diff_freq : [usize; 31] = [0; 31];
    let mut output = try!(File::create(output_file));
    let mut pac_stat = HashMap::new();
    for (_, pac_info) in results {

        let ref pac = pac_info.0;
        let ref qual = pac_info.1;

        // mutation statistics
        let mut index = 0;
        let mut distance = 0;

        for (i, j) in pac.chars().zip(wt_pac.chars()) {
            if qual[index] > 63 && i != j {
                diff_freq[index] += 1;
                distance += 1;
            }
            index += 1;
        }
        diff_counts[distance] += 1;
        if distance > 8 {
            println!("# {} {}", distance,  pac);
        }

        // pac sites statistics
        if pac_stat.contains_key(pac) {
            *pac_stat.get_mut(pac).unwrap() += 1;
        }
        else {
            pac_stat.insert(pac, 1);
        }

    }
    println!("# Overall statistics:");
    for i in 0..31 {
        println!("# {}\t{}", i, diff_counts[i]);
    }

    println!("# Per-base statistics:");
    for i in 0..31 {
        println!("# {}\t{}", i, diff_freq[i]);
    }

    //try!(write!(output, "{}", "# pac counts:\n"));

    for (pac, counts) in &pac_stat {
        let mut hamming = 0;
        for (i, j) in pac.chars().zip(wt_pac.chars()) {
            if i != j {
                hamming += 1;
            }
        }
        try!(write!(output, "{} {} {}\n", pac, hamming, counts));
    }

    Ok("Done".into())
}

fn main() {
    let args : Vec<String> = env::args().collect();

    let file1 = fastq::Reader::from_file(&args[1]).unwrap();
    let file2 = fastq::Reader::from_file(&args[2]).unwrap();

    let mut num_records = 0;
    let mut num_duplicates = 0;
    let mut num_qual_skip = 0;
    let mut results : HashMap<String, (String, Vec<u8>)>= HashMap::new();


    let wt_read1 = if &args[3] == "M" {b"ACTAAGTGAGATGAATATGGCGGCACCAAAGGGCAACCGATTTTGGGAGGCCCGCAGTAGTCATGGGCGAAATCCTAAATTCGAATCGCCTGAGGCGCTGTGGGCTGCTTGTTGTGAA"}
    else {b"AAGTGAGATGAATATGGCGGCACCAAAGGGCAACCGATTTTGGGAGGCCCGCAGTAGTCATGGGCGAAATCCTAAATTCGAATCGCCTGAGGCGCTGTGGGCTGCTTGTTGTGAATAC"};

    for (record1, record2) in file1.records().zip(file2.records()) {

        // take read1, filter low quality reads
        let read1 = record1.unwrap();
        let desc = read1.id().unwrap().split(":").skip(5).collect::<Vec<&str>>();
        let description = desc[0].to_string() + ":" + desc[1];
        let mut trim = 124;
        let mut am = "".to_string();

        for i in 0..120 {
            if qual_check(&read1.qual()[i .. i+5], &[63, 63, 63, 63, 63]) {
                trim = i+1;
                println!("# {} {}: Read 1 trimmed at {}.", num_records, description, trim);
                break;
            }
        }

        if trim < 18 {
            println!("# {}: Useful read too short. Skipping. L = {}", num_records, trim);
            num_qual_skip += 1;
            num_records += 1;
            continue;
        }

        // check if the read is the right read
        let seq1 = String::from_utf8_lossy(&read1.seq()[0 .. trim]);
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(seq1.len(), wt_read1.len(), -5, -1, &score);
        let alignment = aligner.global(&seq1[8..seq1.len()].as_bytes(), wt_read1);
        if alignment.score < (2 * trim as i32 - 133 - 30) {
            println!("# {} {}: wrong read 1 skipping", num_records, description);
            println!("# {} {}", &seq1[8..seq1.len()], alignment.score);
            num_records += 1;
            num_qual_skip += 1;
            continue;
        }

        if &args[3] == "M" {
            if trim < 33 {
                println!("# {}: Useful read too short for M. Skipping. L = {}", num_records, trim);
                num_qual_skip += 1;
                num_records += 1;
                continue;
            }

            for i in 0 .. cmp::min(trim-6, 110) {
                if &seq1[i .. i+5] == "GCGGC" {
                    match &seq1[i+5 .. i+6] {
                        "A" => am = "WT".to_string(),
                        "G" => am = "AM".to_string(),
                        _ => am = "NA".to_string(),
                    }
                    println!("# am_codon = {}", &seq1[i .. i+6]);
                    break;
                }
            }
        }

        // average quality filtering
        //let avg_qual = read1.qual().iter().fold(0, |a, &b| a as u32 + b as u32);
        //if avg_qual < (125 * 30) { // corresponding to an average quality of 20
        //    println!("# low quality read 1 skipping: {}", avg_qual);
        //    continue;
        //}


        // now deal with read2
        let read2 = record2.unwrap();
        // average quality filtering
        //let avg_qual = read2.qual().iter().fold(0, |a, &b| a as u32 + b as u32);
        //if avg_qual < 125*30 {
        //    println!("# {}: low quality read 2 skipping: {}", num_records, avg_qual);
        //    num_qual_skip += 1;
        //    continue;
        //}
        trim = 124;
        for i in 0..119 {
            if qual_check(&read2.qual()[i .. i+5], &[63, 63, 63, 63, 63]) {
                trim = i+1;
                println!("# {} {}: Read 2 trimmed at {}.", num_records, description, i);
                break;
            }
        }

        if trim < 80 {
            println!("# {}: Useful read too short. Skipping. L = {}", num_records, trim);
            num_qual_skip += 1;
            num_records += 1;
            continue;
        }

        let seq2 =String::from_utf8_lossy(&read2.seq()[0 .. trim]);


        // extract barcodes
        let bc1 = &seq1[0..8];
        let bc2 = &seq2[0..8];
        let bc = bc1.to_string() + bc2;

        // check the pac sequences
        let wt_pac = b"AGAGAAGATTTATCTGAAGTCGTTACGCGAG";
        let seq2_rc = reverse_complement(&seq2);
        let qual : Vec<u8> = read2.qual().iter().cloned().rev().collect();
        let mut aligner = Aligner::with_capacity(wt_pac.len(), seq2.len(), -5, -1, &score);
        let alignment = aligner.local(wt_pac, &seq2_rc.as_bytes());

        let mut pac_start = alignment.ystart;
        if alignment.yend - pac_start < 25 {
            println!("# {} incomplete pac skip.", description);
            num_qual_skip += 1;
            num_records += 1;
            continue;
        }

        match &seq2_rc[pac_start .. pac_start+5] {
            "GAGAA" => pac_start = if pac_start >= 1 {pac_start - 1} else {0},
            "AGAAG" => pac_start = if pac_start >= 2 {pac_start - 2} else {0},
            "GAAGA" => pac_start = if pac_start >= 3 {pac_start - 3} else {0},
            "AAGAT" => pac_start = if pac_start >= 4 {pac_start - 4} else {0},
            "AGATT" => pac_start = if pac_start >= 5 {pac_start - 5} else {0},
            _ => pac_start += 0,
        }
        let pac_end = if pac_start + 31 < trim {pac_start + 31} else { alignment.yend };

        let pac_qual_avg : f32 = qual[pac_start .. pac_end].iter().cloned().map(|x| x as f32).sum::<f32>() / (pac_end - pac_start) as f32;

        let pac = String::from_utf8_lossy(&seq2_rc[pac_start .. pac_end]
                                          .as_bytes()).into_owned();

        if pac_qual_avg < 63.0  {
            println!("# {} {}: pac quality too low ({}).", num_records, description, pac_qual_avg);
            if &args[3] == "M" {
                println!("{} {} {} {}", num_records, bc, pac, am);
            } else {
                println!("{} {} {}", num_records, bc, pac);
            }
            num_records += 1;
            num_qual_skip += 1;
            continue;
        }


        if &args[3] == "M" {
            println!("{} {} {} {}", num_records, bc, pac, am);
        } else {
            println!("{} {} {}", num_records, bc, pac);
        }

        if results.contains_key(&bc) {
            if results[&bc].0 == pac {
                println!("# {}: duplicate found", num_records);
                num_duplicates += 1;
            }
            else {
                println!("# {}: possible sequencing error? {} {}", num_records, &pac, String::from_utf8_lossy(wt_pac));
            }
        }
        else {
            results.insert(bc, (pac, qual.clone()));
        }
        num_records += 1;
    }
    println!("# {} records processed;", num_records);
    println!("# {} low quality reads;", num_qual_skip);
    println!("# {} possible duplicates.", num_duplicates);

    data_stat(&results, &args[4]);
}
