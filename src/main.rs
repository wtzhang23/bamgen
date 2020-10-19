use std::path::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::LineWriter;
use std::cell::RefCell;
use rust_htslib::bam::header::*;
use rust_htslib::bam::Reader;
use rust_htslib::bam::Writer;
use rust_htslib::bam::Record;
use rust_htslib::bam::Format;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Read;
use rand::prelude::*;
use clap::Clap;
use bio_types::genome::AbstractInterval;
use serde::{Serialize, Deserialize};

#[derive(Clap)]
#[clap(version = "1.0", author = "William Z. <wtzhang23@gmail.com>")]
enum Opts {
    Construct(Construct),
    Diff(Diff),
}

impl Opts {
    fn run(self) {
        match self {
            Opts::Construct(construct) => construct.run(),
            Opts::Diff(diff) => diff.run()
        }
    }
}

#[derive(Serialize, Deserialize, Eq, PartialEq, Hash, Debug)]
struct Loci(String, i64);

#[derive(Serialize, Deserialize)]
struct GroundTruth {
    cpl_prior_duplication: HashMap::<Loci, usize>,
    total_reads: usize,
    reads_per_qname: HashMap::<Vec<u8>, usize>,
}

impl GroundTruth {
    fn new() -> Self {
        Self {
            cpl_prior_duplication: HashMap::new(),
            total_reads: 0,
            reads_per_qname: HashMap::new(),
        }
    }
}

#[derive(Clap)]
struct Construct {
    #[clap(long, default_value="10")]
    n_ref: usize,
    #[clap(long, default_value="20")]
    min_read_len: usize,
    #[clap(long, default_value="30")]
    max_read_len: usize,
    #[clap(long, default_value="50")]
    read_depth: usize,
    #[clap(long, default_value="10000")]
    read_count: usize,
    #[clap(long)]
    paired: bool,
    #[clap(long, default_value=".9")]
    percent_with_umi: f64,
    #[clap(long, default_value="12")]
    umi_len: usize,
    #[clap(long, default_value="5")]
    pcr_steps: usize,
    #[clap(long, default_value=".0001")]
    mut_rate: f32,
    #[clap(long, default_value="1000")]
    min_ref_len: usize,
    #[clap(long, default_value="2000")]
    max_ref_len: usize,
    #[clap(long, default_value="OX")]
    umi_tag: String,
    out_path: PathBuf,
    ground_truth: PathBuf
}

impl Construct {
    fn run(self) {
        // ground truth
        let ground_truth = RefCell::new(GroundTruth::new());

        // rng
        let mut rng = rand::thread_rng();

        let umi_tag = self.umi_tag.as_bytes().to_owned();
        
        // construct header
        let mut header = Header::new();

        let mut lengths = vec![0; self.n_ref];
        for ri in 0..self.n_ref {
            let name = format!("ref{}", ri);
            lengths[ri] = rng.gen_range(self.min_ref_len, self.max_ref_len);
            let length = lengths[ri].to_string();

            let mut record = HeaderRecord::new(b"SQ");
            record.push_tag(b"SN", &name);
            record.push_tag(b"LN", &length);
            header.push_record(&record);
        }

        // get read distribution
        let wumi_count = (self.read_count as f64 * self.percent_with_umi).floor() as usize;

        let mut indices = rand::seq::index::sample(&mut rng, self.read_count - 2, self.n_ref - 1).into_iter().map(|x| x + 1).collect::<Vec<usize>>();
        indices.sort();

        let mut counts = vec![0; self.n_ref];
        let mut last = 0;
        for (i, idx) in indices.into_iter().enumerate() {
            assert!(idx > last);
            counts[i] = idx - last;
            last = idx;
        }
        counts[self.n_ref - 1] = self.read_count - last;

        // writer
        let mut writer = Writer::from_path(self.out_path, &header, Format::BAM).unwrap();
        let mut write = |record: &Record| {
            writer.write(record).expect("Failed to write record.");
            let mut gt = ground_truth.borrow_mut();
            gt.total_reads += 1;
            let qname = record.qname().to_owned();

            // for counts per qname
            if let Some(count) = gt.reads_per_qname.get_mut(&qname) {
                *count += 1;
            } else {
                gt.reads_per_qname.insert(qname, 1);
            }
        };

        // sequence generator
        let sample_nucleotide = |rng: &mut ThreadRng| {
            let base = rng.gen::<usize>() % 4;
            if base == 0 {
                "A".to_owned()
            } else if base == 1 {
                "T".to_owned()
            } else if base == 2 {
                "G".to_owned()
            } else {
                "C".to_owned()
            }
        };

        let gen_seq = |len: usize, rng: &mut ThreadRng| {
            let mut seq = String::new();
            
            for _i in 0..len {
                seq.push_str(&sample_nucleotide(rng));
            }
            seq
        };

        let mut read_id = 0;
        let mut record = Record::new();
        for (ri, (reference_len, read_count)) in lengths.into_iter().zip(counts.into_iter()).enumerate() {
            let has_umi = ri < wumi_count;
            assert!(self.min_ref_len <= reference_len && self.max_ref_len >= reference_len);
            
            // construct reference
            let ref_name = format!("ref{}", ri);
            let reference = gen_seq(reference_len, &mut rng);
            record.set_tid(ri as i32);

            // generate batches
            let mut batches = Vec::new();
            let mut count = 0;

            let max = if self.read_depth > read_count {
                0
            } else {
                read_count - self.read_depth
            };

            // construct batches per loci
            while count < max {
                let c = (rng.gen::<usize>() % (self.read_depth - 1)) + 1;
                batches.push(c);
                count += c;
            }
            batches.push(read_count - count);

            // generate indices for batch
            let mut indices = rand::seq::index::sample(&mut rng, reference_len - self.min_read_len, batches.len()).into_vec();
            indices.sort();

            if has_umi {
                let mut umis = Box::new(HashMap::new());
                for (batch_size, batch_pos) in batches.into_iter().zip(indices) {
                    record.set_pos(batch_pos as i64);

                    for _r in 0..batch_size { // generate batch_size # of reads
                        let read_len = rng.gen_range(self.min_read_len, std::cmp::min(reference_len - batch_pos, self.max_read_len));
                        assert!(self.min_read_len <= read_len && self.max_read_len >= read_len);

                        let read_seq = &reference[batch_pos..(batch_pos + read_len)];

                        umis.insert(gen_seq(self.umi_len, &mut rng), 1);

                        // simulate PCR
                        for _step in 0..self.pcr_steps {
                            let mut new_umis = Box::new(HashMap::<String, usize>::new());
                            for (umi, count) in umis.drain() {
                                for _i in 0..(count << 1) {
                                    let mut new_umi = String::new();
                                    for c in umi.chars() {
                                        let to_mut = rng.gen::<f32>() <= self.mut_rate;
                                        if to_mut {
                                            new_umi.push_str(&sample_nucleotide(&mut rng));
                                        } else {
                                            new_umi.push_str(&c.to_string());
                                        }
                                    }

                                    if new_umis.contains_key(&new_umi) {
                                        let count = new_umis.get(&new_umi).unwrap().to_owned();
                                        new_umis.insert(new_umi, count + 1);
                                    } else {
                                        new_umis.insert(new_umi, 1);
                                    }
                                }
                            }
                            umis = new_umis;
                        }

                        // write reads
                        for (umi, count) in umis.drain() {
                            for _i in 0..count {
                                let qname = format!("read{}", read_id);

                                record.set(qname.as_bytes(), None, read_seq.as_bytes(), &vec![255 as u8; read_seq.len()]);
                                record.remove_aux(&umi_tag);
                                record.push_aux(&umi_tag, &Aux::String(umi.as_bytes()));

                                write(&record);

                                if self.paired { // for pair
                                    let pair_pos = rng.gen_range(0, reference_len - self.min_read_len);
                                    let len = rng.gen_range(self.min_read_len, std::cmp::min(reference_len - pair_pos, self.max_read_len));
                                    let pair_seq = &reference[pair_pos..pair_pos + len];

                                    record.set(qname.as_bytes(), None, pair_seq.as_bytes(), &vec![255 as u8; pair_seq.len()]);
                                    record.remove_aux(&umi_tag);
                                    write(&record);
                                }
                                read_id += 1;
                            }
                        }
                        // write reads per loci
                        let mut gt = ground_truth.borrow_mut();
                        let loci = Loci(ref_name.clone(), batch_pos as i64);
                        gt.cpl_prior_duplication.insert(loci, batch_size);
                    }
                }
            } else {
                for (batch_size, batch_pos) in batches.into_iter().zip(indices) {
                    record.set_pos(batch_pos as i64);

                    for _r in 0..batch_size { // generate batch_size # of reads
                        let read_len = rng.gen_range(self.min_read_len, std::cmp::min(reference_len - batch_pos, self.max_read_len));
                        assert!(self.min_read_len <= read_len && self.max_read_len >= read_len);

                        let read_seq = &reference[batch_pos..(batch_pos + read_len)];
                        let dup_count = 1 << self.pcr_steps; // skip simulating pcr; just write read
                        
                        // write reads
                        for _i in 0..dup_count {
                            let qname = format!("read{}", read_id);

                            record.set(qname.as_bytes(), None, read_seq.as_bytes(), &vec![255 as u8; read_seq.len()]);
                            record.remove_aux(&umi_tag);

                            write(&record);

                            if self.paired { // for pair
                                let pair_pos = rng.gen_range(0, reference_len - self.min_read_len);
                                let len = rng.gen_range(self.min_read_len, std::cmp::min(reference_len - pair_pos, self.max_read_len));
                                let pair_seq = &reference[pair_pos..pair_pos + len];

                                record.set(qname.as_bytes(), None, pair_seq.as_bytes(), &vec![255 as u8; pair_seq.len()]);
                                assert!(record.aux(&umi_tag) == None); // should have been deleted above
                                write(&record);
                            }
                            read_id += 1;
                        }
                    }

                    // write reads per loci
                    let mut gt = ground_truth.borrow_mut();
                    let loci = Loci(ref_name.clone(), batch_pos as i64);
                    gt.cpl_prior_duplication.insert(loci, batch_size);
                }
            }
        }

        // write ground truth
        let mut file = File::create(&self.ground_truth).unwrap();
        file.write_all(ron::to_string(&ground_truth).unwrap().as_bytes()).unwrap();

        println!("Generated {} records\n", ground_truth.borrow().total_reads);
    }
}

#[derive(Clap)]
struct Diff {
    output_bam: PathBuf,
    ground_truth: PathBuf,
    diff: PathBuf,
    #[clap(long)]
    paired: bool,
    #[clap(long, default_value="MI")]
    group_tag: String,
}

impl Diff {
    fn run(self) {
        let contents = std::fs::read_to_string(&self.ground_truth).unwrap();
        let mut ground_truth: GroundTruth = ron::from_str(&contents).unwrap();
        let mut predicted_cpl = HashMap::<Loci, i64>::new();
        let mut num_per_qname = HashMap::<Vec<u8>, usize>::new();
        let mut group_ids = HashMap::<Vec<u8>, Option<i32>>::new();

        let mut reader = Reader::from_path(&self.output_bam).unwrap();
        let mut count = 0;
        for record in reader.records().map(|r| r.unwrap()) {
            count += 1; // will keep track of number of records in file
            
            if let Some(count) = num_per_qname.get_mut(record.qname()) {
                *count += 1;
            } else {
                num_per_qname.insert(record.qname().to_owned(), 1);
            }

            // get predicted counts per position
            if !record.is_duplicate() {
                let loci = Loci(record.contig().to_owned(), record.pos());
                if let Some(count) = predicted_cpl.get_mut(&loci) {
                    *count += 1;
                } else {
                    predicted_cpl.insert(loci, 1);
                }
            }

            // check paired to see if group ids agree
            if self.paired {
                if let Some(group_id) = record.aux(self.group_tag.as_bytes()) {
                    let group_id = match group_id {
                        Aux::String(i) => String::from_utf8(i.to_owned()).unwrap().parse::<i32>().unwrap(),
                        _ => panic!("group id is supposed to be an integer.")
                    };
                    if let Some(other_id) = group_ids.get(record.qname()) {
                        let other_id = other_id.expect("Read w/ same qname had no group id beforehand while this read does not.\n");
                        assert_eq!(other_id, group_id, "Failed to transfer group id to paired read.");
                    } else {
                        group_ids.insert(record.qname().to_owned(), Some(group_id));
                    }
                } else {
                    if let Some(_id) = group_ids.insert(record.qname().to_owned(), None) {
                        panic!("Read q/ same qname had a group id beforehand while this read does not.\n");
                    }
                }
            }
        }
        assert_eq!(count, ground_truth.total_reads, "Deduplication left out/added too much reads.");

        // check that qnames match counts
        for (qname, count) in num_per_qname.drain() {
            let true_qname_count = ground_truth.reads_per_qname.remove(&qname).expect(&format!("qname {} does not exist in ground truth", &String::from_utf8(qname.clone()).unwrap()));
            assert_eq!(count, true_qname_count, "qname {} lost some or gained some reads", &String::from_utf8(qname).unwrap());
        }
        assert!(ground_truth.reads_per_qname.is_empty());

        // align predicted w/ ground truths
        let mut aligned_counts = Vec::new();
        for (loci, count) in predicted_cpl.drain() {
            let true_count = ground_truth.cpl_prior_duplication.remove(&loci).expect(&format!("Loci {:#?} not found.", loci));
            aligned_counts.push((count, true_count, count as i64 - true_count as i64));
        }
        assert!(ground_truth.cpl_prior_duplication.is_empty());

        let file = File::create(&self.diff).unwrap();
        let mut file = LineWriter::new(file);
        file.write_all(format!("{}\t{}\n", "output", "ground truth").as_bytes()).unwrap();

        for (count, gt_count, diff) in aligned_counts.into_iter() {
            file.write_all(format!("{}\t{}\t{}\n", count, gt_count, diff).as_bytes()).unwrap();
        }
    }
}

fn main() {
    let opts: Opts = Opts::parse();
    opts.run();
}