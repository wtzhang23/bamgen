use std::path::*;
use rust_htslib::bam::header::*;
use rust_htslib::bam::Reader;
use rust_htslib::bam::Writer;
use rust_htslib::bam::Record;
use rust_htslib::bam::Format;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Read;
use rand::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::LineWriter;
use clap::Clap;
use bio_types::genome::AbstractInterval;

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
    out_path: PathBuf,
    ground_truth: PathBuf
}

impl Construct {
    fn run(self) {
        // rng
        let mut rng = rand::thread_rng();
        
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

        let mut ground_truth = HashMap::<String, HashMap<usize, usize>>::new();
        let mut record = Record::new();
        for (ri, (len, read_count)) in lengths.into_iter().zip(counts.into_iter()).enumerate() {
            assert!(self.min_ref_len <= len && self.max_ref_len >= len);
            
            let ref_name = format!("ref{}", ri);
            let reference = gen_seq(len, &mut rng);
            record.set_tid(ri as i32);

            // generate batches
            let mut batches = Vec::new();
            let mut count = 0;

            let max = if self.read_depth > read_count {
                0
            } else {
                read_count - self.read_depth
            };

            while count < max {
                let c = (rng.gen::<usize>() % (self.read_depth - 1)) + 1;
                batches.push(c);
                count += c;
            }
            batches.push(read_count - count);

            // generate indices for batch
            let mut indices = rand::seq::index::sample(&mut rng, len - self.min_read_len, batches.len()).into_vec();
            indices.sort();

            let mut umis = Box::new(HashMap::new());
            for (group, pos) in batches.into_iter().zip(indices) {
                record.set_pos(pos as i64);
                for _r in 0..group {
                    let len = rng.gen_range(self.min_read_len, std::cmp::min(len - pos, self.max_read_len));
                    assert!(self.min_read_len <= len && self.max_read_len >= len);

                    let seq = &reference[pos..(pos + len)];

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
                        for _i in 0..(count) {
                            let qname = format!("read{}", read_id);

                            record.set(qname.as_bytes(), None, seq.as_bytes(), &vec![255 as u8; seq.len()]);
                            record.remove_aux(b"OX");
                            record.push_aux(b"OX", &Aux::String(umi.as_bytes()));
                            read_id += 1;

                            writer.write(&record).expect("Failed to write record.");
                        }
                    }
                }

                if let Some(pos_map) = ground_truth.get_mut(&ref_name) {
                    assert!(!pos_map.contains_key(&pos));
                    pos_map.insert(pos, group);
                } else {
                    let mut pos_map = HashMap::new();
                    pos_map.insert(pos, group);
                    ground_truth.insert(ref_name.clone(), pos_map);
                }
                assert_eq!(*ground_truth.get(&ref_name).unwrap().get(&pos).unwrap(), group);
            }
        }

        // write ground truth
        let mut file = File::create(&self.ground_truth).unwrap();
        file.write_all(serde_json::to_string(&ground_truth).unwrap().as_bytes()).unwrap();
    }
}

#[derive(Clap)]
struct Diff {
    output_bam: PathBuf,
    ground_truth: PathBuf,
    diff: PathBuf
}

impl Diff {
    fn run(self) {
        let contents = std::fs::read_to_string(&self.ground_truth).unwrap();
        let ground_truth: HashMap<String, HashMap<usize, usize>> = serde_json::from_str(&contents).unwrap();
        let mut to_compare = HashMap::<String, HashMap<usize, usize>>::new();

        // get predicted counts per position
        let mut reader = Reader::from_path(&self.output_bam).unwrap();
        for record in reader.records().map(|r| r.unwrap()).filter(|r| !r.is_duplicate()) {
            let pos = record.pos() as usize;
            let reference = record.contig();
            if let Some(pos_map) = to_compare.get_mut(reference) {
                if let Some(count) = pos_map.get_mut(&pos) {
                    *count += 1;
                } else {
                    pos_map.insert(pos, 1);
                }
            } else {
                let mut pos_map = HashMap::new();
                pos_map.insert(pos, 1);
                to_compare.insert(reference.to_owned(), pos_map);
            }
        }

        // align predicted w/ ground truths
        let mut aligned_counts = Vec::new();
        for (reference, mut pos_map) in to_compare.drain() {
            for (pos, count) in pos_map.drain() {
                assert!(ground_truth.contains_key(&reference), "Reference {} not found.", reference);
                let corr_pos_map = ground_truth.get(&reference).unwrap();
                assert!(corr_pos_map.contains_key(&pos), "Position {} in Reference {} not found.", pos, reference);
                let corresponding = corr_pos_map.get(&pos).unwrap();
                aligned_counts.push((count, corresponding));
            }
        }

        let file = File::create(&self.diff).unwrap();
        let mut file = LineWriter::new(file);
        file.write_all(format!("{}\t{}\n", "output", "ground truth").as_bytes()).unwrap();

        for (count, gt_count) in aligned_counts.into_iter() {
            file.write_all(format!("{}\t{}\n", count, gt_count).as_bytes()).unwrap();
        }
    }
}

fn main() {
    let opts: Opts = Opts::parse();
    opts.run();
}