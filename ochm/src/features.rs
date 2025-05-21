use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::ops::ControlFlow;

use anyhow::{anyhow, bail, Context};
use derive_new::new;
use itertools::Itertools;
use ndarray::{Array2, Axis};
use rayon::prelude::*;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

use mod_kit::mod_bam::{BaseModCall, ModBaseInfo};
use mod_kit::mod_base_code::{DnaBase, ModCodeRepr, ParseChar};
use mod_kit::monoid::Moniod;
use mod_kit::threshold_mod_caller::MultipleThresholdModCaller;
use mod_kit::util::{get_query_name_string, record_is_primary, Region, Strand};

// There are 5 "regular" features:
//
const N_BASE_FEATURES: usize = 5;
const A_OFFSET: usize = 0;
const C_OFFSET: usize = 1;
const G_OFFSET: usize = 2;
const T_OFFSET: usize = 3;
const DELETE_OFFSET: usize = 4;

#[derive(Debug)]
pub(crate) struct MedakaFeatures {
    pub matrix: Array2<f32>,
    pub positions: Vec<(u32, MedakaFeatureType)>,
}

#[derive(new, Hash, Eq, PartialEq, Ord, PartialOrd, Debug)]
pub(crate) struct ModifiedPrimaryBase {
    pub dna_base: DnaBase,
    pub mod_code_repr: ModCodeRepr,
}

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Copy, Clone)]
pub enum MedakaFeatureType {
    Major,
    Minor(usize),
}

#[derive(Debug, Hash, Eq, PartialEq, Copy, Clone)]
pub(crate) enum FeatureCode {
    Dna(DnaBase),
    ModBase { dna_base: DnaBase, mod_code: ModCodeRepr },
    // Filtered(DnaBase),
    Delete,
}

#[derive(Debug)]
struct MedakaFeatureAccumulator {
    counts: BTreeMap<
        (u32, MedakaFeatureType),
        FxHashMap<Strand, FxHashMap<FeatureCode, f32>>,
    >,
    n_reads_with_mods: usize,
    total_reads: usize,
}

impl MedakaFeatureAccumulator {
    fn add_delete(&mut self, pos: u32, strand: Strand) {
        self.add_feature(
            pos,
            FeatureCode::Delete,
            1.0f32,
            MedakaFeatureType::Major,
            strand,
        )
    }

    #[inline]
    fn add_feature(
        &mut self,
        pos: u32,
        feature_code: FeatureCode,
        prob: f32,
        medaka_feature_type: MedakaFeatureType,
        strand: Strand,
    ) {
        let k = (pos, medaka_feature_type);
        *self
            .counts
            .entry(k)
            .or_insert(FxHashMap::default())
            .entry(strand)
            .or_insert(FxHashMap::default())
            .entry(feature_code)
            .or_insert(0f32) += prob;
        match feature_code {
            FeatureCode::ModBase { dna_base, .. } => {
                let non_mod_prob = 1f32 - prob;
                self.add_feature(
                    pos,
                    FeatureCode::Dna(dna_base),
                    non_mod_prob,
                    medaka_feature_type,
                    strand,
                );
            }
            _ => {}
        }
    }

    fn encode(
        self,
        code_lookup: &HashMap<ModifiedPrimaryBase, usize>,
    ) -> MedakaFeatures {
        let n_mods = code_lookup.len();
        let n_features = 5usize + 5usize + n_mods;
        let mut matrix: Array2<f32> =
            Array2::zeros((self.counts.len(), n_features));
        let mut positions = Vec::new();

        for (idx, ((pos, feature_type), strand_counts)) in
            self.counts.iter().enumerate()
        {
            positions.push((*pos, *feature_type));
            for (strand, feature_counts) in strand_counts {
                let offset = if *strand == Strand::Positive {
                    0usize
                } else {
                    N_BASE_FEATURES + n_mods
                };
                for (code, count) in feature_counts {
                    let j = match code {
                        FeatureCode::Dna(base) => match base {
                            DnaBase::A => Some(A_OFFSET + offset),
                            DnaBase::C => Some(C_OFFSET + offset),
                            DnaBase::G => Some(G_OFFSET + offset),
                            DnaBase::T => Some(T_OFFSET + offset),
                        },
                        FeatureCode::Delete => Some(DELETE_OFFSET + offset),
                        FeatureCode::ModBase { dna_base, mod_code } => {
                            if let Some(idx) = code_lookup.get(
                                &ModifiedPrimaryBase::new(*dna_base, *mod_code),
                            ) {
                                Some(*idx)
                            } else {
                                match dna_base {
                                    DnaBase::A => Some(A_OFFSET + offset),
                                    DnaBase::C => Some(C_OFFSET + offset),
                                    DnaBase::G => Some(G_OFFSET + offset),
                                    DnaBase::T => Some(T_OFFSET + offset),
                                }
                            }
                        }
                    };
                    if let Some(j) = j {
                        matrix[[idx, j]] = *count;
                    }
                }
            }
        }

        MedakaFeatures { matrix, positions }
    }

    fn into_contiguous(
        mut self,
        code_lookup: &HashMap<ModifiedPrimaryBase, usize>,
    ) -> Vec<MedakaFeatures> {
        if self.counts.is_empty() {
            vec![]
        } else {
            let mut acc = Vec::new();
            let mut current = BTreeMap::new();
            let (k, dat) = self.counts.pop_first().unwrap();
            let mut pointer = k.0;
            current.insert(k, dat);

            while let Some((k, dat)) = self.counts.pop_first() {
                let pos = k.0;
                if pos.checked_sub(pointer).expect("should be increasing") > 1 {
                    let mut next_chunk = BTreeMap::new();
                    next_chunk.insert(k, dat);
                    let done = std::mem::replace(&mut current, next_chunk);
                    acc.push(done);
                } else {
                    current.insert(k, dat);
                }
                pointer = pos;
            }

            if !current.is_empty() {
                acc.push(current);
            }

            acc.into_par_iter()
                .map(|counts| Self {
                    counts,
                    n_reads_with_mods: 0,
                    total_reads: 0,
                })
                .map(|x| x.encode(code_lookup))
                .collect()
        }
    }

    fn set_read_count(&mut self, mod_calls_read_cache: &ModCallsReadCache) {
        self.n_reads_with_mods = mod_calls_read_cache.cache.len();
        self.total_reads =
            self.n_reads_with_mods + mod_calls_read_cache.skip_list.len();
    }
}

impl Moniod for MedakaFeatureAccumulator {
    fn zero() -> Self {
        Self { counts: BTreeMap::new(), n_reads_with_mods: 0, total_reads: 0 }
    }

    fn op(self, other: Self) -> Self {
        let mut this = self;
        this.op_mut(other);
        this
    }

    fn op_mut(&mut self, other: Self) {
        for (k, v) in other.counts.into_iter() {
            self.counts.insert(k, v);
        }
        self.n_reads_with_mods += other.n_reads_with_mods;
        self.total_reads += other.total_reads;
    }

    fn len(&self) -> usize {
        self.counts.len()
    }
}

struct ModCallsReadCache<'a> {
    cache: FxHashMap<String, FxHashMap<usize, BaseModCall>>,
    caller: &'a MultipleThresholdModCaller,
    skip_list: FxHashSet<String>,
}

impl<'a> ModCallsReadCache<'a> {
    fn new(caller: &'a MultipleThresholdModCaller) -> Self {
        Self {
            cache: FxHashMap::default(),
            skip_list: FxHashSet::default(),
            caller,
        }
    }

    fn add_record(&mut self, record: &bam::Record, record_name: &str) {
        match ModBaseInfo::new_from_record(record) {
            Ok(info) => {
                let mut read_pos_to_calls = FxHashMap::default();
                for (base, probs) in info.pos_seq_base_mod_probs.iter() {
                    for (pos, prob) in probs.pos_to_base_mod_probs.iter() {
                        let call = self.caller.call(base, &prob);
                        assert!(
                            read_pos_to_calls.insert(*pos, call).is_none(),
                            "multiple calls for a read pos"
                        );
                    }
                }
                self.cache.insert(record_name.to_owned(), read_pos_to_calls);
            }
            Err(_e) => {
                self.skip_list.insert(record_name.to_string());
            }
        }
    }

    fn get_mod_calls(
        &mut self,
        record: &bam::Record,
        record_name: &str,
        forward_read_position: usize,
    ) -> Option<BaseModCall> {
        if self.skip_list.contains(record_name) {
            None
        } else {
            if let Some(calls) = self.cache.get(record_name) {
                calls.get(&forward_read_position).copied()
            } else {
                self.add_record(record, record_name);
                assert!(
                    self.cache.contains_key(record_name)
                        || self.skip_list.contains(record_name)
                );
                self.get_mod_calls(record, record_name, forward_read_position)
            }
        }
    }
}

#[derive(new)]
struct PileupIter<'a> {
    pileups: bam::pileup::Pileups<'a, bam::IndexedReader>,
    start_pos: u32,
    end_pos: u32,
}

impl<'a> Iterator for PileupIter<'a> {
    type Item = bam::pileup::Pileup;

    fn next(&mut self) -> Option<Self::Item> {
        let mut pileup: Option<Self::Item> = None;
        while let Some(Ok(plp)) = self.pileups.next() {
            let pos = plp.pos();
            let off_end = pos >= self.end_pos;
            if off_end {
                return None;
                // we're done
            } else if pos < self.start_pos {
                // advance into region we're looking at
                continue;
            } else {
                pileup = Some(plp);
                break;
            }
        }
        pileup
    }
}

pub(crate) fn get_features_for_region(
    region: &Region,
    reader: &mut bam::IndexedReader,
    caller: &MultipleThresholdModCaller,
    mod_lookup: &HashMap<ModifiedPrimaryBase, usize>,
) -> anyhow::Result<Vec<MedakaFeatures>> {
    get_features(region, reader, caller, &mod_lookup)
        .map(|acc| acc.into_contiguous(&mod_lookup))
}

fn get_features(
    region: &Region,
    reader: &mut bam::IndexedReader,
    caller: &MultipleThresholdModCaller,
    known_mods: &HashMap<ModifiedPrimaryBase, usize>,
) -> anyhow::Result<MedakaFeatureAccumulator> {
    let header = reader.header();
    let tid = header.tid(region.name.as_bytes()).ok_or_else(|| {
        anyhow!("BAM header does not contain {}", region.name)
    })?;
    reader.fetch(FetchDefinition::Region(
        tid as i32,
        region.start as i64,
        region.end as i64,
    ))?;

    let pileup_iter: PileupIter =
        PileupIter::new(reader.pileup(), region.start, region.end);

    let mut features = MedakaFeatureAccumulator::zero();
    let mut cache = ModCallsReadCache::new(caller);
    for plp in pileup_iter {
        let position = plp.pos();
        let alignments = plp
            .alignments()
            // no ref-skip
            .filter(|aln| !aln.is_refskip())
            // take only primary
            .filter(|aln| record_is_primary(&aln.record()));

        for alignment in alignments {
            // todo remove unwrap here
            let record_name =
                get_query_name_string(&alignment.record()).unwrap();
            let strand = if alignment.record().is_reverse() {
                Strand::Negative
            } else {
                Strand::Positive
            };
            if alignment.is_del() {
                // info!("{record_name}, {position}, DELETE");
                features.add_delete(position, strand);
                continue;
            }
            let q_pos = alignment.qpos().unwrap();
            let forward_seq_position = if alignment.record().is_reverse() {
                alignment
                    .qpos()
                    .and_then(|p| alignment.record().seq_len().checked_sub(p))
                    .and_then(|p| p.checked_sub(1))
                    .unwrap()
            } else {
                q_pos
            };

            let get_read_base = |p: usize| -> DnaBase {
                DnaBase::parse_char(alignment.record().seq()[p] as char)
                    .unwrap()
            };

            let read_base = get_read_base(q_pos);
            let (major_feature, p) = get_feature(
                &record_name,
                &alignment.record(),
                &mut cache,
                read_base,
                forward_seq_position,
                known_mods,
            );

            features.add_feature(
                position,
                major_feature,
                p,
                MedakaFeatureType::Major,
                strand,
            );

            match alignment.indel() {
                Indel::Ins(ins_size) => {
                    for offset in (0..ins_size).map(|x| (x as usize) + 1) {
                        let query_position = if alignment.record().is_reverse()
                        {
                            forward_seq_position.checked_sub(offset).unwrap()
                        } else {
                            forward_seq_position.saturating_add(offset)
                        };
                        let read_base =
                            get_read_base(q_pos.saturating_add(offset));
                        let (feature, p) = get_feature(
                            &record_name,
                            &alignment.record(),
                            &mut cache,
                            read_base,
                            query_position,
                            known_mods,
                        );
                        features.add_feature(
                            position,
                            feature,
                            p,
                            MedakaFeatureType::Minor(offset),
                            strand,
                        );
                    }
                }
                _ => {}
            }
        }
    }
    features.set_read_count(&cache);
    Ok(features)
}

#[inline]
fn get_feature(
    record_name: &str,
    record: &bam::Record,
    mod_calls_read_cache: &mut ModCallsReadCache,
    read_base: DnaBase,
    forward_read_position: usize,
    known_mods: &HashMap<ModifiedPrimaryBase, usize>,
) -> (FeatureCode, f32) {
    // let debug = record_name == "d1780f84-9c61-4897-986c-9296135fbdba";
    match mod_calls_read_cache.get_mod_calls(
        record,
        record_name,
        forward_read_position,
    ) {
        Some(BaseModCall::Canonical(_)) => {
            (FeatureCode::Dna(read_base), 1.0f32)
        }
        Some(BaseModCall::Modified(p, code)) => {
            if known_mods
                .contains_key(&ModifiedPrimaryBase::new(read_base, code))
            {
                (
                    FeatureCode::ModBase {
                        dna_base: read_base,
                        mod_code: code,
                    },
                    p,
                )
            } else {
                (FeatureCode::Dna(read_base), 1.0f32)
            }
        }
        Some(BaseModCall::Filtered) => {
            unreachable!("should never have filtered")
        }
        None => (FeatureCode::Dna(read_base), 1f32),
    }
}

pub(crate) fn get_code_lookup_for_mods(
    known_mods: &[(DnaBase, ModCodeRepr)],
) -> HashMap<ModifiedPrimaryBase, usize> {
    known_mods
        .iter()
        .flat_map(|(primary_base, mod_code)| {
            [
                ModifiedPrimaryBase::new(*primary_base, *mod_code),
                ModifiedPrimaryBase::new(primary_base.complement(), *mod_code),
            ]
        })
        .collect::<BTreeSet<ModifiedPrimaryBase>>()
        .into_iter()
        .sorted_by(|a, b| a.dna_base.cmp(&b.dna_base))
        .enumerate()
        .map(|(i, k)| (k, i + 5usize)) // there are 5 "regular" features
        .collect::<HashMap<ModifiedPrimaryBase, usize>>()
}

// pub type LabeledCountsFeatures = CountsFeatures<u64>;
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CountsFeatures<L: Serialize + Copy> {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub counts: Vec<f32>,
    pub coverages: Vec<u32>,
    pub n_features: usize,
    pub width: usize,
    pub label: L,
}

impl<T: Serialize + Copy> CountsFeatures<T> {
    pub fn check_has_informative_positions(
        self,
        code_lookup: &HashMap<ModifiedPrimaryBase, usize>,
    ) -> anyhow::Result<Self> {
        let offset = code_lookup.len() + 5usize;
        let indices_to_check = code_lookup
            .keys()
            .flat_map(|modified_base| {
                let i = match &modified_base.dna_base {
                    DnaBase::A => 0usize,
                    DnaBase::C => 1usize,
                    DnaBase::G => 2usize,
                    DnaBase::T => 3usize,
                };
                [i, i + offset]
            })
            .chain(code_lookup.values().copied())
            .collect::<Vec<usize>>();
        let has_informative_positions =
            self.counts.chunks(self.n_features).any(|chunk| {
                indices_to_check.iter().map(|i| chunk[*i]).sum::<f32>() > 0f32
            });
        if has_informative_positions {
            Ok(self)
        } else {
            bail!(
                "chunk {}:{}-{} has doesn't have any informative positions",
                self.chrom,
                self.start,
                self.end
            )
        }
    }

    pub fn has_min_coverage_at_all_positions(
        self,
        min_coverage: u32,
    ) -> anyhow::Result<Self> {
        if self.coverages.iter().all(|c| *c >= min_coverage) {
            Ok(self)
        } else {
            bail!(
                "chunk {}:{}-{} does not have minimum coverage at all \
                 positions",
                self.chrom,
                self.start,
                self.end
            )
        }
    }

    pub fn from_medaka_major_features(
        region: &Region,
        medaka_features: MedakaFeatures,
        chunk_size: u64,
        label: T,
    ) -> anyhow::Result<Self> {
        let n_features = medaka_features.matrix.shape()[1];
        let major_matrix = medaka_features
            .matrix
            .rows()
            .into_iter()
            .zip(medaka_features.positions.iter())
            .filter_map(|(row, (_, feature_type))| match feature_type {
                MedakaFeatureType::Major => Some(row.to_vec()),
                MedakaFeatureType::Minor(_) => None,
            })
            .collect::<Vec<_>>();
        let n_positions = major_matrix.len();
        let feature_matrix = Array2::from_shape_vec(
            [n_positions, n_features],
            major_matrix.into_iter().flatten().collect(),
        )
        .context("failed creation of filtered array")?;

        if feature_matrix.sum() == 0f32 {
            bail!("zero coverage over region {region}")
        }

        let coverages = feature_matrix.sum_axis(Axis(1));

        let mut feature_mat = feature_matrix.mapv(|x| x);
        feature_mat.rows_mut().into_iter().zip(coverages.iter()).for_each(
            |(mut row, &cov)| {
                if cov > 0f32 {
                    row /= cov;
                }
            },
        );
        assert!(
            !feature_mat.iter().any(|x| x.is_nan()),
            "detected NaN in features"
        );

        let n_features = feature_mat.ncols();
        let width = feature_mat.nrows();
        let counts = feature_mat.into_raw_vec();

        if width < chunk_size as usize {
            bail!("width < chunk size")
        }

        let coverages = coverages.mapv(|x| x as u32).to_vec();
        Ok(Self {
            chrom: region.name.clone(),
            start: region.start,
            end: region.end,
            coverages,
            counts,
            n_features,
            width,
            label,
        })
    }

    pub(crate) fn into_chunked(
        self,
        chunk_size: u32,
        overlap: u32,
    ) -> Vec<Self> {
        assert_eq!(self.counts.len() % self.n_features, 0);
        let chunk_size = chunk_size as usize;
        let overlap = overlap as usize;
        let feature_width = chunk_size * self.n_features;
        let step_size_bp = chunk_size.checked_sub(overlap).unwrap();
        let step_size = step_size_bp * self.n_features;
        let chunks = (0..self.counts.len())
            .step_by(step_size)
            .enumerate()
            .try_fold(Vec::new(), |mut agg, (i, start)| {
                let end = start + feature_width;
                if end <= self.counts.len() {
                    let counts = self.counts[start..end].to_vec();
                    let cov_start = i * step_size_bp;
                    let cov_end = cov_start + chunk_size;
                    let coverages = self.coverages[cov_start..cov_end].to_vec();

                    let start = self.start as usize + (step_size_bp * i);
                    let end = start + chunk_size;
                    let chunk = Self {
                        chrom: self.chrom.clone(),
                        start: start as u32,
                        end: end as u32,
                        counts,
                        n_features: self.n_features,
                        coverages,
                        width: chunk_size,
                        label: self.label,
                    };
                    agg.push(chunk);
                    ControlFlow::Continue(agg)
                } else {
                    ControlFlow::Break(agg)
                }
            });
        match chunks {
            ControlFlow::Continue(x) => x,
            ControlFlow::Break(x) => x,
        }
    }
}
