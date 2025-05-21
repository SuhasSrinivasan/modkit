use std::path::PathBuf;

use anyhow::{anyhow, bail, Context};
use log::{debug, error};
use rayon::prelude::*;
use rust_htslib::bam;

use mod_kit::mod_base_code::{DnaBase, ModCodeRepr};
use mod_kit::threshold_mod_caller::MultipleThresholdModCaller;
use mod_kit::util::Region;

use crate::features;
use crate::features::CountsFeatures;
use crate::util::{InferenceMessage, ModelConfiguration};

pub(crate) struct FeatureLoader {
    bam: PathBuf,
    caller: MultipleThresholdModCaller,
    known_mods: Vec<(DnaBase, ModCodeRepr)>,
    chunk_size: u32,
    overlap: u32,
    min_coverage: u32,
}

impl FeatureLoader {
    pub(crate) fn new(
        bam: &PathBuf,
        caller: MultipleThresholdModCaller,
        model_config: &ModelConfiguration,
        width: u32,
        step_size: u32,
        min_coverage: u32,
    ) -> anyhow::Result<Self> {
        let overlap = width.checked_sub(step_size).ok_or_else(|| {
            anyhow!("illegal step size, must be less than {width}")
        })?;
        if width <= overlap {
            bail!("width cannot be <= overlap")
        }
        bam::IndexedReader::from_path(bam)
            .context("failed to create indexed reader")?;
        Ok(Self {
            bam: bam.clone(),
            caller,
            known_mods: model_config.get_known_mods(),
            chunk_size: width,
            overlap,
            min_coverage,
        })
    }
    pub(crate) fn run(
        &self,
        rcv: crossbeam::channel::Receiver<InferenceMessage<Vec<Region>>>,
        sender: crossbeam::channel::Sender<
            InferenceMessage<Vec<CountsFeatures<()>>>,
        >,
    ) {
        let code_lookup = features::get_code_lookup_for_mods(&self.known_mods);
        loop {
            match rcv.recv() {
                Ok(InferenceMessage::Work(regions)) => regions.into_par_iter().for_each(|region| {
                    if let Ok(mut reader) = bam::IndexedReader::from_path(&self.bam) {
                        let features = features::get_features_for_region(
                            &region,
                            &mut reader,
                            &self.caller,
                            &code_lookup,
                        )
                            .map(|medaka_features| {
                                medaka_features
                                    .into_iter()
                                    .filter_map(|f| {
                                        match CountsFeatures::from_medaka_major_features(
                                            &region,
                                            f,
                                            self.chunk_size as u64,
                                            (),
                                        ) {
                                            Ok(counts_features) => Some(counts_features),
                                            Err(e) => {
                                                debug!("failed to make medaka features for region {region}, {e}");
                                                None
                                            }
                                        }
                                    })
                                    .collect::<Vec<_>>()
                            })
                            .map(|cs| {
                                cs.into_iter()
                                    .flat_map(|x| {
                                        x.into_chunked(self.chunk_size, self.overlap)
                                    })
                                    .filter_map(|x| {
                                        let x = x.check_has_informative_positions(&code_lookup)
                                            .and_then(|x| x.has_min_coverage_at_all_positions(self.min_coverage));
                                        match x {
                                            Ok(f) => Some(f),
                                            Err(_e) => {
                                                // todo may want to remove this here and just filter when writing..?
                                                //  also may want to skip when not sufficient coverage?
                                                // debug!("skipping chunk, {e}");
                                                None
                                            }
                                        }
                                    })
                                    .collect::<Vec<_>>()
                            });
                        match features {
                            Ok(xs) => {
                                sender.send(InferenceMessage::Work(xs)).unwrap()
                            }
                            Err(e) => {
                                debug!("failed to make features, {e}");
                            }
                        }
                    } else {
                        error!("failed to get reader");
                    }
                }),
                Ok(InferenceMessage::Done) => {
                    debug!("feature gen received Done, propagating..");
                    sender.send(InferenceMessage::Done).unwrap()
                }
                Err(e) => {
                    debug!("Feature loader got {e}, stopping.");
                    break;
                }
            }
        }
        drop(sender);
    }
}
