use derive_new::new;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::path::Path;

use mod_kit::mod_base_code::{DnaBase, ModCodeRepr};

pub(crate) mod filenames {
    pub const WEIGHTS_FN: &'static str = "model.mpk";
    pub const CONFIG_FN: &'static str = "model_config.json";
}

#[derive(Deserialize, Serialize, new)]
pub(crate) struct ModelConfiguration {
    pub num_features: usize,
    pub num_classes: usize,
    pub hidden_size: usize,
    pub chunk_size: u32,
    pub modified_bases: HashMap<char, Vec<String>>,
}

impl ModelConfiguration {
    pub(crate) fn modified_bases(&self) -> HashMap<DnaBase, Vec<ModCodeRepr>> {
        self.modified_bases
            .iter()
            .map(|(b, cs)| {
                (
                    DnaBase::parse(*b).unwrap(),
                    cs.iter().map(|x| ModCodeRepr::parse(x).unwrap()).collect(),
                )
            })
            .collect()
    }

    pub(crate) fn get_known_mods(&self) -> Vec<(DnaBase, ModCodeRepr)> {
        self.modified_bases().into_iter().fold(
            Vec::new(),
            |mut acc, (base, codes)| {
                for code in codes {
                    acc.push((base, code));
                }
                acc
            },
        )
    }

    pub(crate) fn from_path<T: AsRef<Path>>(p: T) -> anyhow::Result<Self> {
        let model_config: ModelConfiguration =
            serde_json::from_reader(File::open(p)?)?;
        Ok(model_config)
    }
}

#[derive(Debug)]
pub(crate) enum InferenceMessage<T: Debug> {
    Work(T),
    Done,
}
