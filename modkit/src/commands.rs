use anyhow::Result as AnyhowResult;
use clap::Subcommand;
use mod_kit::bedmethyl_util::subcommands::EntryBedMethyl;
use mod_kit::dmr::subcommands::BedMethylDmr;
use mod_kit::entropy::subcommand::MethylationEntropy;
use mod_kit::extract::subcommand::ExtractMods;
use mod_kit::localise::subcommand::EntryLocalize;
use mod_kit::modbam_util::subcommands::{
    Adjust, CallMods, EntryModBam, ModSummarize, SampleModBaseProbs, Update,
};
use mod_kit::motifs::subcommand::{EntryFindMotifs, EntryMotifs};
use mod_kit::pileup::subcommand::{DuplexModBamPileup, ModBamPileup};
use mod_kit::repair_tags::RepairTags;
use mod_kit::stats::subcommand::EntryStats;
use mod_kit::validate::subcommand::ValidateFromModBam;

#[derive(Subcommand)]
pub enum Commands {
    /// Tabulates base modification calls across genomic positions. This
    /// command produces a bedMethyl formatted file. Schema and description
    /// of fields can be found in the README.
    Pileup(ModBamPileup),
    /// Performs various operations on BAM files containing base modification
    /// information, such as converting base modification codes and ignoring
    /// modification calls. Produces a BAM output file.
    #[clap(hide = true)]
    AdjustMods(Adjust),
    /// Renames Mm/Ml to tags to MM/ML. Also allows changing the mode flag
    /// from silent '.' to explicitly '?' or '.'.
    #[clap(hide = true)]
    UpdateTags(Update),
    /// Calculate an estimate of the base modification probability
    /// distribution.
    #[clap(hide = true)]
    SampleProbs(SampleModBaseProbs),
    /// Summarize the mod tags present in a BAM and get basic statistics. The
    /// default output is a totals table (designated by '#' lines) and a
    /// modification calls table. Descriptions of the columns can be found
    /// in the README.
    #[clap(hide = true)]
    Summary(ModSummarize),
    /// Call mods from a modbam, creates a new modbam with probabilities set to
    /// 100% if a base modification is called or 0% if called canonical.
    #[clap(hide = true)]
    CallMods(CallMods),
    /// Extract read-level base modification information from a modBAM into a
    /// tab-separated values table.
    #[clap(subcommand)]
    Extract(ExtractMods),
    /// Repair MM and ML tags in one bam with the correct tags from another. To
    /// use this command, both modBAMs _must_ be sorted by read name. The
    /// "donor" modBAM's reads must be a superset of the acceptor's reads.
    /// Extra reads in the donor are allowed, and multiple reads with the
    /// same name (secondary, etc.) are allowed in the acceptor. Reads with
    /// an empty SEQ field cannot be repaired and will be rejected. Reads
    /// where there is an ambiguous alignment of the acceptor to the
    /// donor will be rejected (and logged). See the full documentation for
    /// details.
    #[clap(hide = true)]
    Repair(RepairTags),
    /// Perform DMR test on a set of regions. Output a BED file of regions
    /// with the score column indicating the magnitude of the difference. Find
    /// the schema and description of fields can in the README as well as a
    /// description of the model and method. See subcommand help for
    /// additional details.
    #[clap(subcommand)]
    Dmr(BedMethylDmr),
    /// Tabulates double-stranded base modification patters (such as
    /// hemi-methylation) across genomic motif positions. This command
    /// produces a bedMethyl file, the schema can be found in the online
    /// documentation.
    PileupHemi(DuplexModBamPileup),
    /// Validate results from a set of mod-BAM files and associated BED files
    /// containing the ground truth modified base status at reference
    /// positions.
    Validate(ValidateFromModBam),
    #[clap(hide = true)]
    FindMotifs(EntryFindMotifs),
    /// Various commands to search for, evaluate, or further regine sequence
    /// motifs enriched for base modification. Also can generate BED files of
    /// motif locations.
    #[clap(subcommand)]
    Motif(EntryMotifs),
    /// Use a mod-BAM to calculate methylation entropy over genomic windows.
    Entropy(MethylationEntropy),
    /// Investigate patterns of base modifications, by aggregating pileup
    /// counts "localized" around genomic features of interest.
    #[clap(alias = "localise")]
    Localize(EntryLocalize),
    /// Calculate base modification levels over regions.
    Stats(EntryStats),
    /// Utilities to work with bedMethyl files
    #[clap(subcommand)]
    #[command(name = "bedmethyl", alias = "bm")]
    BedMethyl(EntryBedMethyl),
    /// Utilities to work with modBAM files
    #[clap(subcommand)]
    #[command(name = "modbam", alias = "mb")]
    ModBam(EntryModBam),
    /// Identify regions of open chromatin based on exogenous 6mA signal
    #[clap(subcommand)]
    #[command(name = "open-chromatin", alias = "oc")]
    OpenChromatin(ochm::subcommand::OpenChromatin),
    #[command(hide = true)]
    GenerateShellCompletions {
        /// The shell to generate the completions for
        #[arg(value_enum)]
        shell: clap_complete_command::Shell,
    },
}

impl Commands {
    pub fn run(&self) -> AnyhowResult<()> {
        match self {
            Self::AdjustMods(x) => x.run(),
            Self::Pileup(x) => x.run(),
            Self::SampleProbs(x) => x.run(),
            Self::Summary(x) => x.run(),
            Self::UpdateTags(x) => x.run(),
            Self::CallMods(x) => x.run(),
            Self::Extract(x) => x.run(),
            Self::Repair(x) => x.run(),
            Self::Dmr(x) => x.run(),
            Self::PileupHemi(x) => x.run(),
            Self::Validate(x) => x.run(),
            Self::FindMotifs(x) => x.run(),
            Self::Motif(x) => x.run(),
            Self::Entropy(x) => x.run(),
            Self::Localize(x) => x.run(),
            Self::Stats(x) => x.run(),
            Self::BedMethyl(x) => x.run(),
            Self::ModBam(x) => x.run(),
            Self::OpenChromatin(x) => x.run(),
            Self::GenerateShellCompletions { shell } => {
                use clap::{CommandFactory as _, Parser};

                #[derive(Parser)]
                #[command(name = "modkit")]
                struct Cli {
                    #[command(subcommand)]
                    command: Commands,
                }
                shell.generate(&mut Cli::command(), &mut std::io::stdout());
                Ok(())
            }
        }
    }
}
