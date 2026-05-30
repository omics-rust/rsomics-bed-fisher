use std::fs::File;
use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_fisher::{fisher, write_fisher};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-fisher", disable_help_flag = true)]
pub struct Cli {
    /// Query BED file A (required).
    #[arg(short = 'a', long = "query", required = true)]
    pub query: PathBuf,

    /// Reference BED file B (required).
    #[arg(short = 'b', long = "reference", required = true)]
    pub reference: PathBuf,

    /// Genome sizes file (chrom\tlength per line) (required).
    #[arg(short = 'g', long = "genome", required = true)]
    pub genome: PathBuf,

    /// Merge overlapping intervals in A and B before testing.
    #[arg(short = 'm', long = "merge")]
    pub merge: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let fa = File::open(&self.query).map_err(RsomicsError::Io)?;
        let fb = File::open(&self.reference).map_err(RsomicsError::Io)?;
        let fg = File::open(&self.genome).map_err(RsomicsError::Io)?;
        let result = fisher(fa, fb, fg, self.merge)?;
        let stdout = io::stdout();
        let mut out = stdout.lock();
        write_fisher(&result, &mut out)
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Compute Fisher's exact test for interval overlap significance (bedtools fisher equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["-a <QUERY> -b <REFERENCE> -g <GENOME> [OPTIONS]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('a'),
                long: "query",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Query BED file (A)",
                why_default: None,
            },
            FlagSpec {
                short: Some('b'),
                long: "reference",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Reference BED file (B)",
                why_default: None,
            },
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Genome sizes file (chrom TAB length)",
                why_default: None,
            },
            FlagSpec {
                short: Some('m'),
                long: "merge",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Merge overlapping intervals before testing",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Test overlap significance between ChIP peaks and gene promoters",
        command: "rsomics-bed-fisher -a peaks.bed -b promoters.bed -g hg38.genome",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
