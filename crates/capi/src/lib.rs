use std::{mem::transmute, path::PathBuf};

use ggcat_api::{ExtraElaboration, GGCATConfig, GGCATInstance, GeneralSequenceBlockData};

#[repr(transparent)]
struct GGCATInstanceFFI(GGCATInstance);

fn ggcat_create(config: ffi::GGCATConfigFFI) -> &'static GGCATInstanceFFI {
    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: if config.use_temp_dir {
            Some(PathBuf::from(config.temp_dir))
        } else {
            None
        },
        memory: config.memory,
        prefer_memory: config.prefer_memory,
        total_threads_count: config.total_threads_count,
        intermediate_compression_level: if config.intermediate_compression_level != u32::MAX {
            Some(config.intermediate_compression_level)
        } else {
            None
        },
        stats_file: if config.use_stats_file {
            Some(PathBuf::from(config.stats_file))
        } else {
            None
        },
    });
    unsafe { std::mem::transmute(instance) }
}

fn ggcat_build_from_files(
    instance: &'static GGCATInstanceFFI,
    // The input files
    input_files: &[String],

    // The output file
    output_file: String,

    // The names of the colors, ordered by color index
    color_names: &[String],

    // Specifies the k-mers length
    kmer_length: usize,
    // The threads to be used
    threads_count: usize,
    // Treats reverse complementary kmers as different
    forward_only: bool,
    // Overrides the default m-mers (minimizers) length
    minimizer_length: usize,

    // Enable colors
    colors: bool,

    // Minimum multiplicity required to keep a kmer
    min_multiplicity: usize,

    // Extra elaboration step
    extra_elab: usize,
) -> String {
    const EXTRA_ELABORATION_STEP_NONE: usize = 0;
    const EXTRA_ELABORATION_STEP_UNITIG_LINKS: usize = 1;
    const EXTRA_ELABORATION_STEP_GREEDY_MATCHTIGS: usize = 2;
    const EXTRA_ELABORATION_STEP_EULERTIGS: usize = 3;
    const EXTRA_ELABORATION_STEP_PATHTIGS: usize = 4;

    instance
        .0
        .build_graph(
            input_files
                .iter()
                .map(|f| {
                    if f.ends_with(".gfa") {
                        todo!("GFA support is not implemented yet");
                        // GeneralSequenceBlockData::GFA() // PathBuf::from(f))
                    } else {
                        GeneralSequenceBlockData::FASTA(PathBuf::from(f))
                    }
                })
                .collect(),
            PathBuf::from(output_file),
            if color_names.len() > 0 {
                Some(color_names)
            } else {
                None
            },
            kmer_length,
            threads_count,
            forward_only,
            if minimizer_length == usize::MAX {
                None
            } else {
                Some(minimizer_length)
            },
            colors,
            min_multiplicity,
            match extra_elab {
                EXTRA_ELABORATION_STEP_NONE => ExtraElaboration::None,
                EXTRA_ELABORATION_STEP_UNITIG_LINKS => ExtraElaboration::UnitigLinks,
                EXTRA_ELABORATION_STEP_GREEDY_MATCHTIGS => ExtraElaboration::GreedyMatchtigs,
                EXTRA_ELABORATION_STEP_EULERTIGS => ExtraElaboration::Eulertigs,
                EXTRA_ELABORATION_STEP_PATHTIGS => ExtraElaboration::Pathtigs,
                _ => panic!("Invalid extra_elab value: {}", extra_elab),
            },
        )
        .to_str()
        .unwrap()
        .to_string()
}

/// Queries a (optionally) colored graph with a specific set of sequences as queries
fn ggcat_query_graph(
    instance: &'static GGCATInstanceFFI,

    // The input graph
    input_graph: String,
    // The input query as a .fasta file
    input_query: String,

    // The output file
    output_file_prefix: String,

    // Specifies the k-mers length
    kmer_length: usize,
    // The threads to be used
    threads_count: usize,
    // Treats reverse complementary kmers as different
    forward_only: bool,
    // Overrides the default m-mers (minimizers) length
    minimizer_length: usize,

    // Enable colors
    colors: bool,

    // Query output format
    color_output_format: usize,
) -> String {
    const COLORED_QUERY_OUTPUT_FORMAT_JSON_LINES_WITH_NUMBERS: usize = 0;
    const COLORED_QUERY_OUTPUT_FORMAT_JSON_LINES_WITH_NAMES: usize = 1;

    instance
        .0
        .query_graph(
            PathBuf::from(input_graph),
            PathBuf::from(input_query),
            PathBuf::from(output_file_prefix),
            kmer_length,
            threads_count,
            forward_only,
            if minimizer_length == usize::MAX {
                None
            } else {
                Some(minimizer_length)
            },
            colors,
            match color_output_format {
                COLORED_QUERY_OUTPUT_FORMAT_JSON_LINES_WITH_NUMBERS => {
                    ggcat_api::ColoredQueryOutputFormat::JsonLinesWithNumbers
                }
                COLORED_QUERY_OUTPUT_FORMAT_JSON_LINES_WITH_NAMES => {
                    ggcat_api::ColoredQueryOutputFormat::JsonLinesWithNames
                }
                _ => panic!("Invalid color_output_format value: {}", color_output_format),
            },
        )
        .to_str()
        .unwrap()
        .to_string()
}

/// Obtains the standard colormap file path from a graph file path
pub fn ggcat_get_colormap_file(graph_file: String) -> String {
    GGCATInstance::get_colormap_file(PathBuf::from(graph_file))
        .to_str()
        .unwrap()
        .to_string()
}

/// Returns a vector of color names in the given graph.
/// The color indexes returned from the dump_unitigs function
/// can be used to index this vector.
pub fn ggcat_dump_colors(
    // The input colormap
    input_colormap: String,
) -> Vec<String> {
    GGCATInstance::dump_colors(input_colormap).collect()
}

/// Dumps the unitigs of the given graph, optionally with colors
/// It's not guaranteed that maximal unitigs are returned, as only kmers with the same colors subset
/// are returned as whole unitigs to speedup colormap reading times
fn ggcat_dump_unitigs(
    instance: &'static GGCATInstanceFFI,
    // The input graph
    graph_input: String,
    // Specifies the k-mers length
    kmer_length: usize,
    // Overrides the default m-mers (minimizers) length
    minimizer_length: usize,
    // Enable colors
    colors: bool,
    // The threads to be used
    threads_count: usize,

    output_function_context: usize,
    output_function_ptr: usize,
) {
    let output_function: extern "C" fn(usize, usize, usize, usize, usize, bool) =
        unsafe { transmute(output_function_ptr) };

    instance.0.dump_unitigs(
        PathBuf::from(graph_input),
        kmer_length,
        if minimizer_length == usize::MAX {
            None
        } else {
            Some(minimizer_length)
        },
        colors,
        threads_count,
        |sequence, colors, same_colors| {
            output_function(
                output_function_context,
                sequence.as_ptr() as usize,
                sequence.len(),
                colors.as_ptr() as usize,
                colors.len(),
                same_colors,
            );
        },
    )
}

#[cxx::bridge]
mod ffi {
    /// Main config of GGCAT. This config is global and should be passed to GGCATInstance::create
    pub struct GGCATConfigFFI {
        /// If false, a memory only mode is attempted. May crash for large input data if there is no enough RAM memory.
        pub use_temp_dir: bool,

        /// Directory for temporary files
        pub temp_dir: String,

        /// Maximum suggested memory usage (GB)
        /// The tool will try use only up to this GB of memory to store temporary files
        /// without writing to disk. This usage does not include the needed memory for the processing steps.
        /// GGCAT can allocate extra memory for files if the current memory is not enough to complete the current operation
        pub memory: f64,

        /// Use all the given memory before writing to disk
        pub prefer_memory: bool,

        /// The total threads to be used
        pub total_threads_count: usize,

        /// The default lz4 compression level for the intermediate files, -1 to use default values
        pub intermediate_compression_level: u32,

        /// True if the stats file should be created
        pub use_stats_file: bool,
        /// The path to an optional json-formatted real time stats file
        pub stats_file: String,
    }

    extern "Rust" {
        type GGCATInstanceFFI;

        /// Creates a new GGCATInstance. If an instance already exists, it will be returned, ignoring the new config.
        fn ggcat_create(config: GGCATConfigFFI) -> &'static GGCATInstanceFFI;

        /// Builds a new graph from the given input files, with the specified parameters
        fn ggcat_build_from_files(
            instance: &'static GGCATInstanceFFI,
            // The input files
            input_files: &[String],

            // The output file
            output_file: String,

            // The names of the colors, ordered by color index
            color_names: &[String],

            // Specifies the k-mers length
            kmer_length: usize,
            // The threads to be used
            threads_count: usize,
            // Treats reverse complementary kmers as different
            forward_only: bool,
            // Overrides the default m-mers (minimizers) length
            minimizer_length: usize,

            // Enable colors
            colors: bool,

            // Minimum multiplicity required to keep a kmer
            min_multiplicity: usize,

            // Extra elaboration step
            extra_elab: usize,
        ) -> String;

        /// Queries a (optionally) colored graph with a specific set of sequences as queries
        fn ggcat_query_graph(
            instance: &'static GGCATInstanceFFI,

            // The input graph
            input_graph: String,
            // The input query as a .fasta file
            input_query: String,

            // The output file
            output_file_prefix: String,

            // Specifies the k-mers length
            kmer_length: usize,
            // The threads to be used
            threads_count: usize,
            // Treats reverse complementary kmers as different
            forward_only: bool,
            // Overrides the default m-mers (minimizers) length
            minimizer_length: usize,

            // Enable colors
            colors: bool,

            // Query output format
            color_output_format: usize,
        ) -> String;

        fn ggcat_get_colormap_file(graph_file: String) -> String;

        /// Returns a vector of color names in the given graph.
        /// The color indexes returned from the dump_unitigs function
        /// can be used to index this vector.
        pub fn ggcat_dump_colors(
            // The input colormap
            input_colormap: String,
        ) -> Vec<String>;

        /// Dumps the unitigs of the given graph, optionally with colors
        /// It's not guaranteed that maximal unitigs are returned, as only kmers with the same colors subset
        /// are returned as whole unitigs to speedup colormap reading times
        fn ggcat_dump_unitigs(
            instance: &'static GGCATInstanceFFI,
            // The input graph
            graph_input: String,
            // Specifies the k-mers length
            kmer_length: usize,
            // Overrides the default m-mers (minimizers) length
            minimizer_length: usize,
            // Enable colors
            colors: bool,
            // The threads to be used
            threads_count: usize,
            output_function_context: usize,
            // extern "C" fn(context: usize, seq_ptr: usize, seq_len: usize, col_ptr: usize, col_len: usize, same_colors: bool),
            output_function_ptr: usize,
        );
    }
}
