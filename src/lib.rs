use ggcat_api::{
    ColoredQueryOutputFormat, ExtraElaboration, GGCATConfig, GGCATInstance,
    GeneralSequenceBlockData,
};
use itertools::Itertools;
use pyo3::prelude::*;
use pyo3::types::PyString;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::Write;
use std::{path::PathBuf, sync::Mutex};

#[pyfunction]
fn build_graph(paths: Vec<&str>, output: &str, k: usize, threads_count: usize) -> PyResult<()> {
    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(PathBuf::from("/tmp")),
        memory: 2.0,
        prefer_memory: true,
        total_threads_count: threads_count,
        intermediate_compression_level: None,
        stats_file: None,
    });
    let graph_file = PathBuf::from("/tmp/sal-dbg.fa");

    let mut transformed: Vec<GeneralSequenceBlockData> = vec![];
    for p in paths.iter() {
        transformed.push(GeneralSequenceBlockData::FASTA(PathBuf::from(p)))
    }
    let cpaths: Vec<_> = (0..paths.len()).map(|i| format!("{}", i)).collect();
    let graph_file = instance.build_graph(
        transformed,
        graph_file,
        Some(&cpaths),
        k,
        threads_count,
        false,
        None,
        true,
        1,
        ExtraElaboration::UnitigLinks,
    );
    let _print_kmer_lock = Mutex::new(());
    let color_names: Vec<_> =
        GGCATInstance::dump_colors(GGCATInstance::get_colormap_file(&graph_file)).collect();
    instance.dump_unitigs(
        graph_file,
        k,
        None,
        true,
        threads_count,
        false,
        // WARNING: this function is called asynchronously from multiple threads, so it must be thread-safe.
        // Also the same_colors boolean is referred to the previous call of this function from the current thread
        |read, colors, _same_colors| {
            let mut file = OpenOptions::new()
                .write(true)
                .append(true)
                .create(true)
                .open(PathBuf::from(output))
                .unwrap();

            let data = format!(
                ">C:{}\n{}",
                colors.iter().map(|c| &color_names[*c as usize]).join(","),
                std::str::from_utf8(read).unwrap()
            );

            if let Err(e) = writeln!(file, "{}", data) {
                eprintln!("Couldn't write to file: {}", e);
            }
        },
    );

    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn ggcat_vizitig(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(build_graph, m)?)?;
    Ok(())
}
