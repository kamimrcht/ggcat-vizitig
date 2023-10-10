use ggcat_api::{
    ColoredQueryOutputFormat, ExtraElaboration, GGCATConfig, GGCATInstance,
    GeneralSequenceBlockData,
};
use itertools::Itertools;
use std::fs::OpenOptions;
use std::io::Write;
use std::{path::PathBuf, sync::Mutex};
use std::io::BufRead;

fn main() {
    let instance = GGCATInstance::create(GGCATConfig {
        temp_dir: Some(PathBuf::from("/tmp")),
        memory: 2.0,
        prefer_memory: true,
        total_threads_count: 16,
        intermediate_compression_level: None,
        stats_file: None,
    });

    let graph_file = PathBuf::from("/tmp/sal-dbg.fa");
    
    // let k = 31;
    let args: Vec<String> = std::env::args().collect();
    let k = &args[1].parse::<usize>().expect("no k value");

    let fof_path = std::env::args().nth(2).expect("no path given");  // Replace this with the path to your file
    let fof = std::fs::File::open(&fof_path).expect("Could not open file");
    let reader = std::io::BufReader::new(fof);
    let in_files: Vec<String> = reader.lines().filter_map(Result::ok).collect();
    let transformed: Vec<GeneralSequenceBlockData> = in_files.iter()
        .map(|s| GeneralSequenceBlockData::FASTA(PathBuf::from(s)))
        .collect();
    let transformed_length = transformed.len();
    let files_strings: Vec<String> = (1..=transformed_length)
            .map(|i| format!("{}", i))
            .collect();
    
    let result: Option<&[String]> = Some(&files_strings[..]);
    let threads_count = 16;

    // Example building of a colored graph from three FASTA files
    // building also bcalm2-style links across maximal unitigs
    let graph_file = instance.build_graph(
        // vec![
        //     GeneralSequenceBlockData::FASTA(PathBuf::from("../../../example-inputs/sal1.fa")),
        //     GeneralSequenceBlockData::FASTA(PathBuf::from("../../../example-inputs/sal2.fa")),
        //     GeneralSequenceBlockData::FASTA(PathBuf::from("../../../example-inputs/sal3.fa")),
        // ],
        transformed,
        graph_file.clone(),
        // Some(&["sal1".to_string(), "sal2".to_string(), "sal3".to_string()]),
        result,
        *k,
        threads_count,
        false,
        None,
        true,
        1,
        ExtraElaboration::UnitigLinks,
    );

    // let input_query = PathBuf::from("../../../example-inputs/query.fa");

    // let output_query = instance.query_graph(
    //     graph_file.clone(),
    //     input_query,
    //     PathBuf::from("/tmp/query-results"),
    //     k,
    //     threads_count,
    //     false,
    //     None,
    //     true,
    //     ColoredQueryOutputFormat::JsonLinesWithNames,
    // );

    // println!("Output query file: {:?}", output_query.display());

    let _print_kmer_lock = Mutex::new(());

    let color_names: Vec<_> =
        GGCATInstance::dump_colors(GGCATInstance::get_colormap_file(&graph_file)).collect();

    instance.dump_unitigs(
        graph_file,
        *k,
        None,
        true,
        threads_count,
        false,
        // WARNING: this function is called asynchronously from multiple threads, so it must be thread-safe.
        // Also the same_colors boolean is referred to the previous call of this function from the current thread
        |read, colors, _same_colors| {

            // let path = "output.txt";
            let path = std::env::args().nth(3).expect("no path given");

            // let mut output_file: File = File::create(path).unwrap();


            // let _lock = print_kmer_lock.lock().unwrap();
            // println!(
            //     ">C:{:?}",
            //     colors.iter().map(|c| &color_names[*c as usize]).format(" ")
            // );
            // println!("{}", std::str::from_utf8(read).unwrap());


            let mut file = OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open(path)
            .unwrap();
            
            let data = format!(">C:{:?}\n{}", colors.iter().map(|c| &color_names[*c as usize]).format(" "), std::str::from_utf8(read).unwrap());

            if let Err(e) = writeln!(file, "{}", data) {
                eprintln!("Couldn't write to file: {}", e);
            }


            
            // output_file.write_all(data.as_bytes()).unwrap();
        },
    );
}
