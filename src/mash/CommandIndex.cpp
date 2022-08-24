// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#include "CommandIndex.h"
#include "Sketch.h"
#include "sketchParameterSetup.h"
#include <iostream>

#include <chopper/layout/configuration.hpp>
#include <chopper/layout/data_store.hpp>
#include <chopper/layout/hierarchical_binning.hpp>

using std::cerr;
using std::endl;
using std::string;
using std::vector;

namespace mash {

CommandIndex::CommandIndex()
: Command()
{
    name = "Index";
    summary = "Create Indexes (reduced representations for fast operations).";
    description = "Create a Index file, which is a reduced representation of a sequence or set of sequences (based on min-hashes) that can be used for fast distance estimations. Inputs can be fasta or fastq files (gzipped or not), and \"-\" can be given to read from standard input. Input files can also be files of file names (see -l). For output, one Index file will be generated, but it can have multiple Indexes within it, divided by sequences or files (see -i). By default, the output file name will be the first input file with a '.msh' extension, or 'stdin.msh' if standard input is used (see -o).";
    argumentString = "<input> [<input>] ...";

    useOption("help");
    addOption("list", Option(Option::Boolean, "l", "Input", "List input. Lines in each <input> specify paths to sequence files, one per line.", ""));
    addOption("prefix", Option(Option::File, "o", "Output", "Output prefix (first input file used if unspecified). The suffix '.msh' will be appended.", ""));
    // addOption("id", Option(Option::File, "I", "Index", "ID field for Index of reads (instead of first sequence ID).", ""));
    // addOption("comment", Option(Option::File, "C", "Index", "Comment for a Index of reads (instead of first sequence comment).", ""));
    // addOption("counts", Option(Option::Boolean, "M", "Index", "Store multiplicity of each k-mer in each Index.", ""));
    useIndexOptions();
}

int CommandIndex::run() const
{
    if ( arguments.size() == 0 || options.at("help").active )
    {
        print();
        return 0;
    }

    int verbosity = 1;//options.at("silent").active ? 0 : options.at("verbose").active ? 2 : 1;
    bool list = options.at("list").active;

    for ( int i = 0; i < arguments.size(); i++ )
    {
        if ( false && hasSuffix(arguments[i], suffixIndex) )
        {
            cerr << "ERROR: " << arguments[i] << " looks like it is already Indexed." << endl;
            exit(1);
        }
    }

    string prefix;

    if ( options.at("prefix").argument.length() > 0 )
    {
        prefix = options.at("prefix").argument;
    }
    else
    {
        if ( arguments[0] == "-" )
        {
            prefix = "stdin";
        }
        else
        {
            prefix = arguments[0];
        }
    }

    // -------------------------------------------------------------------------
    // load mash sketch
    // -------------------------------------------------------------------------
	Sketch sketch;
    Sketch::Parameters parameters;
    vector<string> refArgVector;
    refArgVector.push_back(arguments[0]);
    sketch.initFromFiles(refArgVector, parameters);

    // fill chopper data for layouting
    chopper::layout::data_store chopper_data;

    // load mash sketch file
    string hll_sketches_dir{prefix + "_sketches"};
    std::filesystem::create_directory(hll_sketches_dir);

    for ( int i = 0; i < sketch.getReferenceCount(); i++ )
	{
		const HashList & hashes = sketch.getReference(i).hashesSorted;

        // compute and write hll-sketch to file for chopper-layout
        chopper::sketch::hyperloglog sketch(12/* default sketch_bits of chopper*/);
        for (auto k_hash : hashes)
            sketch.add(reinterpret_cast<char*>(&(k_hash->first)), sizeof((k_hash->first)));
        std::filesystem::path path = hll_sketches_dir / std::to_string(i);
        path += ".hll";
        std::ofstream hll_fout(path, std::ios::binary);
        sketch.dump(hll_fout);

        // fill layout data
		chopper_data.filenames.push_back(std::to_string(i));
        chopper_data.kmer_sizes.push_back(hashes.size());
	}

    // -------------------------------------------------------------------------
    // Compute Chopper Layout
    // -------------------------------------------------------------------------
    chopper::layout::configuration config;

    config.output_filename = prefix + ".layout";
    config.prefix = prefix;
    config.sketch_directory = hll_sketches_dir;
    config.tmax = std::sqrt(sketch.getReferenceCount());
    config.num_hash_functions = 4;
    config.false_positive_rate = 0.05;
    donfig.alpha = 1.2;
    donfig.max_rearrangement_ratio =  0.5;
    donfig.threads = 8;
    donfig.estimate_union = true;
    donfig.rearrange_user_bins = true;
    donfig.determine_best_tmax = false;
    donfig.force_all_binnings = false;
    donfig.output_verbose_statistics = false;
    donfig.debug = false;

    chopper_data.compute_fp_correction(config.false_positive_rate, config.num_hash_functions, config.tmax);

    std::stringstream output_buffer;
    std::stringstream header_buffer;

    chopper_data.output_buffer = &output_buffer;
    chopper_data.header_buffer = &header_buffer;
    chopper_data.false_positive_rate = config.false_positive_rate;

    size_t max_hibf_id;

    chopper::layout::hibf_statistics global_stats{config, chopper_data.fp_correction, chopper_data.kmer_counts};
    chopper_data.stats = &global_stats.top_level_ibf;
    size_t dummy{};

    max_hibf_id = chopper::layout::hierarchical_binning{chopper_data, config}.execute(); // compute layout

    if (config.output_verbose_statistics)
    {
        global_stats.print_header();
        global_stats.print_summary(dummy);
    }

    // brief Write the output to the layout file.
    {
        std::ofstream fout{config.output_filename};
        write_layout_header_to(config, max_hibf_id, header_buffer.str(), fout);
        fout << output_buffer.str();
    }

    // collect input files
    // vector<string> files;
    // for ( int i = 0; i < arguments.size(); i++ )
    // {
    //     if ( list )
    //     {
    //         splitFile(arguments[i], files);
    //     }
    //     else
    //     {
    //         files.push_back(arguments[i]);
    //     }
    // }

    // -------------------------------------------------------------------------
    // Build index on layout
    // -------------------------------------------------------------------------



    std::string index_file{prefix + ".indec"}
    cerr << "Writing to " << index_file << "..." << endl;

    // write index

    return 0;
}

} // namespace mash
