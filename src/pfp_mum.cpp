/*
 * File: pfp_mum.cpp
 * Description: Heavily adapted from docprofiles: https://github.com/oma219/docprofiles/tree/main
 *              Main file for building the prefix-free parse, and outputing MUMs along with r-index data structures (BWT/SA/LCP). 
 *              This workflow is based on doc profiles by Omar Ahmed, 
 *              which is in turn based on the pfp_lcp.hpp developed by Massimiliano Rossi.
 *
 * Date: December 20th, 2023
 */

#include <pfp_mum.hpp> 
#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <filesystem>
#include <ref_builder.hpp>
#include <vector>
#include <pfp.hpp>
#include <pfp_lcp_mum.hpp>
#include <getopt.h>
#include <queue>

int build_main(int argc, char** argv) {
    /* main method for build the document profiles */
    if (argc == 1) return pfpdoc_build_usage();

    // grab the command-line options, and validate them
    PFPDocBuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    build_opts.validate();

    // determine output path for reference, and print all info
    build_opts.output_ref.assign(build_opts.output_prefix + ".fna");
    print_build_status_info(&build_opts);

    // Build the input reference file, and bitvector labeling the end for each doc
    STATUS_LOG("build_main", "building the reference file based on file-list");
    auto start = std::chrono::system_clock::now();
    
    RefBuilder ref_build(build_opts.input_list, build_opts.output_prefix, build_opts.use_rcomp);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // Make sure that document numbers can be store in 2 bytes, and 
    // it makes sense with respect to k
    if (ref_build.num_docs >= MAXDOCS)
        FATAL_ERROR("An index cannot be build over %ld documents, "
                    "please reduce to a max of 65,535 docs.", ref_build.num_docs);
    // if ((build_opts.use_taxcomp || build_opts.use_topk) 
    //     && ref_build.num_docs < build_opts.numcolsintable)
    //     FATAL_ERROR("the k provided is larger than the number of documents.");

    // Check and make sure the document to extract is a valid id
    // if (build_opts.doc_to_extract > ref_build.num_docs) {
    //     FATAL_ERROR("Document #%d was requested to be extracted,"
    //                 " but there are only %d documents", build_opts.doc_to_extract, ref_build.num_docs);
    // }

    // Determine the paths to the BigBWT executables
    HelperPrograms helper_bins;
    if (!std::getenv("PFPDOC_BUILD_DIR")) {FATAL_ERROR("Need to set PFPDOC_BUILD_DIR environment variable.");}
    helper_bins.build_paths((std::string(std::getenv("PFPDOC_BUILD_DIR")) + "/bin/").data());
    helper_bins.validate();

    // Parse the input text with BigBWT, and load it into pf object
    STATUS_LOG("build_main", "generating the prefix-free parse for given reference");
    start = std::chrono::system_clock::now();

    run_build_parse_cmd(&build_opts, &helper_bins);
    DONE_LOG((std::chrono::system_clock::now() - start));

    STATUS_LOG("build_main", "building the parse and dictionary objects");
    start = std::chrono::system_clock::now();

    pf_parsing pf(build_opts.output_ref, build_opts.pfp_w);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // Print info regarding the compression scheme being used
    std::cerr << "\n";
    // if (build_opts.use_taxcomp)
    //     FORCE_LOG("build_main", "taxonomic compression of the doc profiles will be used");
    // else if (build_opts.use_topk)
    //     FORCE_LOG("build_main", "top-k compression of the doc profile will be used");
    // else   
    //     FORCE_LOG("build_main", "no compression scheme will be used for the doc profiles");

    // Builds the BWT, SA, LCP, and document array profiles and writes to a file
    STATUS_LOG("build_main", "building bwt and doc profiles based on pfp");
    start = std::chrono::system_clock::now();

    pfp_lcp lcp(pf, build_opts.output_ref, &ref_build, build_opts.missing_genomes);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // Print stats before closing out
    FORCE_LOG("build_main", "finished building");
    std::cerr << "\n";
    
    return 0;
}

void run_build_parse_cmd(PFPDocBuildOptions* build_opts, HelperPrograms* helper_bins) {
    // Generates and runs the command-line for executing the PFP of the reference 
    std::ostringstream command_stream;
    // if (build_opts->threads > 0) {
    //     std::string curr_exe = "";
    //     if (build_opts->is_fasta) {curr_exe.assign(helper_bins->parse_fasta_bin);}
    //     else {curr_exe.assign(helper_bins->parse_bin);}

    //     command_stream << curr_exe << " "; //<< " -i ";
    //     command_stream << build_opts->output_ref << " ";
    //     command_stream << "-w " << build_opts->pfp_w;
    //     command_stream << " -p " << build_opts->hash_mod;
    //     // command_stream << " -t " << build_opts->threads;
    // }
    // else {
        std::string curr_exe = "";
        command_stream << helper_bins->parseNT_bin << " "; // << " -i ";
        command_stream << build_opts->output_ref << " ";
        command_stream << "-w " << build_opts->pfp_w;
        command_stream << " -p " << build_opts->hash_mod;
    // }
    if (build_opts->is_fasta) {command_stream << " -f";}

    // std::cout << command_stream.str() << std::endl;
    //LOG(build_opts->verbose, "build_parse", ("Executing this command: " + command_stream.str()).data());
    auto parse_log = execute_cmd(command_stream.str().c_str());
    //OTHER_LOG(parse_log.data());
}

std::string execute_cmd(const char* cmd) {
    std::array<char, 256> buffer{};
    std::string output = "";

    std::string cmd_plus_stderr = std::string(cmd) + " 2>&1";
    FILE* pipe = popen(cmd_plus_stderr.data(), "r"); // Extract stderr as well
    if (!pipe) {FATAL_ERROR("popen() failed!");}

    try {
        std::size_t bytes;
        while ((bytes = fread(buffer.data(), sizeof(char), sizeof(buffer), pipe))) {
            output += std::string(buffer.data(), bytes);
        }
    } catch (...) {
        pclose(pipe);
        FATAL_ERROR("Error occurred while reading popen() stream.");
    }

    // Check if the command failed
    size_t exit_code = WEXITSTATUS(pclose(pipe));
    if (exit_code != 0) {
        std::cout << "\n";
        FATAL_ERROR("external command failed ... here is the error message:\n%s", output.data());
    }

    return output;
}

bool endsWith(const std::string& str, const std::string& suffix) {
    // Checks if the string ends the suffix
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

bool is_integer(const std::string& str) {
    /* Checks if string passed is an integer */
    std::string::const_iterator iter = str.begin();
    while (iter != str.end() && std::isdigit(*iter)) {++iter;}
    return !str.empty() && iter == str.end();
}

std::vector<std::string> split(std::string input, char delim) {
    /* Takes in a string, and splits it based on delimiters */
    std::vector<std::string> word_list;
    std::string curr_word = "";

    for (char ch: input) {
        if (ch == delim && curr_word.length()) {word_list.push_back(curr_word); curr_word = "";}
        else if (ch != delim) {curr_word += ch;}
    }
    if (curr_word.length()) {word_list.push_back(curr_word);}
    return word_list;
}

int is_file(std::string path) {
    /* Checks if the path is a valid file-path */
    std::ifstream test_file(path.data());
    if (test_file.fail()) {return 0;}
    
    test_file.close();
    return 1;
}

int is_dir(std::string path) {
    /* Checks if the directory is a valid path */
    return std::filesystem::exists(path);
}

void print_build_status_info(PFPDocBuildOptions* opts) {
    /* prints out the information being used in the current run */
    std::fprintf(stderr, "\nOverview of Parameters:\n");
    std::fprintf(stderr, "\tInput file-list: %s\n", opts->input_list.data());
    std::fprintf(stderr, "\tOutput ref path: %s\n", opts->output_ref.data());
    std::fprintf(stderr, "\tPFP window size: %d\n", opts->pfp_w);
    std::fprintf(stderr, "\tInclude rev-comp?: %d\n", opts->use_rcomp);
    std::fprintf(stderr, "\tfinding multi-MUMs present in N - %d genomes\n", opts->missing_genomes);
    // std::fprintf(stderr, "\tUse heuristics?: %d\n\n", opts->use_heuristics);
}

void parse_build_options(int argc, char** argv, PFPDocBuildOptions* opts) {
    /* parses the arguments for the build sub-command, and returns a struct with arguments */

    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"filelist",   required_argument, NULL,  'f'},
        {"output",       required_argument, NULL,  'o'},
        {"revcomp",   no_argument, NULL,  'r'},
        {"missing-genomes",   optional_argument, NULL,  'k'},
        // {"taxcomp",   no_argument, NULL,  't'},
        // {"num-col",   required_argument, NULL,  'k'},
        // {"top-k",   no_argument, NULL,  'p'},
        // {"print-doc", required_argument, NULL, 'e'},
        // {"no-heuristic", no_argument, NULL, 'n'},
        {"modulus", required_argument, NULL, 'm'},
        {0, 0, 0,  0}
    };

    int c = 0;
    int long_index = 0;
    while ((c = getopt_long(argc, argv, "hf:o:w:rtk:pe:nm:", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': pfpdoc_build_usage(); std::exit(1);
            case 'f': opts->input_list.assign(optarg); break;
            case 'o': opts->output_prefix.assign(optarg); break;
            case 'w': opts->pfp_w = std::atoi(optarg); break;
            case 'r': opts->use_rcomp = true; break;
            // case 't': opts->use_taxcomp = true; break;
            // case 'p': opts->use_topk = true; break;
            case 'k': opts->missing_genomes = std::atoi(optarg); break;
            // case 'e': opts->doc_to_extract = std::atoi(optarg); break;
            // case 'n': opts->use_heuristics = false; break;
            case 'm': opts->hash_mod = std::atoi(optarg); break;
            default: pfpdoc_build_usage(); std::exit(1);
        }
    }
}

int pfpdoc_build_usage() {
    /* prints out the usage information for the build method */
    std::fprintf(stderr, "\npfp_mum build - find mums using PFP.\n");
    std::fprintf(stderr, "Usage: pfp_mum build [options]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-28sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-18s%-10spath to a file-list of genomes to use\n", "-f, --filelist", "[FILE]");
    std::fprintf(stderr, "\t%-18s%-10soutput prefix path if using -f option\n", "-o, --output", "[arg]");
    std::fprintf(stderr, "\t%-28sinclude the reverse-complement of sequence (default: false)\n\n", "-r, --revcomp");

    std::fprintf(stderr, "\t%-28sfind multi-MUMs in at least N - k genomes (default: 0, strict multi-MUM)\n\n", "-k, --missing-genomes");

    // std::fprintf(stderr, "\t%-28suse taxonomic compression of the document array (default: false)\n", "-t, --taxcomp");
    // std::fprintf(stderr, "\t%-28suse top-k compression of the document array (default: false)\n", "-p, --top-k");
    // std::fprintf(stderr, "\t%-18s%-10snumber of columns to include in the main table (default: 7)\n\n", "-k, --num-col", "[INT]");
    
    // std::fprintf(stderr, "\t%-18s%-10sdocument number whose profiles to extract\n", "-e, --print-doc", "[INT]");
    // std::fprintf(stderr, "\t%-28sturn off any heuristics used to build profiles\n\n", "-n, --no-heuristic");

    std::fprintf(stderr, "\t%-18s%-10swindow size used for pfp (default: 10)\n", "-w, --window", "[INT]");
    std::fprintf(stderr, "\t%-18s%-10shash-modulus used for pfp (default: 100)\n\n", "-m, --modulus", "[INT]");

    return 0;
}

int pfpdoc_usage() {
    /* Prints the usage information for pfp_doc */
    std::fprintf(stderr, "\npfp_mum has different sub-commands to run:\n");
    std::fprintf(stderr, "Usage: pfp_mum <sub-command> [options]\n\n");

    std::fprintf(stderr, "Commands:\n");
    std::fprintf(stderr, "\tbuild\tbuilds BWT/SA/LCP, and computes mums\n");
    // std::fprintf(stderr, "\trun\truns queries with respect to the document array structure\n");
    // std::fprintf(stderr, "\tinfo\tprint out information regarding this index and document array\n\n");
    return 0;
}

int main(int argc, char** argv) {
    /* main method for pfp_doc */
    std::fprintf(stderr, "\033[1m\033[31m\npfp-mum version: %s\033[m\033[0m\n", PFPDOC_VERSION);

    if (argc > 1) {
        if (std::strcmp(argv[1], "build") == 0) 
            return build_main(argc-1, argv+1);
    }
    return pfpdoc_usage();
}