/*
 * File: pfp_mum.cpp
 * Description: Structure adapted from docprofiles: https://github.com/oma219/docprofiles/tree/main
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
#include <mem_finder.hpp>
#include <mum_finder.hpp>
#include <getopt.h>
#include <queue>

int build_main(int argc, char** argv, bool mum_mode) {
    /* main method for finding matches */
    if (argc == 1) return mumemto_build_usage();

    // grab the command-line options, and validate them
    BuildOptions build_opts;
    parse_build_options(argc, argv, &build_opts);
    // print_build_status_info(&build_opts);
    build_opts.validate();

    // determine output path for reference, and print all info
    build_opts.output_ref.assign(build_opts.output_prefix + ".fna");
    print_build_status_info(&build_opts, mum_mode);

    if (build_opts.input_list.length() == 0) {build_opts.input_list = make_filelist(build_opts.files, build_opts.output_prefix);}
    // Build the input reference file, and bitvector labeling the end for each doc
    STATUS_LOG("build_main", "building the reference file based on file-list");
    auto start = std::chrono::system_clock::now();
    
    RefBuilder ref_build(build_opts.input_list, build_opts.output_prefix, build_opts.use_rcomp);
    DONE_LOG((std::chrono::system_clock::now() - start));

    // Make sure that document numbers can be store in 2 bytes
    if (ref_build.num_docs >= MAXDOCS)
        FATAL_ERROR("An index cannot be build over %ld documents, "
                    "please reduce to a max of 65,535 docs.", ref_build.num_docs);

    // Determine the paths to the BigBWT executables
    HelperPrograms helper_bins;
    // if (!std::getenv("PFPMUM_BUILD_DIR")) {FATAL_ERROR("Need to set PFPMUM_BUILD_DIR environment variable.");}
    std::filesystem::path path;
    if (!std::getenv("PFPMUM_BUILD_DIR")) {
        path = std::filesystem::canonical("/proc/self/exe").parent_path();
    }
    else {
        path = std::filesystem::path(std::string(std::getenv("PFPMUM_BUILD_DIR")));
    }
    helper_bins.build_paths((path / "bin/").string());
    helper_bins.validate();
    if (!build_opts.from_parse){
        // Parse the input text with BigBWT, and load it into pf object
        STATUS_LOG("build_main", "generating the prefix-free parse for given reference");
        start = std::chrono::system_clock::now();

        run_build_parse_cmd(&build_opts, &helper_bins);
        DONE_LOG((std::chrono::system_clock::now() - start));
    }
    STATUS_LOG("build_main", "building the parse and dictionary objects");
    start = std::chrono::system_clock::now();

    pf_parsing pf(build_opts.output_ref, build_opts.pfp_w);
    DONE_LOG((std::chrono::system_clock::now() - start));

    std::cerr << "\n";

    // Builds the BWT, SA, LCP, and document array profiles and writes to a file

    if (build_opts.missing_genomes >= ref_build.num_docs - 1)
    {
        std::string match_type = mum_mode ? "MUMs" : "MEMs";
        std::string message = "Too few number of sequences, defaulting to multi-" + match_type + " in 2 or more sequences";
        FORCE_LOG("build_main", message.c_str());
        build_opts.missing_genomes = ref_build.num_docs - 2;
    }
    
    start = std::chrono::system_clock::now();

    pfp_lcp lcp(pf, build_opts.output_ref, &ref_build);
    if (mum_mode){
        STATUS_LOG("build_main", "finding multi-MUMs from pfp");
        mum_finder match_finder(build_opts.output_ref, &ref_build, build_opts.min_match_len, build_opts.missing_genomes + 1, build_opts.overlap);
        lcp.process(match_finder);
        match_finder.close();
    }
    else {
        STATUS_LOG("build_main", "finding multi-MEMs from pfp");
        mem_finder match_finder(build_opts.output_ref, &ref_build, build_opts.min_match_len, ref_build.num_docs - build_opts.missing_genomes);
        lcp.process(match_finder);
        match_finder.close();
    }

    lcp.close();

    DONE_LOG((std::chrono::system_clock::now() - start));

    // Print stats before closing out
    FORCE_LOG("build_main", "finished computing matches");
    std::cerr << "\n";
    
    return 0;
}

void run_build_parse_cmd(BuildOptions* build_opts, HelperPrograms* helper_bins) {
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
    // std::cout << "Executing this command: " << command_stream.str().c_str() << std::endl;
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

void print_build_status_info(BuildOptions* opts, bool mum_mode) {
    /* prints out the information being used in the current run */
    std::fprintf(stderr, "\nOverview of Parameters:\n");
    if (opts->input_list.length())
        std::fprintf(stderr, "\tInput file-list: %s\n", opts->input_list.data());
    else {
        std::fprintf(stderr, "\tInput files: ");
        for (int i = 0; i < opts->files.size() - 1; i++)
            {
                std::fprintf(stderr, "%s,", opts->files.at(i).data());
            }
            std::fprintf(stderr, "%s\n", opts->files.at(opts->files.size() - 1).data());
    }
    std::string match_type = mum_mode ? "MUM" : "MEM";
    std::fprintf(stderr, "\tOutput ref path: %s\n", opts->output_ref.data());
    std::fprintf(stderr, "\tPFP window size: %d\n", opts->pfp_w);
    std::fprintf(stderr, "\tMinimum %s length: %d\n", match_type, opts->min_match_len);
    std::fprintf(stderr, "\tInclude reverse complement?: %d\n", opts->use_rcomp);
    std::fprintf(stderr, "\tfinding multi-%ss present in N - %d genomes\n", match_type, opts->missing_genomes);
    if (opts->overlap && opts->missing_genomes > 0 && mum_mode)
        std::fprintf(stderr, "\t\t- including overlapping multi-MUMs\n");
}

std::string make_filelist(std::vector<std::string> files, std::string output_prefix) {
    std::string fname = output_prefix + "_filelist.txt";
    std::ofstream outfile(fname);
    for (size_t i = 0; i < files.size(); ++i) {
        outfile << files[i] << " " << i + 1 << std::endl;
    }
    outfile.close();
    return fname;
}


void parse_build_options(int argc, char** argv, BuildOptions* opts) {
    /* parses the arguments for the build sub-command, and returns a struct with arguments */

    static struct option long_options[] = {
        {"help",      no_argument, NULL,  'h'},
        {"filelist",   required_argument, NULL,  'f'},
        {"output",       required_argument, NULL,  'o'},
        {"revcomp",   no_argument, NULL,  'r'},
        {"missing-genomes",   optional_argument, NULL,  'k'},
        {"no-overlap",   no_argument, NULL,  'p'},
        {"modulus", required_argument, NULL, 'm'},
        {"from-parse",   no_argument, NULL,  's'},
         {"min-match-len",   optional_argument, NULL,  'l'},
        {0, 0, 0,  0}
    };
    int c = 0;
    int long_index = 0;
    while ((c = getopt_long(argc, argv, "hf:o:w:sl:rk:p:m:", long_options, &long_index)) >= 0) {
        switch(c) {
            case 'h': mumemto_build_usage(); std::exit(1);
            case 'f': opts->input_list.assign(optarg); break;
            case 'o': opts->output_prefix.assign(optarg); break;
            case 'w': opts->pfp_w = std::atoi(optarg); break;
            case 'r': opts->use_rcomp = true; break;
            case 'p': opts->overlap = false; break;
            case 'k': opts->missing_genomes = std::atoi(optarg); break;
            case 'm': opts->hash_mod = std::atoi(optarg); break;
            case 's': opts->from_parse = true; break;
            case 'l': opts->min_match_len = std::atoi(optarg); break;
            default: mumemto_build_usage(); std::exit(1);
        }
    }

    for (int i = optind; i < argc; ++i) {
        opts->files.push_back(argv[i]);
    }

}

int mumemto_build_usage() {
    /* prints out the usage information for the build method */
    std::fprintf(stderr, "\nmumemto - find maximal [unique | exact] matches using PFP.\n");
    std::fprintf(stderr, "Usage: mumemto [mum | mem] [options] [input_fasta [...]]\n\n");

    std::fprintf(stderr, "Options:\n");
    std::fprintf(stderr, "\t%-28sprints this usage message\n", "-h, --help");
    std::fprintf(stderr, "\t%-18s%-10spath to a file-list of genomes to use (overrides positional args)\n", "-f, --filelist", "[FILE]");
    std::fprintf(stderr, "\t%-18s%-10soutput prefix path if using -f option\n", "-o, --output", "[arg]");
    std::fprintf(stderr, "\t%-28sinclude the reverse-complement of sequence (default: false)\n\n", "-r, --revcomp");
    std::fprintf(stderr, "\t%-28sminimum MUM or MEM length (default: 20)\n\n", "-l, --min-match-len");

    std::fprintf(stderr, "\t%-28sfind multi-MUMs or multi-MEMs in at least N - k genomes\n\t%-28s(default: 0, match must occur in all sequences, i.e. strict multi-MUM/MEM)\n\n", "-k, --missing-genomes","");
    std::fprintf(stderr, "\t%-28soutput subset multi-MUMs that overlap shorter, more complete multi-MUMs\n\t%-28s(not applicable in MEM mode) (default: true w/ -k)\n", "-p, --no-overlap","");

    std::fprintf(stderr, "PFP options:\n");
    std::fprintf(stderr, "\t%-18s%-10swindow size used for pfp (default: 10)\n", "-w, --window", "[INT]");
    std::fprintf(stderr, "\t%-18s%-10shash-modulus used for pfp (default: 100)\n\n", "-m, --modulus", "[INT]");
    std::fprintf(stderr, "\t%-28suse pre-computed pf-parse\n\n", "-s, --from-parse");

    return 0;
}

int mumemto_usage() {
    /* Prints the usage information for mumemto */
    std::fprintf(stderr, "\nmumemto - find maximal [unique | exact] matches using PFP.\n");
    std::fprintf(stderr, "Usage: mumemto [mum | mem] [options] [input_fasta [...]]\n\n");

    std::fprintf(stderr, "\nmumemto has different modes to run:\n");
    std::fprintf(stderr, "\tmum\tcomputes maximal unique matches (MUMs) in the collection of sequences\n");
    std::fprintf(stderr, "\tmem\tcomputes maximal exact matches (MEMs) in the collection of sequences\n");
    return 0;
}

int main(int argc, char** argv) {
    /* main method for mumemto */
    std::fprintf(stderr, "\033[1m\033[31m\nmumemto version: %s\033[m\033[0m\n", PFPMUM_VERSION);
    bool mum_mode;
    if (argc > 1) {
        if (std::strcmp(argv[1], "mum") == 0)
            mum_mode = true;
        else if (std::strcmp(argv[1], "mem") == 0)
            mum_mode = false;
        else
        {
            std::fprintf(stderr, "\nOne of [mum | mem] mode selection required!\n");
            return mumemto_usage();
        }
        return build_main(argc-1, argv+1, mum_mode);
    }
    return mumemto_usage();
}
