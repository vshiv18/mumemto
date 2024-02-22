/*
 * File: pfp_mum.hpp
 * Description: Header file for pfp_mum, contains the
 *              function headers and other needed structs 
 *              for pfp_mum.cpp
 *              Heavily adapted from docprofiles by Omar Ahmed: https://github.com/oma219/docprofiles/tree/main
 * Date: December 20th, 2023
 */

#ifndef PFP_MUM_H
#define PFP_MUM_H

#include <string>
#include <iostream>
#include <filesystem>
#include <vector>

/* Useful MACROs */
#define FATAL_ERROR(...) do {std::fprintf(stderr, "\n\033[31mError: \033[m"); std::fprintf(stderr, __VA_ARGS__);\
                              std::fprintf(stderr, "\n\n"); std::exit(1);} while(0)
#define ASSERT(condition, msg) do {if (!condition){std::fprintf(stderr, "\n\n\033[31mAssertion Failed:\033[m %s\n\n", msg); \
                                                   std::exit(1);}} while(0)
#define STATUS_LOG(x, ...) do {std::fprintf(stderr, "\033[32m[%s] \033[0m", x); std::fprintf(stderr, __VA_ARGS__ ); \
                               std::fprintf(stderr, " ... ");} while(0)
#define DONE_LOG(x) do {auto sec = std::chrono::duration<double>(x); \
                        std::fprintf(stderr, "done.  (%.3f sec)\n", sec.count());} while(0)
#define FORCE_LOG(func, ...)  do {std::fprintf(stderr, "\033[32m[%s] \033[m", func); \
                                  std::fprintf(stderr, __VA_ARGS__); \
                                  std::fprintf(stderr, "\n");} while (0)

// Defintions
#define PFPMUM_VERSION "1.0.0"

#define MAXLCPVALUE 65535 // 2^16 - 1
#define MAXDOCS 65535

#define AVX2_PRESENT __AVX2__ 
#define AVX512BW_PRESENT __AVX512BW__ 


/* Function declations */
int mumemto_usage();
int build_main(int argc, char** argv, bool mum_mode);
int mumento_build_usage();
int is_file(std::string path);
int is_dir(std::string path);
std::string  make_filelist(std::vector<std::string> files, std::string output_prefix);
std::vector<std::string> split(std::string input, char delim);
bool is_integer(const std::string& str);
bool endsWith(const std::string& str, const std::string& suffix);
std::string execute_cmd(const char* cmd);


struct BuildOptions {
    public:
        std::string input_list = "";
        std::string output_prefix = "";
        std::string output_ref = "";
        std::vector<std::string> files;
        bool use_rcomp = false;
        size_t pfp_w = 10;
        size_t hash_mod = 100;
        // size_t threads = 0;
        bool is_fasta = true;
        size_t missing_genomes = 0;
        bool overlap = true;
        bool from_parse = false;
        size_t min_match_len = 20;

        void validate() {
            /* checks the arguments and make sure they are valid */
            if (input_list.length() && !is_file(input_list)) // provided a file-list
                FATAL_ERROR("The provided file-list is not valid.");
            else if (input_list.length() && is_file(input_list) && (files.size() > 0)) {
                FORCE_LOG("build_main", "Using filelist, ignoring positional args");
                files.clear();
            }
            else if (input_list.length() == 0 && (files.size() == 0))
                FATAL_ERROR("Need to provide a file-list or files as positional args for processing.");
            
            for (auto f : files) {
                if (!is_file(f)) {
                    FATAL_ERROR(("The following file path is not valid: " + f).c_str());
                }
            }
            std::filesystem::path p (output_prefix);
            if (!is_dir(p.parent_path().string()))
                FATAL_ERROR("Output path prefix is not in a valid directory."); 
   
        }
};

struct HelperPrograms {
  /* Contains paths to run helper programs */
  std::string base_path = "";
  std::string parseNT_bin = "newscanNT.x";
  std::string parse_fasta_bin = "newscan.x";
  std::string parse_bin = "pscan.x";
  
public:
  void build_paths(std::string base) {
      /* Takes the base path, and combines it with names to generate executable paths */
      base_path.assign(base);
      parseNT_bin.assign(base_path + parseNT_bin);
      parse_fasta_bin.assign(base_path + parse_fasta_bin);
      parse_bin.assign(base_path + parse_bin);
  }

  void validate() const {
      /* Makes sure that each path for an executable is valid */
      bool invalid_path = !is_file(parseNT_bin) | !is_file(parse_fasta_bin) | !is_file(parse_bin);
      if (invalid_path) {FATAL_ERROR("One or more of helper program paths are invalid.");}
  }
};

/* Function Declartions involving structs */
void parse_build_options(int argc, char** argv, BuildOptions* opts);
void print_build_status_info(BuildOptions* opts, bool mum_mode);
void run_build_parse_cmd(BuildOptions* build_opts, HelperPrograms* helper_bins);

#endif /* End of PFP_MUM_H */