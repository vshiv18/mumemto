/*
 * File: ref_builder.hpp
 * Description: Header for ref_builder.cpp that includes 
 *              definition of the RefBuilder class.
 * Date: August 31, 2022
 */

#ifndef REF_BUILD_H
#define REF_BUILD_H

#include <string>
#include <sdsl/bit_vectors.hpp>

class RefBuilder {
public:
    std::string input_file = "";
    std::string output_ref = "";
    bool use_revcomp = false;

    sdsl::bit_vector doc_ends;
    sdsl::rank_support_v<1> doc_ends_rank;
    size_t num_docs = 0;
    size_t total_length = 0;
    
    RefBuilder(std::string input_data, std::string output_prefix, bool use_rcomp);

}; // end of RefBuilder class


#endif /* end of REF_BUILD_H */