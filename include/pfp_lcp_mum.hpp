/* pfp_lcp - lcp from prefix free parsing 
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file pfp_lcp.hpp
   \brief pfp_lcp.hpp define and build the lcp from the prefix-free parsing.
            Adapted by Vikram Shivakumar to compute Maximal Unique Matches (MUM) between sequences (12/20/2023)
   \author Massimiliano Rossi
   \date 01/07/2020
*/


#ifndef _LCP_PFP_HH
#define _LCP_PFP_HH

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C"
{
#include <gsacak.h>
}

#include <ref_builder.hpp>
#include <pfp.hpp>

#include <unordered_map>
#include <queue>
#include <boost/circular_buffer.hpp>

class pfp_lcp{
public:

    pf_parsing& pf;
    std::vector<size_t> min_s; // Value of the minimum lcp_T in each run of BWT_T
    std::vector<size_t> pos_s;    // Position of the minimum lcp_T in each run of BWT_T

    const size_t MIN_MUM_LENGTH = 20;

    uint8_t head;
    size_t length = 0; // Length of the current run of BWT_T
    size_t num_docs = 0;
    // std::vector<uint8_t> heads;


    pfp_lcp(pf_parsing &pfp_, std::string filename, RefBuilder* ref_build) : 
                pf(pfp_),
                min_s(1, pf.n),
                pos_s(1,0),
                num_docs(ref_build->num_docs),
                head(0),
                bwt_window(num_docs),
                doc_window(num_docs),
                lcp_window(num_docs),
                sa_window(num_docs)
                // heads(1, 0)
    {
        // Opening output files
        std::string outfile = filename + std::string(".lcp");
        if ((lcp_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".ssa");
        if ((ssa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".esa");
        if ((esa_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".bwt");
        if ((bwt_file = fopen(outfile.c_str(), "w")) == nullptr)
            error("open() file " + outfile + " failed");

        outfile = filename + std::string(".mums");
        mum_file.open(outfile);

        //instantiate doc tracker
        window_docs = new int[num_docs]();

        assert(pf.dict.d[pf.dict.saD[0]] == EndOfDict);

        phrase_suffix_t curr;
        phrase_suffix_t prev;

        inc(curr);
        while (curr.i < pf.dict.saD.size())
        {

            if(is_valid(curr)){
                // Compute the next character of the BWT of T
                std::vector<phrase_suffix_t> same_suffix(1, curr);  // Store the list of all phrase ids with the same suffix.

                phrase_suffix_t next = curr;

                while (inc(next) && (pf.dict.lcpD[next.i] >= curr.suffix_length))
                {

                    assert(next.suffix_length >= curr.suffix_length);
                    assert((pf.dict.b_d[next.sn] == 0 && next.suffix_length >= pf.w) || (next.suffix_length != curr.suffix_length));
                    if (next.suffix_length == curr.suffix_length)
                    {
                        same_suffix.push_back(next);
                    }

                }

                // Hard case
                int_t lcp_suffix = compute_lcp_suffix(curr,prev);

                typedef std::pair<int_t *, std::pair<int_t *, uint8_t>> pq_t;

                // using lambda to compare elements.
                auto cmp = [](const pq_t &lhs, const pq_t &rhs) {
                    return *lhs.first > *rhs.first;
                };

                std::priority_queue<pq_t, std::vector<pq_t>, decltype(cmp)> pq(cmp);
                for (auto s: same_suffix)
                {
                    size_t begin = pf.pars.select_ilist_s(s.phrase + 1);
                    size_t end = pf.pars.select_ilist_s(s.phrase + 2);
                    pq.push({&pf.pars.ilist[begin], {&pf.pars.ilist[end], s.bwt_char}});
                }

                size_t prev_occ;
                bool first = true;
                while (!pq.empty())
                {
                    auto curr_occ = pq.top();
                    pq.pop();

                    if (!first)
                    {
                        // Compute the minimum s_lcpP of the the current and previous occurrence of the phrase in BWT_P
                        lcp_suffix = curr.suffix_length + min_s_lcp_T(*curr_occ.first, prev_occ);
                    }
                    first = false;
                    // Update min_s
                    print_lcp(lcp_suffix, j);

                    update_ssa(curr, *curr_occ.first);

                    update_bwt(curr_occ.second.second, 1);

                    update_esa(curr, *curr_occ.first);

                    // Start of MUM computation code
                    uint8_t curr_bwt_ch = curr_occ.second.second;
                    right_lcp = lcp_suffix;
                    size_t sa_i = ssa;
                    size_t doc_i = ref_build->doc_ends_rank(ssa);

                    bool valid_window = sa_window.size() == num_docs;

                    if(valid_window && is_mum())
                    {
                        // std::cout << "MUM FOUND!" << std::endl;
                        write_mum();
                        // skip = num_docs - 1;
                    }
                    // else if (skip > 0)
                    // {
                    //     skip--;
                    //     total_skips++;
                    // }
                    
                    // update the window datastructures (i.e. slide over)
                    update_lcp_window(right_lcp, valid_window);
                    update_sa_window(sa_i, valid_window);
                    update_bwt_window(curr_bwt_ch, valid_window);
                    update_doc_window(doc_i, valid_window);


                    // assert(bwt_window.size() == num_docs);
                    // assert(doc_window.size() == num_docs);
                    // assert(sa_window.size() == num_docs);
                    // assert(lcp_window.size() == num_docs);

                    // End of MUM computation code

                    
                    // Update prevs
                    prev_occ = *curr_occ.first;

                    // Update pq
                    curr_occ.first++;
                    if (curr_occ.first != curr_occ.second.first)
                        pq.push(curr_occ);

                    j += 1;
                }

                prev = same_suffix.back();
                
                curr = next;
                
            }
            else
            {
                inc(curr);
            }
        }
        // print last BWT char and SA sample
        print_sa();
        print_bwt();

        // Close output files
        fclose(ssa_file);
        fclose(esa_file);
        fclose(bwt_file);
        fclose(lcp_file);

        delete[] window_docs;

        mum_file.close();
        // std::cout << "Skipped checking " << total_skips << " entries (" << (total_skips * 100.0 / j) << "%)";
    }

private:
    typedef struct
    {
        size_t i = 0; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t phrase = 0;
        size_t suffix_length = 0;
        int_da sn = 0;
        uint8_t bwt_char = 0;
    } phrase_suffix_t;

    
    size_t j = 0;

    size_t ssa = 0;
    size_t esa = 0;


    // Helper functions and variables to compute MUMs
    //
    // variables needed to track mums
    size_t left_lcp = 0; // lcp of window and suffix preceding window
    size_t right_lcp = 0; // lcp of window and suffix succeeding window

    // represent RMQ lcp of window as an ordered_set
    // std::multiset<size_t> rmq_window;
    

    // try the linear RMQ algorithm
    std::deque<std::pair<size_t, size_t>> lcp_pq;
    
    // try circular buffers instead of deques!
    boost::circular_buffer<size_t> bwt_window;
    boost::circular_buffer<size_t> doc_window;
    boost::circular_buffer<size_t> lcp_window;
    boost::circular_buffer<size_t> sa_window;

    // for window_docs, define an update that decrements the outgoing and increments the incoming doc
    // track docs present in window using frequency map, removing keys if 0
    // std::unordered_map<size_t, int> window_docs_vec;
    int* window_docs;
    size_t total_unique_docs = 0;

    // track BWT chars present in window using frequency map, removing keys if 0
    // checks if MUM is left extendable (set of BWT chars is the same)
    // std::unordered_map<uint8_t, int> window_bwt_vec;
    int window_bwt[4] = {0};
    size_t total_unique_bwt = 0;
    const std::map<uint8_t,int> nucMap = {
        {'A', 0},
        {'C', 1},
        {'G', 2},
        {'T', 3},
        {0, 4} //null char
    };

    // stores the current RMQ (avoiding recomputing)
    size_t mum_length = 0;

    // min number of iterations to avoid checking for a mum
    // int skip = 0;
    // int total_skips = 0;

    inline void add_doc(size_t d)
    {
        if(window_docs[d])
        {
            window_docs[d]++;
        }
        else
        {
            window_docs[d] = 1;
            total_unique_docs++;
        }
    }
    inline void remove_doc(size_t d)
    {
        assert(window_docs[d] > 0);
        if(window_docs[d] == 1)
        {
            window_docs[d]--;
            total_unique_docs--;
        }
        else
        {
            window_docs[d]--;
        }
    }
    inline void add_bwt(uint8_t bwt_c)
    {
        int idx = nucMap.at(bwt_c);
        if(window_bwt[idx])
        {
            window_bwt[idx]++;
        }
        else
        {
            window_bwt[idx] = 1;
            total_unique_bwt++;
        }
    }
    inline void remove_bwt(uint8_t bwt_c)
    {
        int idx = nucMap.at(bwt_c);
        
        assert(window_bwt[idx] > 0);
        if(window_bwt[idx] == 1)
        {
            window_bwt[idx]--;
            total_unique_bwt--;
        }
        else
        {
            window_bwt[idx]--;
        }
    }

    inline void update_bwt_window(uint8_t bwt_c, bool valid_window)
    {
        add_bwt(bwt_c);
        if(valid_window)
        {
            remove_bwt(bwt_window.front());
        }
        bwt_window.push_back(bwt_c);
    }

    inline void update_doc_window(size_t doc, bool valid_window)
    {
        add_doc(doc);
        if(valid_window)
        {
            remove_doc(doc_window.front());
        }
        doc_window.push_back(doc);
        // if(num_docs - window_docs.size() > skip)
        //     skip = num_docs - window_docs.size();
    }

    inline void update_sa_window(size_t sa_entry, bool valid_window)
    {
        // slide over sa window
        sa_window.push_back(sa_entry);
    }

    inline void update_lcp_window(size_t lcp, bool valid_window)
    {
        // get next lcp from queue, and remove it from set too
        // add lcp to the queue to maintain window order
        lcp_window.push_back(lcp);
        left_lcp = lcp_window.front();

        // update pq
        while(!lcp_pq.empty() && lcp_pq.front().second <= (j + 1 - num_docs))
            lcp_pq.pop_front();
        while(!lcp_pq.empty() && lcp_pq.back().first > lcp)
            lcp_pq.pop_back();
        lcp_pq.push_back(std::pair<size_t, size_t>(lcp, j));
    }

    inline size_t rmq_of_window()
    {
        // get min LCP in window from ordered set
        // *rmq_window.begin();
        // return *std::min_element(std::next(lcp_window.begin()), lcp_window.end());
        // if(lcp_pq.front().first != *std::min_element(std::next(lcp_window.begin()), lcp_window.end()))
        // {
        //     std::cout << "Position: "<< j << ", Actual min: " << (*std::min_element(std::next(lcp_window.begin()), lcp_window.end())) << std::endl;
        //     for(auto it = lcp_pq.begin(); it != std::prev(lcp_pq.end()); ++it)
        //     {
        //         std::cout << "(" << (*it).first << "," << (*it).second << ")" << ',';
        //     }
        //     std::cout << "(" << lcp_pq.back().first << "," << lcp_pq.back().second << ")" << std::endl;
        //     std::cout << "actual window: ";
        //     for(auto it = lcp_window.begin(); it != std::prev(lcp_window.end()); ++it)
        //     {
        //         std::cout << *it << ',';
        //     }
        //     std::cout << lcp_window.back() << std::endl;
        // }
        return lcp_pq.front().first;
    }
    inline bool is_mum()
    {
        // Check each condition: (check the fast conditions first, then compute RMQ if needed)
        // Check that every doc appears once
        if(total_unique_docs != num_docs)
            return false;
        // Check BWT chars in that range are not all identical (i.e. can be left extended by 1, not maximal)
        if(total_unique_bwt == 1)
            return false;
        // check RMQ LCP of window > min_mum
        mum_length = rmq_of_window();
        // long enough MUM
        if(mum_length < MIN_MUM_LENGTH)
            return false;
        // Suffix preceding and succeding window don't share long enough prefix (mum is not unique!)
        if(left_lcp >= mum_length || right_lcp >= mum_length)
            return false;
        
        return true;

    }

    inline void write_mum()
    {
        mum_file << std::to_string(mum_length) << '\t';
        for(auto it = sa_window.begin(); it != std::prev(sa_window.end()); ++it)
        {
            mum_file << *it << ',';
        }
        mum_file << std::to_string(sa_window.back()) << std::endl;
        // std::cout << std::to_string(mum_length) << '\t';
        // std::ostream_iterator<size_t> output_iterator(std::cout, ",");
        // std::copy(sa_window.begin(), std::prev(sa_window.end()), output_iterator);
        // std::cout << std::to_string(sa_window.back()) << std::endl;        
    }


    FILE *lcp_file;

    FILE *bwt_file;

    FILE *ssa_file;
    FILE *esa_file;

    std::ofstream mum_file;

    inline bool inc(phrase_suffix_t& s)
    {
        s.i++;
        if (s.i >= pf.dict.saD.size())
            return false;
        s.sn = pf.dict.saD[s.i];
        s.phrase = pf.dict.rank_b_d(s.sn);
        // s.phrase = pf.dict.daD[s.i] + 1; // + 1 because daD is 0-based
        assert(!is_valid(s) || (s.phrase > 0 && s.phrase < pf.pars.ilist.size()));
        s.suffix_length = pf.dict.select_b_d(pf.dict.rank_b_d(s.sn + 1) + 1) - s.sn - 1;
        if(is_valid(s))
            s.bwt_char = (s.sn == pf.w ? 0 : pf.dict.d[s.sn - 1]);
        return true;
    }

    inline bool is_valid(phrase_suffix_t& s)
    {
        // avoid the extra w # at the beginning of the text
        if (s.sn < pf.w)
            return false;
        // Check if the suffix has length at least w and is not the complete phrase.
        if (pf.dict.b_d[s.sn] != 0 || s.suffix_length < pf.w)
            return false;
        
        return true;
    }
    
    inline int_t min_s_lcp_T(size_t left, size_t right)
    {
        // assume left < right
        if (left > right)
            std::swap(left, right);

        assert(pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] >= pf.w);

        return (pf.s_lcp_T[pf.rmq_s_lcp_T(left + 1, right)] - pf.w);
    }

    inline int_t compute_lcp_suffix(phrase_suffix_t& curr, phrase_suffix_t& prev)
    {
        int_t lcp_suffix = 0;

        if (j > 0)
        {
            // Compute phrase boundary lcp
            lcp_suffix = pf.dict.lcpD[curr.i];
            for (size_t k = prev.i + 1; k < curr.i; ++k)
            {
                lcp_suffix = std::min(lcp_suffix, pf.dict.lcpD[k]);
            }

            if (lcp_suffix >= curr.suffix_length && curr.suffix_length == prev.suffix_length)
            {
                // Compute the minimum s_lcpP of the phrases following the two phrases
                // we take the first occurrence of the phrase in BWT_P
                size_t left = pf.pars.ilist[pf.pars.select_ilist_s(curr.phrase + 1)]; //size_t left = first_P_BWT_P[phrase];
                // and the last occurrence of the previous phrase in BWT_P
                size_t right = pf.pars.ilist[pf.pars.select_ilist_s(prev.phrase + 2) - 1]; //last_P_BWT_P[prev_phrase];
                
                lcp_suffix += min_s_lcp_T(left,right);
            }
        }

        return lcp_suffix;
    }

    inline void print_lcp(int_t val, size_t pos)
    {
        size_t tmp_val = val;
        if (fwrite(&tmp_val, THRBYTES, 1, lcp_file) != 1)
            error("LCP write error 1");
    }

    // We can put here the check if we want to store the LCP or stream it out
    inline void new_min_s(int_t val, size_t pos)
    {
        min_s.push_back(val);
        pos_s.push_back(j);
    }

    inline void update_ssa(phrase_suffix_t &curr, size_t pos)
    {
        ssa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(ssa < (pf.n - pf.w + 1ULL));
    }

    inline void update_esa(phrase_suffix_t &curr, size_t pos)
    {
        esa = (pf.pos_T[pos] - curr.suffix_length) % (pf.n - pf.w + 1ULL); // + pf.w;
        assert(esa < (pf.n - pf.w + 1ULL));
    }

    inline void print_sa()
    {
        if (j < (pf.n - pf.w + 1ULL))
        {
            size_t pos = j;
            if (fwrite(&pos, SSABYTES, 1, ssa_file) != 1)
                error("SA write error 1");
            if (fwrite(&ssa, SSABYTES, 1, ssa_file) != 1)
                error("SA write error 2");
        }

        if (j > 0)
        {
            size_t pos = j - 1;
            if (fwrite(&pos, SSABYTES, 1, esa_file) != 1)
                error("SA write error 1");
            if (fwrite(&esa, SSABYTES, 1, esa_file) != 1)
                error("SA write error 2");
        }
    }

    inline void print_bwt()
    {   
        for (size_t i = 0; i < length; ++i)
        {
            if (fputc(head, bwt_file) == EOF)
                error("BWT write error 1");
        }
    }

    inline void update_bwt(uint8_t next_char, size_t length_)
    {
        if (head != next_char)
        {
            print_sa();
            print_bwt();

            head = next_char;

            length = 0;
        }
        length += length_;

    }

};

#endif /* end of include guard: _LCP_PFP_HH */
