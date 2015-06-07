#ifndef  _EXTERNAL_SORT
#define  _EXTERNAL_SORT

#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "heapsort.hh"
#include "linereader.hh"
#include "utilities.hh"
/*
All functions required for the initial sequence sorting and blocking
*/

using namespace std;

// Function object for sorting sequence id/length pairs.
// On ties, use increasing order of sequence id numbers.
struct sort_pred {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        if (left.second == right.second) {
			return right.first > left.first;	
		} else {
			return left.second > right.second;
		}
	}
};


int disk_sort_file(string outputdir, string tobe_sorted_file_name, string sorted_file_name,\
     int chunk_size, string(*key_extractor)(const string &) ) ;

int sort_and_create_blocks(string outputdir, string faa_file, float block_mb, map<int, int> &block_active_seqs, map<string, int> &seq_lengths); 

int merge_sorted_files_create_blocks(vector<string>& filenames, float block_mb, string outputdir, string sorted_filename); 

void write_sorted_sequences(vector<Line *>& lines, string filename); 

void remove_file(string filename); 
#endif // _EXTERNAL_SORT
