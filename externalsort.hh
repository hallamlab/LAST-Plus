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

typedef unsigned long long countT;
typedef Line *LINE;

int disk_sort_file(string outputdir, string tobe_sorted_file_name, string sorted_file_name,
     countT chunk_size, string(*key_extractor)(const string &) ) ;

int merge_sorted_files_create_blocks(vector<string> &filenames, string sorted_file_name);

void write_sorted_sequences(vector<Line *>& lines, string filename); 

void remove_file(string filename); 
#endif // _EXTERNAL_SORT
