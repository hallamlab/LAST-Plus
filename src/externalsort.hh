#ifndef  _EXTERNAL_SORT
#define  _EXTERNAL_SORT

#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>

#include "heapsort.hh"
#include "linereader.hh"
#include "utilities.hh"
#include "lastal.hh"
#include "tempfiles.hh"

#define STR_BUFFER_SIZE 1000
#define MERGE_SIZE 100
/*
   All functions required for the initial sequence sorting and blocking
   */

using namespace std;

typedef unsigned long long countT;
typedef Line *LINE;

int disk_sort_file(const string &outputdir, 
    const string &tobe_sorted_file_name, 
    const string &sorted_file_name, 
    countT chunk_size, 
    string(*key_extractor)(const string &), 
    const std::vector<std::string> &mergelist);

std::vector<std::string> merge_some_files(const std::vector<std::string> &mergelist, 
    std::vector<TEMPFILES*> &directories,
    const string &tmpdir);

int merge_sorted_files(const vector<string> &filenames, 
    const string &sorted_file_name,
    const string &tmpdir);

void remove_file(const string &filename); 

string generate_directory_name(const string &tmpdir);

void write_sorted_sequences(vector<Line *> &lines, 
    const string &filename);

#endif // _EXTERNAL_SORT
