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

int disk_sort_file(string outputdir, string tobe_sorted_file_name, string sorted_file_name, countT chunk_size, string(*key_extractor)(const string &), const std::vector<std::string> &mergelist);
int merge_sorted_files(const vector<string> &filenames, string sorted_file_name);
void remove_file(string filename); 
string generate_directory_name();

#endif // _EXTERNAL_SORT
