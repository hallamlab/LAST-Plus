#ifndef ___UTILITIES___
#define ___UTILITIES___
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <regex.h>

using namespace std;


void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');

char * split_n_pick(const string  &strn,  char *buf, char d, unsigned int n);

void split_seq_name(const std::string  &strn, std::vector<char *> &v, char *buf);

std::string get_orf_name(std::string & strn, std::vector<char *> &v, char *buf);

bool matchString(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

std::string extract_sequence_name(const std::string &name);

string to_string(unsigned long i);

string ShortenORFId(const string &s, regex_t *r);

string ShortenORFId(const string &s);

#define MAX_ERROR_MSG 0x1000
int compile_regex(regex_t * r, const char * regex_text);

string getpattern(regex_t *r , const char *to_match, unsigned int no );

string getECNo(const char *str, unsigned int d);

string getCOGID(const char *str);
string getKEGGID(const char *str);
string getSEEDID(const char *str);

void topHits(std::string filename, int maxHits);

string getTaxonomyFromProduct(const char *str);

string orf_extractor_from_blast(const string &line);

double evalue_extractor_from_blast(const string &line);

unsigned int hashIntoBucket(const char *str, unsigned int index) ;

unsigned long long hashStringIntoBucket(const char *str, unsigned long long index) ;

string function_extractor_from_list(const string & line);

string to_upper(const string &str);

#endif //_UTILITIES

