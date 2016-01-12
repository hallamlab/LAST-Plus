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

std::string random_str(const int len);

std::string extract_sequence_name(const std::string &name);

void topHits(std::string filename, int maxHits);
void topHitsVector(std::vector<std::string> &outputVector, int maxHits);

string orf_extractor_from_blast(const string &line);

double evalue_extractor_from_blast(const string &line);

double bit_score_extractor_from_blast(const string &line);

// Niels' debug functions
void testFunction();
void printOutSequence(const unsigned char* field, const unsigned char decode[]);

#endif //_UTILITIES
