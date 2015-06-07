#ifndef  __HEAPSORT
#define  __HEAPSORT

#include <vector>
#include <math.h>
#include <utility>
#include "linereader.hh"
using namespace std;

/*
Heap-sort a vector of key-value pairs in descending order of their values.
*/

void print_heap(vector<pair<int, Line *> >& A);
void heapify(vector< pair<int, Line *> >& A, int i, int S); 
void build_heap(int S, vector<pair<int, Line *> >& A);

#endif   // __HEAPSORT
