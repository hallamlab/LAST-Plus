#include "heapsort.hh"
#include <iostream>

using namespace std;

void print_heap(vector<pair<int, Line *> >& A) {
  unsigned int j; 
  for (j = 0; j<A.size(); j++) {
    cout << "(" << A[j].first << ", " << A[j].second << ") ";
  }
  cout << "\n";
}

void heapify(vector<pair<int, Line *> >& A, int i, int S) {
  while (true) {
    int l = (2*i) + 1;
    int r = l + 1;
    int max = i;


    if (l < S && A[l].second->orfid < A[max].second->orfid) {
      max = l;
    } else if (l < S && A[l].second->orfid == A[max].second->orfid) {
      if (l < S && A[l].second->evalue < A[max].second->evalue) {
        max = l;
      }
    }

    if (r < S && A[r].second->orfid < A[max].second->orfid) {
      max = r;
    } else if (r < S && A[r].second->orfid == A[max].second->orfid) {
      if (r < S && A[r].second->evalue < A[max].second->evalue) {
        max = r;
      }
    }

    if (max != i && i < S) {
      pair<int, Line *> temp = A[i];
      A[i] = A[max];
      A[max] = temp;
    } else {
      break;
    }
    i = max;
  }
}

void build_heap(int S, vector<pair<int, Line *> >& A) {
  int i = floor(S/2);
  while (i >= 0) {
    heapify(A, i, S);
    i--;
  }
}
