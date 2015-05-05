#ifndef __UTILS__
#define __UTILS__

#include <iostream>
#include <string>
#include <fstream>

#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "io.hh"

struct SequenceStatistics{
  double letterCount;
  double sequenceCount;
  std::string alphabet;
  std::vector<double> letterProbs;

  SequenceStatistics()
      : letterCount(0), sequenceCount(0), alphabet(""), letterProbs(0) {}
};

void makeAlphabet( cbrc::Alphabet & alph, bool isProtein=false );

void fastaFileSequenceStats( std::string fastaFile, SequenceStatistics *stats );

bool __isDubiousDna( const cbrc::Alphabet& alph, const cbrc::MultiSequence& multi );




#endif  // __UTILS__ guard end
