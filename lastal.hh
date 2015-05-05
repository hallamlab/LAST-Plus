#ifndef __LASTAL_HH
#define __LASTAL_HH

#include "LastalArguments.hh"
#include "QualityPssmMaker.hh"
#include "OneQualityScoreMatrix.hh"
#include "TwoQualityScoreMatrix.hh"
#include "qualityScoreUtil.hh"
#include "LambdaCalculator.hh"
#include "GeneticCode.hh"
#include "SubsetSuffixArray.hh"
#include "Centroid.hh"
#include "GappedXdropAligner.hh"
#include "AlignmentPot.hh"
#include "Alignment.hh"
#include "SegmentPairPot.hh"
#include "SegmentPair.hh"
#include "ScoreMatrix.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "DiagonalTable.hh"
#include "GeneralizedAffineGapCosts.hh"
#include "gaplessXdrop.hh"
#include "gaplessPssmXdrop.hh"
#include "gaplessTwoQualityXdrop.hh"
#include "io.hh"
#include "stringify.hh"
#include "lastex.hh"
#include "LastexArguments.hh"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE


#include <pthread.h>
#include <semaphore.h>
#include <cstring> //memcpy

using namespace cbrc;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

struct threadData{
  //LastalArguments args;
  Alphabet alph;
  Alphabet queryAlph;  // for translated alignment
  GeneticCode geneticCode;
  //const unsigned maxNumOfIndexes = 16;
  SubsetSuffixArray suffixArrays[16];
  //ScoreMatrix scoreMatrix;
  //GeneralizedAffineGapCosts gapCosts;
  GappedXdropAligner gappedXdropAligner;
  Centroid centroid( GappedXdropAligner );
  //Centroid centroid;
  //LambdaCalculator lambdaCalculator;
  MultiSequence query;  // sequence that hasn't been indexed by lastdb
  MultiSequence text;  // sequence that has been indexed by lastdb
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  OneQualityScoreMatrix oneQualityScoreMatrix;
  OneQualityScoreMatrix oneQualityScoreMatrixMasked;
  OneQualityExpMatrix oneQualityExpMatrix;
  QualityPssmMaker qualityPssmMaker;
  //sequenceFormat::Enum referenceFormat;  // defaults to 0
  TwoQualityScoreMatrix twoQualityScoreMatrix;
  TwoQualityScoreMatrix twoQualityScoreMatrixMasked;
  //int minScoreGapless;
  //int isCaseSensitiveSeeds = -1;  // initialize it to an "error" value
  //unsigned numOfIndexes = 1;  // assume this value, if unspecified
};

void prepareThreadData();

#endif
