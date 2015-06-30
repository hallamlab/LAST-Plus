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
#include <string>
#include <vector>

#include <queue>

#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include "semaphores.hh"
#include "externalsort.hh"
#include "utilities.hh"

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args->verbosity > 0 ) std::cerr << "lastal: " << x << '\n'

#define INPUT_SIZE 25000

using namespace cbrc;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

namespace {

  LastalArguments *args;
  LambdaCalculator *lambdaCalculator;
  int minScoreGapless;
  int isCaseSensitiveSeeds = -1;  // initialize it to an "error" value
  unsigned numOfIndexes = 1;  // assume this value, if unspecified
  sequenceFormat::Enum referenceFormat = sequenceFormat::fasta;
  SubsetSuffixArray *suffixArrays;
  MultiSequence *text;
  indexT minSeedLimit;

  Alphabet *alph;
  Alphabet *queryAlph;  // for translated alignment
  GeneticCode *geneticCode;

  GeneralizedAffineGapCosts gapCosts;
  ScoreMatrix *scoreMatrix;
  OneQualityScoreMatrix *oneQualityScoreMatrix;
  OneQualityScoreMatrix *oneQualityScoreMatrixMasked;
  OneQualityExpMatrix *oneQualityExpMatrix;
  QualityPssmMaker *qualityPssmMaker;
  TwoQualityScoreMatrix *twoQualityScoreMatrix;
  TwoQualityScoreMatrix *twoQualityScoreMatrixMasked;
}

namespace Phase{ 

  enum Enum{ gapless, gapped, final }; 
}

struct threadData{

  GappedXdropAligner *gappedXdropAligner;
  Centroid *centroid;

  MultiSequence *query;  
  
  std::queue<MultiSequence*> *queryQueue; 
  std::vector< std::vector<countT> > *matchCounts;  // used if outputType == 0

  std::vector<std::string> *outputVector;
  std::queue< std::vector<std::string>* > *outputVectorQueue;

  int identifier;
  int round;

  SEM_T readSema;
  SEM_T writeSema;

  // Find query matches to the suffix array, and do gapless extensions
  void alignGapless( SegmentPairPot& gaplessAlns, char strand );
  // Do gapped extensions of the gapless alignments
  void alignGapped( AlignmentPot& gappedAlns, SegmentPairPot& gaplessAlns, Phase::Enum phase );
  // Print the gapped alignments, after optionally calculating match
  // probabilities and re-aligning using the gamma-centroid algorithm
  void alignFinish( const AlignmentPot& gappedAlns, char strand );
  void makeQualityPssm( bool isApplyMasking );
  // Scan one batch of query sequences against one database volume
  void scan( char strand );
  // Scan one batch of query sequences against one database volume,
  // after optionally translating the query
  void translateAndScan( char strand );
  void reverseComplementPssm();
  void reverseComplementQuery();
  // Scan one batch of query sequences against all database volumes
  void scanAllVolumes();
  void prepareThreadData(int identifier );
  void countMatches( char strand );
  // Write match counts for each query sequence
  void writeCounts(std::ostream& out);
};

struct Dispatcher{

  const uchar* a;  // the reference sequence
  const uchar* b;  // the query sequence
  const uchar* i;  // the reference quality data
  const uchar* j;  // the query quality data
  const ScoreMatrixRow* p;  // the query PSSM
  const ScoreMatrixRow* m;  // the score matrix
  const TwoQualityScoreMatrix& t;
  int d;  // the maximum score drop
  int z;
  Alphabet *aa;


  Dispatcher( Phase::Enum e, MultiSequence *text, MultiSequence *query,
      ScoreMatrix *scoreMatrix, TwoQualityScoreMatrix *twoQualityScoreMatrix, 
      TwoQualityScoreMatrix *twoQualityScoreMatrixMasked, 
      sequenceFormat::Enum referenceFormat, Alphabet *alph) :

    a  ( text->seqReader() ),
    b  ( query->seqReader() ),
    i  ( text->qualityReader() ),
    j  ( query->qualityReader() ),
    p  ( query->pssmReader() ),
    m  ( (e < args->maskLowercase) ?
        scoreMatrix->caseSensitive : scoreMatrix->caseInsensitive ),
    t  ( (e < args->maskLowercase) ?
        *twoQualityScoreMatrixMasked : *twoQualityScoreMatrix ),
    d  ( (e == Phase::gapless) ? args->maxDropGapless :
        (e == Phase::gapped ) ? args->maxDropGapped : args->maxDropFinal ),
    z  ( (args->inputFormat == sequenceFormat::fasta) ? 0 :
        (referenceFormat  == sequenceFormat::fasta) ? 1 : 2 ),
    aa ( aa = alph ){}

  // Shrink the SegmentPair to its longest run of identical matches.
  // This trims off possibly unreliable parts of the gapless alignment.
  // It may not be the best strategy for protein alignment with subset
  // seeds: there could be few or no identical matches...
  void shrinkToLongestIdenticalRun(SegmentPair &sp) {
    sp.maxIdenticalRun(a, b, aa->canonical);
    sp.score = gaplessScore(sp.beg1(), sp.end1(), sp.beg2());
  }

  int forwardGaplessScore( indexT x, indexT y ) const{
    if( z==0 ) return forwardGaplessXdropScore( a+x, b+y, m, d );
    if( z==1 ) return forwardGaplessPssmXdropScore( a+x, p+y, d );
    return forwardGaplessTwoQualityXdropScore( a+x, i+x, b+y, j+y, t, d );
  }

  int reverseGaplessScore( indexT x, indexT y ) const{
    if( z==0 ) return reverseGaplessXdropScore( a+x, b+y, m, d );
    if( z==1 ) return reverseGaplessPssmXdropScore( a+x, p+y, d );
    return reverseGaplessTwoQualityXdropScore( a+x, i+x, b+y, j+y, t, d );
  }

  indexT forwardGaplessEnd( indexT x, indexT y, int s ) const{
    if( z==0 ) return forwardGaplessXdropEnd( a+x, b+y, m, s ) - a;
    if( z==1 ) return forwardGaplessPssmXdropEnd( a+x, p+y, s ) - a;
    return forwardGaplessTwoQualityXdropEnd( a+x, i+x, b+y, j+y, t, s ) - a;
  }

  indexT reverseGaplessEnd( indexT x, indexT y, int s ) const{
    if( z==0 ) return reverseGaplessXdropEnd( a+x, b+y, m, s ) - a;
    if( z==1 ) return reverseGaplessPssmXdropEnd( a+x, p+y, s ) - a;
    return reverseGaplessTwoQualityXdropEnd( a+x, i+x, b+y, j+y, t, s ) - a;
  }

  bool isOptimalGapless( indexT x, indexT e, indexT y ) const{
    if( z==0 ) return isOptimalGaplessXdrop( a+x, a+e, b+y, m, d );
    if( z==1 ) return isOptimalGaplessPssmXdrop( a+x, a+e, p+y, d );
    return isOptimalGaplessTwoQualityXdrop( a+x, a+e, i+x, b+y, j+y, t, d );
  }

  int gaplessScore( indexT x, indexT e, indexT y ) const{
    if( z==0 ) return gaplessAlignmentScore( a+x, a+e, b+y, m );
    if( z==1 ) return gaplessPssmAlignmentScore( a+x, a+e, p+y );
    return gaplessTwoQualityAlignmentScore( a+x, a+e, i+x, b+y, j+y, t );
  }
};

void readOuterPrj( const std::string& fileName, unsigned& volumes,
    countT& refSequences, countT& refLetters );
// Read a per-volume .prj file, with info about a database volume
void readInnerPrj( const std::string& fileName, indexT& seqCount, indexT& seqLen );
// Read one database volume
void readVolume( unsigned volumeNumber );
void readIndex( const std::string& baseName, indexT seqCount );
std::istream& appendFromFasta( std::istream& in, MultiSequence *query );

void *writerFunction(void *arguments);
void readerFunction( std::istream& in );
void *threadFunction(void *__threadData);
void writeHeader( std::ostream& out );
void initializeThreads();
void initializeSemaphores();
void initializeEvalueCalulator(const std::string dbPrjFile, ScoreMatrix *scoreMatrix,
    std::string dbfilePrj);

void createStructures(std::string &matrixFile);
  // Set up a scoring matrix, based on the user options
  void makeScoreMatrix( const std::string& matrixFile) ;
  void makeQualityScorers();
  // Calculate statistical parameters for the alignment scoring scheme
  // Meaningless for PSSMs, unless they have the same scale as the score matrix
  void calculateScoreStatistics();

void lastal( int argc, char** argv );

#endif
