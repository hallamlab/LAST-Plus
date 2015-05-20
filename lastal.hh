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
#include <set>


#include <cstdlib>
#include <pthread.h>
#include "semaphores.hh"
#include "outputStruct.hh"
#include "SubsetSuffixArrayUser.hh"


#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastal: " << x << '\n'


using namespace cbrc;


typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

namespace {

  LastalArguments args;
  LambdaCalculator lambdaCalculator;
  const unsigned maxNumOfIndexes = 16;
  int minScoreGapless;
  int isCaseSensitiveSeeds = -1;  // initialize it to an "error" value
  unsigned numOfIndexes = 1;  // assume this value, if unspecified

  SubsetSuffixArray suffixArrays[16];
  MultiSequence text;

}

namespace Phase{ 

  enum Enum{ gapless, gapped, final }; 
}

struct threadData{

  Alphabet alph;
  Alphabet queryAlph;  // for translated alignment
  GeneticCode geneticCode;
  //SubsetSuffixArray suffixArrays[16];
  SubsetSuffixArrayUser subsetUser;
  GappedXdropAligner gappedXdropAligner;
  Centroid *centroid;
  MultiSequence query;  // sequence that hasn't been indexed by lastdb
  MultiSequence text;  // sequence that has been indexed by lastdb
  std::vector< std::vector<countT> > matchCounts;  // used if outputType == 0
  GeneralizedAffineGapCosts gapCosts;
  ScoreMatrix scoreMatrix;
  sequenceFormat::Enum referenceFormat;  // defaults to 0
  OneQualityScoreMatrix oneQualityScoreMatrix;
  OneQualityScoreMatrix oneQualityScoreMatrixMasked;
  OneQualityExpMatrix oneQualityExpMatrix;
  QualityPssmMaker qualityPssmMaker;
  TwoQualityScoreMatrix twoQualityScoreMatrix;
  TwoQualityScoreMatrix twoQualityScoreMatrixMasked;
  outputStruct *output;
  int identifier;


  void alignGapless( SegmentPairPot& gaplessAlns, char strand );
  void alignGapped( AlignmentPot& gappedAlns, SegmentPairPot& gaplessAlns, Phase::Enum phase );
  void alignFinish( const AlignmentPot& gappedAlns, char strand );
  void makeQualityPssm( bool isApplyMasking );
  void scan( char strand );
  void translateAndScan( char strand );
  void reverseComplementPssm();
  void reverseComplementQuery();
  void scanAllVolumes( unsigned volumes );
  void prepareThreadData(std::string matrixFile, int identifier );
  void readIndex( const std::string& baseName, indexT seqCount );
  void readVolume( unsigned volumeNumber );
  void countMatches( char strand );
  void writeCounts();
  std::istream& appendFromFasta( std::istream& in );
  void callReinit();

  void makeScoreMatrix( const std::string& matrixFile) ;
  void makeQualityScorers();
  void calculateScoreStatistics();
  void readOuterPrj( const std::string& fileName, unsigned& volumes, indexT& minSeedLimit,
      countT& refSequences, countT& refLetters );
  void readInnerPrj( const std::string& fileName, indexT& seqCount, indexT& seqLen );
  //void initializeEvalueCalulator(const std::string dbPrjFile, std::string dbfilePrj);
  //void writeHeader( countT refSequences, countT refLetters, std::ostream& out );
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
	MultiSequenceUser user;


  Dispatcher( Phase::Enum e, MultiSequence &text, MultiSequence &query,
      ScoreMatrix &scoreMatrix, TwoQualityScoreMatrix &twoQualityScoreMatrix, 
      TwoQualityScoreMatrix &twoQualityScoreMatrixMasked, 
      sequenceFormat::Enum referenceFormat, Alphabet &alph) :

    a  ( user.seqReader(text) ),
    b  ( query.seqReader() ),
    i  ( user.qualityReader(text) ),
    j  ( query.qualityReader() ),
    p  ( query.pssmReader() ),
    m  ( (e < args.maskLowercase) ?
        scoreMatrix.caseSensitive : scoreMatrix.caseInsensitive ),
    t  ( (e < args.maskLowercase) ?
        twoQualityScoreMatrixMasked : twoQualityScoreMatrix ),
    d  ( (e == Phase::gapless) ? args.maxDropGapless :
        (e == Phase::gapped ) ? args.maxDropGapped : args.maxDropFinal ),
    z  ( (args.inputFormat == sequenceFormat::fasta) ? 0 :
        (referenceFormat  == sequenceFormat::fasta) ? 1 : 2 ),
    aa ( aa = &alph ){}

  void shrinkToLongestIdenticalRun( SegmentPair& sp);

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

//void writeHeader( countT refSequences, countT refLetters, std::ostream& out );

void writerFunction( std::ostream& out );
void readerFunction( std::istream& in );
void finishAlignment( std::ostream& out );
void* threadFunction( void *args ); 
void initializeThreads();
void initializeSemaphores();
void initializeEvalueCalulator(const std::string dbPrjFile, ScoreMatrix &scoreMatrix, std::string dbfilePrj);
void lastal( int argc, char** argv );

#endif
