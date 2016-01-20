#include "utils.hh"

using namespace cbrc;

// Does the first sequence look like it isn't really DNA?
bool __isDubiousDna( const Alphabet& alph, const MultiSequence& multi ){
  const uchar* seq = multi.seqReader() + multi.seqBeg(0);
  unsigned dnaCount = 0;

  for( MultiSequence::indexT i = 0; i < multi.seqLen(0); i++ ){  // look at the first 100 letters
    uchar c = alph.canonical[ seq[i] ];
    if( c == alph.size ) return false;  // we hit the end of the sequence early

    if( c > alph.size || c == alph.encode[ (uchar)'N' ] ) ++dnaCount;
  }

  if( dnaCount < static_cast<unsigned>(0.90*multi.seqLen(0)) ) return true;  // more than 10% unexpected letters
  else return false;
}


// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet( Alphabet& alph, bool isProtein )
{
  if( isProtein )         
    alph.fromString( alph.protein );
  else                              
    alph.fromString( alph.dna );
}

bool  isFastaFileProtein(std::string fastaFile) {
  using namespace cbrc;
  Alphabet alph;
  MultiSequence multi;
  makeAlphabet( alph );
  multi.initForAppending(1);
  alph.tr( multi.seqWriter(), multi.seqWriter() + multi.unfinishedSize() );

  std::vector<unsigned long long> letterCounts( alph.size );
  std::vector<int> letterTotals( alph.size );


  std::ifstream inFileStream;
  std::istream& in = openIn( fastaFile, inFileStream );
  multi.reinitForAppending();

  unsigned int sequenceCount =0;
  unsigned int lettersTotal =0;
  while( sequenceCount < 5 && in.good() ){
    try{
        multi.appendFromFasta(in, 10000000) ;
         alph.tr( multi.seqWriter() + multi.seqBeg(0), multi.seqWriter() + multi.seqEnd(0) );
         alph.count( multi.seqReader() + multi.seqBeg(0), multi.seqReader() + multi.seqEnd(0), &letterCounts[0] );
         lettersTotal += multi.seqLen(0);
         for( unsigned c = 0; c < alph.size; ++c ) 
            letterTotals[c] += letterCounts[c];
	     letterCounts.assign( alph.size, 0 );
	     multi.reinitForAppending();
         sequenceCount++;
    } catch (const std::exception &ex){
        std::cerr << ex.what() << std::endl;
        std::cerr << "Encountered a malformed sequence. Ignoring sequence and continuing" << std::endl;
        multi.printOffensiveName();
	    multi.reinitForAppending();
        //std::cerr << "Caught the exception" << std::endl;
    }
  }

  unsigned int lettersSum = 0;
  for( unsigned c = 0; c < alph.size; ++c ) {
      lettersSum += letterTotals[c];
  }
  if(static_cast<float>(lettersSum)/static_cast<float>(lettersTotal) < 0.9) return true;

  return false;
}

void fastaFileSequenceStats( std::string fastaFile, 
                            SequenceStatistics *stats, 
                            cbrc::sequenceFormat::Enum format){

  using namespace cbrc;

  Alphabet alph;
  MultiSequence multi;
  
  bool isProtein = false;
    if( !isFastq(format) ){
        isProtein = isFastaFileProtein(fastaFile);
    }

  makeAlphabet( alph, isProtein );
  multi.initForAppending(1);
  alph.tr( multi.seqWriter(), multi.seqWriter() + multi.unfinishedSize() );

  int sequenceCount = 0;
  std::vector<unsigned long long> letterCounts( alph.size );
  std::vector<int> letterTotals( alph.size );


  std::ifstream inFileStream;
  std::istream& in = openIn( fastaFile, inFileStream );
  multi.reinitForAppending();

  sequenceCount = 0;
    
    if(!isFastq(format)){
  while(  in.good() ){
    try{
        multi.appendFromFasta(in, 1000000000);
         alph.tr( multi.seqWriter() + multi.seqBeg(0), multi.seqWriter() + multi.seqEnd(0) );
         alph.count( multi.seqReader() + multi.seqBeg(0), multi.seqReader() + multi.seqEnd(0), &letterCounts[0] );
         for( unsigned c = 0; c < alph.size; ++c ) 
            letterTotals[c] += letterCounts[c];
         
	     letterCounts.assign( alph.size, 0 );
	     multi.reinitForAppending();
         sequenceCount++;
    } catch (const std::exception &ex){
        //std::cerr << ex.what() << std::endl;
        //std::cerr << "Encountered a malformed sequence. Ignoring sequence and continuing" << std::endl;
        //multi.printOffensiveName();
	    multi.reinitForAppending();
        //std::cerr << "Caught the exception" << std::endl;
    }
  }

    }else{

  while(  in.good() ) {

    try{

      multi.appendFromFastq(in, 1000000000);
         alph.tr( multi.seqWriter() + multi.seqBeg(0), multi.seqWriter() + multi.seqEnd(0) );
         alph.count( multi.seqReader() + multi.seqBeg(0), multi.seqReader() + multi.seqEnd(0), &letterCounts[0] );
         for( unsigned c = 0; c < alph.size; ++c ) 
            letterTotals[c] += letterCounts[c];
         
	     letterCounts.assign( alph.size, 0 );
	     multi.reinitForAppending();
         sequenceCount++;

    } catch (const std::exception &ex){
        //std::cerr << ex.what() << std::endl;
        //std::cerr << "Encountered a malformed sequence. Ignoring sequence and continuing" << std::endl;
        //multi.printOffensiveName();
	    multi.reinitForAppending();
        //std::cerr << "Caught the exception" << std::endl;
    }
  }

    }


  stats->sequenceCount = sequenceCount;
  stats->letterCount =0;
  for( unsigned c = 0; c < alph.size; ++c ) {
       stats->letterCount += letterTotals[c];
  }
  for( unsigned c = 0; c < alph.size; ++c ) {
     stats->letterProbs.push_back(static_cast<double>(letterTotals[c])/static_cast<double>(stats->letterCount));
  }

  stats->alphabet = alph.letters; 
}

