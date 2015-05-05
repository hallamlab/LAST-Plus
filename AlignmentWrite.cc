// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"

#include <iomanip>
#include <algorithm>
#include <cassert>
#include <iterator>  // ostream_iterator

#include <cmath>
#include "lastex.hh"
#include "gumbel_params/sls_pvalues.hpp"


// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

static int numDigits( size_t x ){
  int n = 0;
  do{
    ++n;
    x /= 10;
  }while(x);
  return n;
}

// Printing with either C++ streams or sprintf can be noticeably slow.
// So the next 3 functions are used instead.

static char* sprintLeft( char* dest, const char* src, int width ){
  char* end = dest + width;
  while( *src ) *dest++ = *src++;
  while( dest < end ) *dest++ = ' ';
  *dest++ = ' ';
  return dest;
}

static char* sprintSize( char* dest, size_t size, int width ){
  char* end = dest + width;
  char* beg = end;

  do{
    --beg;
    *beg = '0' + size % 10;
    size /= 10;
  }while( size );

  while( dest < beg ) *dest++ = ' ';

  *end++ = ' ';
  return end;
}

static char* sprintChar( char* dest, char c ){
  *dest++ = c;
  *dest++ = ' ';
  return dest;
}

// write x - y as a signed integer
static void writeSignedDifference( size_t x, size_t y, std::ostream& os ){
  if( x >= y )  os << x - y;
  else          os << '-' << y - x;
}

void Alignment::write( const MultiSequence& seq1, const MultiSequence& seq2,
		       char strand, bool isTranslated, const Alphabet& alph,
		       int format, std::ostream& os, LastalArguments &args, 
           const AlignmentExtras& extras) const{
  assert( !blocks.empty() );

  if( format == 0 ) 
       writeTab( seq1, seq2, strand, isTranslated, os, extras );
  else if( format == 2 )  {
       writeBlastOutput( seq1, seq2, strand, isTranslated, alph,os, extras, args );
  }
  else 
       writeMaf( seq1, seq2, strand, isTranslated, alph, os, extras );

}

//!!
void Alignment::writeBlastOutput( const MultiSequence& seq1, const MultiSequence& seq2,
              char strand, bool isTranslated, const Alphabet& alph, std::ostream& os,
			  const AlignmentExtras& extras, LastalArguments &args ) const{

  double fullScore = extras.fullScore;

  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);
  size_t b1 = alnBeg1 - seqStart1;
  size_t b2 = alnBeg2 - seqStart2;
  size_t r1 = alnEnd1 - alnBeg1;
  size_t r2 = alnEnd2 - alnBeg2;
  size_t s1 = seq1.seqLen(w1);
  size_t s2 = seq2.seqLen(w2);

  const int nw = std::max( n1.size(), n2.size() );
  const int bw = std::max( numDigits(b1), numDigits(b2) );
  const int rw = std::max( numDigits(r1), numDigits(r2) );
  const int sw = std::max( numDigits(s1), numDigits(s2) );

  char tab = '\t';

  double evalue = 0;
  size_t identities = 0;
  size_t mismatches = 0;
  size_t gaps = 0;
  size_t alignLength = alnEnd1 - alnBeg1;

  size_t headLen = 2 + nw + 1 + bw + 1 + rw + 3 + sw + 1;
  size_t lineLen = headLen + numColumns( frameSize2 ) + 1;
  std::vector<char> lineVector( lineLen );
  char* line = &lineVector[0];
  line[ lineLen - 1 ] = '\n';
  char* dest;

  dest = sprintChar( line, 's' );
  dest = sprintLeft( dest, n1.c_str(), nw );
  dest = sprintSize( dest, b1, bw );
  dest = sprintSize( dest, r1, rw );
  dest = sprintChar( dest, '+' );
  dest = sprintSize( dest, s1, sw );
  
  writeTopSeq( seq1.seqReader(), alph, frameSize2, dest );
  std::string seq1String = getSequence(line, lineLen);
  
  gaps += countGaps(seq1String);
 
  dest = sprintChar( line, 's' );
  dest = sprintLeft( dest, n2.c_str(), nw );
  dest = sprintSize( dest, b2, bw );
  dest = sprintSize( dest, r2, rw );
  dest = sprintChar( dest, strand );
  dest = sprintSize( dest, s2, sw );

  writeBotSeq( seq2.seqReader(), alph, frameSize2, dest );
  std::string seq2String = getSequence(line, lineLen);

  gaps += countGaps(seq2String);
  size_t identityCount = countIdentities(seq1String, seq2String);
  identities = 100*identityCount/seq1String.length();  
  mismatches = alignLength - gaps - identityCount;

  evalue = evalueForSequences(score,s1, s2);

  double lambda = getLambda();
  double k = getK();
  double area = getArea();

  double bitscore = (lambda*score-log(k))/log(2);
  double evalue2 = area*pow(2,-bitscore);

  if(bitscore >= args.scoreCutoff && evalue2 <= args.evalueCutoff){
  
    os << seq2.seqName(w2) << tab
       << seq1.seqName(w1) << tab
       << identities << tab
       << alignLength << tab
       << mismatches << tab
       << gaps << tab
       << (alnBeg2 - seqStart2) << tab
       << (alnEnd2 -seqStart2) << tab
       << (alnBeg1 - seqStart1) << tab
       << (alnEnd1 -seqStart1) << tab
       << evalue2 << tab
       << bitscore;
       os << '\n';
  }

}

size_t Alignment::countIdentities(std::string& seq1String, std::string& seq2String) const{
    assert(seq1String.length() == seq2String.length());

    size_t count = 0;
    for (size_t i=0;i<seq1String.length();i++){
        if(seq1String[i] == seq2String[i]) count++;
    }
    return count;
}

// the crudest possible thing... retrieve the last element in this
// space delimited string
std::string Alignment::getSequence(char*& line, size_t& length) const {
    std::vector<std::string> columnEntries;
    std::string lineString(line, length);

    std::string delimiter = " ";
    size_t pos = 0;
    std::string token;
    while ((pos = lineString.find(delimiter)) != std::string::npos) {
        token = lineString.substr(0, pos);
        columnEntries.push_back(token);
        lineString.erase(0, pos + delimiter.length());
    }
    lineString.erase(std::remove(lineString.begin(), lineString.end(), '\n'), lineString.end());
    return lineString;
}

size_t Alignment::countGaps(std::string& sequence) const {
    return std::count(sequence.begin(), sequence.end(), '-'); 
}

void Alignment::writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated, std::ostream& os,
			  const AlignmentExtras& extras ) const{
  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  os << score << '\t';

  os << seq1.seqName(w1) << '\t'
     << alnBeg1 - seqStart1 << '\t'
     << alnEnd1 - alnBeg1 << '\t'
     << '+' << '\t'
     << seq1.seqLen(w1) << '\t';

  os << seq2.seqName(w2) << '\t'
     << alnBeg2 - seqStart2 << '\t'
     << alnEnd2 - alnBeg2 << '\t'
     << strand << '\t'
     << seq2.seqLen(w2) << '\t';

  std::cout << (alnBeg1 - seqStart1) << "  " << (alnEnd1 -seqStart1) << " --  " << (alnBeg2 - seqStart2) << "  " << (alnEnd2 -seqStart2) << std::endl;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;
      if( j->size ) os << ',';
      size_t gapBeg1 = j->end1();
      size_t gapEnd1 = i->beg1();
      writeSignedDifference( gapEnd1, gapBeg1, os );  // allow -1 frameshift
      os << ':';
      size_t gapBeg2 = aaToDna( j->end2(), frameSize2 );
      size_t gapEnd2 = aaToDna( i->beg2(), frameSize2 );
      writeSignedDifference( gapEnd2, gapBeg2, os );  // allow -1 frameshift
      if( i->size ) os << ',';
    }
    if( i->size ) os << i->size;
  }

  double fullScore = extras.fullScore;
  if( fullScore > 0 ) os << "\tfullScore=" << fullScore;

  os << '\n';
}

void Alignment::writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated, const Alphabet& alph,
			  std::ostream& os,
			  const AlignmentExtras& extras ) const{
  double fullScore = extras.fullScore;
  const std::vector<uchar>& columnAmbiguityCodes = extras.columnAmbiguityCodes;
  const std::vector<double>& expectedCounts = extras.expectedCounts;

  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);
  size_t b1 = alnBeg1 - seqStart1;
  size_t b2 = alnBeg2 - seqStart2;
  size_t r1 = alnEnd1 - alnBeg1;
  size_t r2 = alnEnd2 - alnBeg2;
  size_t s1 = seq1.seqLen(w1);
  size_t s2 = seq2.seqLen(w2);

  const int nw = std::max( n1.size(), n2.size() );
  const int bw = std::max( numDigits(b1), numDigits(b2) );
  const int rw = std::max( numDigits(r1), numDigits(r2) );
  const int sw = std::max( numDigits(s1), numDigits(s2) );

  size_t headLen = 2 + nw + 1 + bw + 1 + rw + 3 + sw + 1;
  size_t lineLen = headLen + numColumns( frameSize2 ) + 1;
  std::vector<char> lineVector( lineLen );
  char* line = &lineVector[0];
  line[ lineLen - 1 ] = '\n';
  char* dest;

  os << "a";
  os << " score=" << score;
  if( fullScore > 0 ) os << " fullScore=" << fullScore;
  os << '\n';

  dest = sprintChar( line, 's' );
  dest = sprintLeft( dest, n1.c_str(), nw );
  dest = sprintSize( dest, b1, bw );
  dest = sprintSize( dest, r1, rw );
  dest = sprintChar( dest, '+' );
  dest = sprintSize( dest, s1, sw );
  
  writeTopSeq( seq1.seqReader(), alph, frameSize2, dest );

  os.write( line, lineLen );
  
  if( seq1.qualsPerLetter() > 0 ){
    dest = sprintChar( line, 'q' );
    dest += nw + 1;
    dest = sprintLeft( dest, "", bw + 1 + rw + 3 + sw );
    writeTopQual( seq1.qualityReader(), seq1.qualsPerLetter(), dest );
    os.write( line, lineLen );
  }

  dest = sprintChar( line, 's' );
  dest = sprintLeft( dest, n2.c_str(), nw );
  dest = sprintSize( dest, b2, bw );
  dest = sprintSize( dest, r2, rw );
  dest = sprintChar( dest, strand );
  dest = sprintSize( dest, s2, sw );
  writeBotSeq( seq2.seqReader(), alph, frameSize2, dest );
  os.write( line, lineLen );

  if( seq2.qualsPerLetter() > 0 ){
    dest = sprintChar( line, 'q' );
    dest += nw + 1;
    dest = sprintLeft( dest, "", bw + 1 + rw + 3 + sw );
    writeBotQual( seq2.qualityReader(), seq2.qualsPerLetter(), dest );
    os.write( line, lineLen );
  }

  if( columnAmbiguityCodes.size() > 0 ){
    os << "p "
       << std::setw( nw + bw + rw + sw + 6 ) << "";
    std::copy( columnAmbiguityCodes.begin(), columnAmbiguityCodes.end(),
               std::ostream_iterator<uchar>(os) );
    os << '\n';
  }

  if( expectedCounts.size() > 0 ){
    os << 'c';
    for( unsigned i = 0; i < expectedCounts.size(); ++i )
      os << ' ' << expectedCounts[i];
    os << '\n';
  }

  os << '\n';  // blank line afterwards
}

size_t Alignment::numColumns( size_t frameSize ) const{
  size_t num = 0;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // length of unaligned chunk of top sequence (gaps in bottom sequence):
      num += i->beg1() - j->end1();

      // length of unaligned chunk of bottom sequence (gaps in top sequence):
      size_t gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) ++num;
      num += gap2;
    }

    num += i->size;  // length of aligned chunk
  }

  return num;
}

static char* writeGaps( char* dest, size_t num ){
  char* end = dest + num;
  while( dest < end ){
    *dest++ = '-';
  }
  return dest;
}

char* Alignment::writeTopSeq( const uchar* seq, const Alphabet& alph,
			      size_t frameSize, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // append unaligned chunk of top sequence:
      dest = alph.rtCopy( seq + j->end1(), seq + i->beg1(), dest );

      // append gaps for unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) {
          *dest++ = '-';
      }
      dest = writeGaps( dest, gap2 );
    }
    // append aligned chunk of top sequence:
    dest = alph.rtCopy( seq + i->beg1(), seq + i->end1(), dest );
  }
  return dest;
}

char* Alignment::writeBotSeq( const uchar* seq, const Alphabet& alph,
			      size_t frameSize, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, i->beg1() - j->end1());

      //append unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 == 1 ) *dest++ = '\\';
      if( frameshift2 == 2 ) *dest++ = '/';
      dest = alph.rtCopy( seq + i->beg2() - gap2, seq + i->beg2(), dest );
    }

    // append aligned chunk of bottom sequence:
    dest = alph.rtCopy( seq + i->beg2(), seq + i->end2(), dest );
  }

  return dest;
}

static char* writeQuals( const uchar* qualities, size_t beg, size_t end,
			 size_t qualsPerBase, char* dest ){
  for( size_t i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    *dest++ = *std::max_element( q, q + qualsPerBase );
  }
  return dest;
}

char* Alignment::writeTopQual( const uchar* qualities,
			       size_t qualsPerBase, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // assume we're not doing translated alignment

      // append qualities for unaligned chunk of top sequence:
      dest = writeQuals( qualities, j->end1(), i->beg1(), qualsPerBase, dest );

      // append gaps for unaligned chunk of bottom sequence:
      dest = writeGaps( dest, i->beg2() - j->end2());
    }

    // append qualities for aligned chunk of top sequence:
    dest = writeQuals( qualities, i->beg1(), i->end1(), qualsPerBase, dest );
  }

  return dest;
}

char* Alignment::writeBotQual( const uchar* qualities,
			       size_t qualsPerBase, char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // assume we're not doing translated alignment

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, i->beg1() - j->end1());

      // append qualities for unaligned chunk of bottom sequence:
      dest = writeQuals( qualities, j->end2(), i->beg2(), qualsPerBase, dest );
    }

    // append qualities for aligned chunk of bottom sequence:
    dest = writeQuals( qualities, i->beg2(), i->end2(), qualsPerBase, dest );
  }

  return dest;
}
