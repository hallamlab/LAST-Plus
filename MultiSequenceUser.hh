#ifndef MULTISEQUENCEUSER_HH
#define MULTISEQUENCEUSER_HH

#include "ScoreMatrixRow.hh"
#include "VectorOrMmap.hh"
#include "MultiSequence.hh"

#include <string>
#include <iosfwd>

namespace cbrc{

class MultiSequenceUser{
 public:
  typedef unsigned indexT;
  typedef unsigned char uchar;

  indexT whichSequence( indexT coordinate, const MultiSequence &which ) const;

  indexT seqBeg( indexT seqNum, const MultiSequence &which ) const;

	indexT seqEnd( indexT seqNum, const MultiSequence &which ) const;

	indexT seqLen( indexT seqNum, const MultiSequence &which ) const;

	std::string seqName( indexT seqNum, const MultiSequence &which ) const;

  const uchar* seqReader(const MultiSequence &which) const;

	const uchar* qualityReader(const MultiSequence &which) const;

};

}
#endif
