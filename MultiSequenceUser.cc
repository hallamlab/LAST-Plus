// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

#include "MultiSequenceUser.hh"
#include "io.hh"
#include <iterator>  // istreambuf_iterator

using namespace cbrc;

MultiSequenceUser::indexT MultiSequenceUser::whichSequence( indexT coordinate, const MultiSequence
&which )
const{

  const indexT* u = std::upper_bound( which.ends.begin(), which.ends.end(), coordinate );
  assert( u != which.ends.begin() && u != which.ends.end() );
  return u - which.ends.begin() - 1;
}

MultiSequenceUser::indexT MultiSequenceUser::seqBeg( indexT seqNum, const MultiSequence &which )
const{
	return which.ends[seqNum];
}

MultiSequenceUser::indexT MultiSequenceUser::seqEnd( MultiSequenceUser::indexT seqNum, const MultiSequence &which )
const{
	return which.ends[seqNum+1] - which.padSize;
}

MultiSequenceUser::indexT MultiSequenceUser::seqLen( MultiSequenceUser::indexT seqNum, const MultiSequence &which )
const{
	return seqEnd(seqNum, which) - which.ends[seqNum];
}

std::string MultiSequenceUser::seqName( indexT seqNum, const MultiSequence &which )
const{
  return std::string( which.names.begin() + which.nameEnds[ seqNum ],
		      which.names.begin() + which.nameEnds[ seqNum + 1 ] );
}

const MultiSequenceUser::uchar* MultiSequenceUser::seqReader(const MultiSequence &which)
const{
	return which.seq.begin();
}

const MultiSequenceUser::uchar* MultiSequenceUser::qualityReader(const MultiSequence &which)
const{
	return which.qualityScores.begin();
}
