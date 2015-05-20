// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

#include "MultiSequenceUser.hh"
#include "io.hh"
#include <iterator>  // istreambuf_iterator

using namespace cbrc;

MultiSequenceUser::indexT MultiSequenceUser::whichSequence( indexT coordinate, MultiSequence &which )
const{

  const indexT* u = std::upper_bound( which.ends.begin(), which.ends.end(), coordinate );
  assert( u != which.ends.begin() && u != which.ends.end() );
  return u - which.ends.begin() - 1;
}

MultiSequenceUser::indexT MultiSequenceUser::seqBeg( indexT seqNum, MultiSequence &which )
const{
	return which.ends[seqNum];
}

MultiSequenceUser::indexT seqEnd( MultiSequenceUser::indexT seqNum, MultiSequence &which )
const{
	return which.ends[seqNum+1] - which.padSize;
}

MultiSequenceUser::indexT seqLen( MultiSequenceUser::indexT seqNum, MultiSequence &which )
const{
	return seqEnd(seqNum, which) - which.ends[seqNum];
}

std::string MultiSequenceUser::seqName( indexT seqNum, MultiSequence &which )
const{
  return std::string( which.names.begin() + which.nameEnds[ seqNum ],
		      which.names.begin() + which.nameEnds[ seqNum + 1 ] );
}

const MultiSequenceUser::uchar* MultiSequenceUser::seqReader(MultiSequence &which)
const{
	return which.seq.begin();
}
