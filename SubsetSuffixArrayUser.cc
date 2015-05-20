#include "SubsetSuffixArrayUser.hh"

using namespace cbrc;

void SubsetSuffixArrayUser::match( const indexT*& beg, const indexT*& end,
                               const uchar* queryPtr, const uchar* text,
                               indexT maxHits, indexT minDepth, SubsetSuffixArray &which ) const{


  indexT depth = 0;
  const uchar* subsetMap = which.seed.firstMap();

  // match using buckets:
  indexT bucketDepth = which.maxBucketPrefix();
  const indexT* bucketPtr = &(which.buckets[0]);
  indexT bucketBeg = 0;
  indexT bucketEnd = which.index.size();

  while( depth < bucketDepth ){
    if( bucketEnd - bucketBeg <= maxHits || depth >= minDepth ){
      beg = &(which.index[0]) + bucketBeg;
      end = &(which.index[0]) + bucketEnd;
      return;
    }
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset < CyclicSubsetSeed::DELIMITER ){
      ++depth;
      indexT step = which.bucketSteps[depth];
      bucketPtr += subset * step;
      bucketBeg = *bucketPtr;
      bucketEnd = *(bucketPtr + step);
      subsetMap = which.seed.nextMap( subsetMap + 64 );
    }else{  // we hit a delimiter in the query, so finish without any matches:
      bucketBeg = bucketEnd;
    }
  }

  // match using binary search:
  beg = &(which.index[0]) + bucketBeg;
  end = &(which.index[0]) + bucketEnd;

  while( true ){
    if( indexT(end - beg) <= maxHits || depth >= minDepth ) return;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset < CyclicSubsetSeed::DELIMITER ){
      equalRange( beg, end, text+depth, subsetMap, subset );
      ++depth;
      subsetMap = which.seed.nextMap( subsetMap + 64 );
    } else {  // we hit a delimiter in the query, so finish without any matches:
      beg = end;
    }
  }
}

void SubsetSuffixArrayUser::countMatches( std::vector<unsigned long long>& counts,
				      const uchar* queryPtr,
				      const uchar* text, SubsetSuffixArray &which ) const{

  indexT depth = 0;
  const uchar* subsetMap = which.seed.firstMap();

  // match using buckets:
  indexT bucketDepth = which.maxBucketPrefix();
  const indexT* bucketPtr = &(which.buckets[0]);
  indexT bucketBeg = 0;
  indexT bucketEnd = which.index.size();

  while( depth < bucketDepth ){
    if( bucketBeg == bucketEnd ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += bucketEnd - bucketBeg;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    ++depth;
    indexT step = which.bucketSteps[depth];
    bucketPtr += subset * step;
    bucketBeg = *bucketPtr;
    bucketEnd = *(bucketPtr + step);
    subsetMap = which.seed.nextMap( subsetMap + 64 );
  }

  // match using binary search:
  const indexT* beg = &(which.index[0]) + bucketBeg;
  const indexT* end = &(which.index[0]) + bucketEnd;

  while( true ){
    if( beg == end ) return;
    if( counts.size() <= depth ) counts.resize( depth+1 );
    counts[depth] += end - beg;
    uchar subset = subsetMap[ queryPtr[depth] ];
    if( subset == CyclicSubsetSeed::DELIMITER ) return;
    equalRange( beg, end, text+depth, subsetMap, subset );
    ++depth;
    subsetMap = which.seed.nextMap( subsetMap + 64 );
  }
}

void SubsetSuffixArrayUser::equalRange( const indexT*& beg, const indexT*& end,
				    const uchar* textBase,
				    const uchar* subsetMap, uchar subset ){
  while( beg < end ){
    const indexT* mid = beg + std::size_t( end - beg ) / 2;
    uchar s = subsetMap[ textBase[ *mid ] ];
    if( s < subset ){
      beg = mid + 1;
    }else if( s > subset ){
      end = mid;
    }else{
      if( subset > 0 )  // this "if" is unnecessary, but makes it a bit faster
	beg = lowerBound( beg, mid, textBase, subsetMap, subset );
      end = upperBound( mid + 1, end, textBase, subsetMap, subset );
      return;
    }
  }
}

const indexT*
SubsetSuffixArrayUser::lowerBound( const indexT* beg, const indexT* end,
			       const uchar* textBase,
			       const uchar* subsetMap, uchar subset ){
  for( ;; ){
    std::size_t size = end - beg;
    if( size <= 4 ) break;  // 3,4 seem good for hg18 chr21 versus itself
    const indexT* mid = beg + size / 2;
    if( subsetMap[ textBase[ *mid ] ] < subset ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }

  while( subsetMap[ textBase[ *beg ] ] < subset ) ++beg;  // linear search

  return beg;
}

const indexT*
SubsetSuffixArrayUser::upperBound( const indexT* beg, const indexT* end,
			       const uchar* textBase,
			       const uchar* subsetMap, uchar subset ){
  for( ;; ){
    std::size_t size = end - beg;
    if( size <= 4 ) break;  // 3,4 seem good for hg18 chr21 versus itself
    const indexT* mid = beg + size / 2;
    if( subsetMap[ textBase[ *mid ] ] <= subset ){
      beg = mid + 1;
    }else{
      end = mid;
    }
  }

  while( subsetMap[ textBase[ *(end-1) ] ] > subset ) --end;  // linear search

  return end;
}
