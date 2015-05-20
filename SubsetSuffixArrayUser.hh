#ifndef SUBSET_SUFFIX_ARRAY_USER_HH
#define SUBSET_SUFFIX_ARRAY_USER_HH

#include "CyclicSubsetSeed.hh"
#include "VectorOrMmap.hh"
#include "SubsetSuffixArray.hh"

#include "io.hh"
#include <cassert>
#include <sstream>
#include <iostream>

typedef unsigned indexT;

namespace cbrc{

  class SubsetSuffixArrayUser{
    public:
      void match( const indexT*& beg, const indexT*& end,
          const uchar* queryPtr, const uchar* text,
          indexT maxHits, indexT minDepth, SubsetSuffixArray &which ) const;

      void countMatches( std::vector<unsigned long long>& counts,
          const uchar* queryPtr, const uchar* text, SubsetSuffixArray &which ) const;

    private:
     
      static void equalRange( const indexT*& beg, const indexT*& end,
          const uchar* textBase,
          const uchar* subsetMap, uchar symbol );
      static const indexT* lowerBound( const indexT* beg, const indexT* end,
          const uchar* textBase,
          const uchar* subsetMap, uchar subset );
      static const indexT* upperBound( const indexT* beg, const indexT* end,
          const uchar* textBase,
          const uchar* subsetMap, uchar subset );

      indexT defaultBucketDepth();

      void makeBucketSteps( indexT bucketDepth );
  };

}  // end namespace
#endif
