cmake_minimum_required(VERSION 3.2)
project(ThreadedLAST)

set(DIRECTORIES gumbel_params gumbel_params/common
gumbel_params/corelib)

include_directories(${DIRECTORIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -w -m64 -pthread")

set(LASTDB_SOURCE_FILES
    Alphabet.cc
    Alphabet.hh
    MultiSequence.cc
    ScoreMatrixRow.hh
    CyclicSubsetSeed.cc
    SubsetSuffixArray.cc
    LastdbArguments.cc
    LastdbArguments.hh
    SequenceFormat.hh
    io.cc
    io.hh
    fileMap.cc
    SubsetSuffixArraySort.cc
    SubsetSuffixArray.hh
    CyclicSubsetSeed.hh
    Mmap.hh
    MultiSequenceQual.cc
    MultiSequence.hh
    stringify.hh
    lastdb.cc
    VectorOrMmap.hh
    fileMap.hh
    qualityScoreUtil.hh
    version.hh
)

set(LASTAL_SOURCE_FILES
	externalsort.hh
	externalsort.cc
	utilities.hh
	utilities.cc
	heapsort.cc
	heapsort.hh
	linereader.hh
	linereader.cc
 Alphabet.cc
 Alphabet.hh
 MultiSequence.cc
 MultiSequence.hh
 fileMap.hh
 io.hh
 CyclicSubsetSeed.cc
 SubsetSuffixArray.cc
 SubsetSuffixArray.hh
 VectorOrMmap.hh
 LastalArguments.cc
 LastalArguments.hh
 SequenceFormat.hh
 io.cc
 fileMap.cc
 ScoreMatrix.cc
 ScoreMatrix.hh
 DiagonalTable.cc
 DiagonalTable.hh
 SegmentPair.cc
 Alignment.cc
 Alphabet.hh
 GeneralizedAffineGapCosts.hh
 GeneticCode.hh
 GappedXdropAligner.cc
 ScoreMatrixRow.hh
 SegmentPairPot.cc
 SegmentPairPot.hh
 SubsetSuffixArray.cc
 SubsetSuffixArray.hh
 CyclicSubsetSeed.hh
 VectorOrMmap.hh
 Mmap.hh
 fileMap.hh
 AlignmentPot.cc
 AlignmentPot.hh
 Alignment.hh
 GeneralizedAffineGapCosts.cc
 Centroid.cc Centroid.hh GappedXdropAligner.hh
 GeneralizedAffineGapCosts.hh SegmentPair.hh
 GappedXdropAlignerInl.hh
 LambdaCalculator.cc LambdaCalculator.hh
 lambda_calculator.hh lambda_calculator.cc
 TwoQualityScoreMatrix.cc
 TwoQualityScoreMatrix.hh
 OneQualityScoreMatrix.cc
 OneQualityScoreMatrix.hh
 QualityPssmMaker.cc QualityPssmMaker.hh
  qualityScoreUtil.hh
 GeneticCode.cc
 gaplessXdrop.cc gaplessXdrop.hh
 io.cc io.hh
 gaplessPssmXdrop.cc gaplessPssmXdrop.hh
 gaplessTwoQualityXdrop.cc
 AlignmentWrite.cc
 stringify.hh
 MultiSequenceQual.cc MultiSequence.hh
 GappedXdropAlignerPssm.cc
 GappedXdropAligner2qual.cc
 GappedXdropAligner3frame.cc
 lastal.cc LastalArguments.hh SequenceFormat.hh
 SubsetSuffixArray.hh
 CyclicSubsetSeed.hh Mmap.hh
 AlignmentPot.hh SegmentPairPot.hh ScoreMatrix.hh
 DiagonalTable.hh
 gaplessTwoQualityXdrop.hh version.hh
 nrutil.hh
 nrutil.cc ludcmp.cc lubksb.cc
 LastexArguments.cc
 lastex.cc LastexArguments.hh
 gumbel_params/mcf_local_alignment_evaluer.hpp
 gumbel_params/sls_pvalues.hpp gumbel_params/sls_normal_distr_array.hpp
 version.hh

 gumbel_params/mcf_local_alignment_evaluer.cpp
 gumbel_params/mcf_local_alignment_evaluer.hpp
 gumbel_params/sls_pvalues.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/corelib/ncbitype.h
 gumbel_params/corelib/ncbi_limits.h
 gumbel_params/sls_normal_distr_array.hpp
 gumbel_params/njn_localmaxstatmatrix.hpp
 gumbel_params/njn_localmaxstat.hpp gumbel_params/sls_alp_sim.hpp
 gumbel_params/sls_alp_data.hpp gumbel_params/util/random_gen.hpp
 gumbel_params/corelib/ncbistd.hpp gumbel_params/corelib/ncbidbg.hpp
 gumbel_params/sls_alp_regression.hpp gumbel_params/sls_alp.hpp
 gumbel_params/njn_dynprogprob.cpp
 gumbel_params/corelib/ncbi_limits.h gumbel_params/njn_dynprogprob.hpp
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/njn_dynprogprobproto.hpp
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
 gumbel_params/njn_memutil.hpp gumbel_params/njn_ioutil.hpp
 gumbel_params/njn_dynprogproblim.cpp
 gumbel_params/ncbi_pch.hpp gumbel_params/njn_dynprogproblim.hpp
 gumbel_params/corelib/ncbitype.h gumbel_params/corelib/ncbi_limits.h
 gumbel_params/njn_dynprogprob.hpp gumbel_params/corelib/ncbidbg.hpp
 gumbel_params/njn_dynprogprobproto.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/njn_memutil.hpp
 gumbel_params/njn_ioutil.hpp
 gumbel_params/njn_dynprogprobproto.cpp
 gumbel_params/ncbi_pch.hpp gumbel_params/njn_dynprogprobproto.hpp
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
 gumbel_params/njn_ioutil.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/corelib/ncbistre.hpp gumbel_params/njn_ioutil.hpp
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
 gumbel_params/njn_localmaxstat.cpp
 gumbel_params/ncbi_pch.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/sls_alp_data.hpp
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h
 gumbel_params/sls_alp_regression.hpp gumbel_params/corelib/ncbi_limits.h
 gumbel_params/njn_localmaxstat.hpp gumbel_params/njn_memutil.hpp
 gumbel_params/njn_ioutil.hpp gumbel_params/njn_dynprogproblim.hpp
 gumbel_params/njn_dynprogprob.hpp gumbel_params/njn_dynprogprobproto.hpp
 gumbel_params/njn_function.hpp gumbel_params/njn_doubletype.hpp
 gumbel_params/njn_integer.hpp gumbel_params/njn_localmaxstatutil.hpp
 gumbel_params/njn_matrix.hpp gumbel_params/njn_approx.hpp
 gumbel_params/njn_vector.hpp
 gumbel_params/njn_localmaxstatmatrix.cpp
 gumbel_params/ncbi_pch.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h
 gumbel_params/njn_localmaxstatmatrix.hpp
 gumbel_params/njn_localmaxstat.hpp gumbel_params/corelib/ncbitype.h
 gumbel_params/njn_localmaxstatutil.hpp
 gumbel_params/corelib/ncbi_limits.h gumbel_params/njn_matrix.hpp
 gumbel_params/njn_approx.hpp gumbel_params/njn_doubletype.hpp
 gumbel_params/njn_ioutil.hpp gumbel_params/njn_vector.hpp
 gumbel_params/njn_memutil.hpp
 gumbel_params/njn_localmaxstatutil.cpp
 gumbel_params/ncbi_pch.hpp gumbel_params/njn_localmaxstatutil.hpp
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
 gumbel_params/corelib/ncbitype.h gumbel_params/corelib/ncbi_limits.h
 gumbel_params/njn_matrix.hpp gumbel_params/njn_approx.hpp
 gumbel_params/njn_doubletype.hpp gumbel_params/njn_ioutil.hpp
 gumbel_params/njn_vector.hpp gumbel_params/njn_dynprogproblim.hpp
 gumbel_params/njn_dynprogprob.hpp gumbel_params/corelib/ncbidbg.hpp
 gumbel_params/njn_dynprogprobproto.hpp gumbel_params/njn_integer.hpp
 gumbel_params/njn_memutil.hpp gumbel_params/njn_root.hpp
 gumbel_params/njn_function.hpp gumbel_params/sls_alp_data.hpp
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp
 gumbel_params/sls_alp_regression.hpp
 gumbel_params/njn_matrix.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/njn_matrix.hpp gumbel_params/corelib/ncbitype.h
 gumbel_params/njn_approx.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/njn_doubletype.hpp
 gumbel_params/njn_ioutil.hpp gumbel_params/njn_vector.hpp
 gumbel_params/random_gen.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h
 gumbel_params/sls_alp.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/sls_alp.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/sls_alp_data.hpp
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h
 gumbel_params/sls_alp_regression.hpp gumbel_params/corelib/ncbi_limits.h
gumbel_params/sls_alp_data.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/sls_alp_data.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/util/random_gen.hpp
 gumbel_params/corelib/ncbistd.hpp gumbel_params/corelib/ncbidbg.hpp
 gumbel_params/corelib/ncbitype.h gumbel_params/sls_alp_regression.hpp
 gumbel_params/corelib/ncbi_limits.h
gumbel_params/sls_alp_regression.cpp
 gumbel_params/ncbi_pch.hpp gumbel_params/sls_alp_regression.hpp
 gumbel_params/corelib/ncbistl.hpp gumbel_params/common/ncbi_export.h
 gumbel_params/corelib/ncbitype.h gumbel_params/corelib/ncbi_limits.h
 gumbel_params/sls_alp_data.hpp gumbel_params/util/random_gen.hpp
 gumbel_params/corelib/ncbistd.hpp gumbel_params/corelib/ncbidbg.hpp
gumbel_params/sls_alp_sim.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/sls_alp_sim.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/sls_alp_data.hpp
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/corelib/ncbitype.h
 gumbel_params/sls_alp_regression.hpp gumbel_params/corelib/ncbi_limits.h
 gumbel_params/sls_alp.hpp
 gumbel_params/sls_pvalues.cpp gumbel_params/ncbi_pch.hpp
 gumbel_params/sls_pvalues.hpp gumbel_params/corelib/ncbistl.hpp
 gumbel_params/common/ncbi_export.h gumbel_params/corelib/ncbitype.h
 gumbel_params/corelib/ncbi_limits.h
 gumbel_params/sls_normal_distr_array.hpp gumbel_params/sls_alp_data.hpp
 gumbel_params/util/random_gen.hpp gumbel_params/corelib/ncbistd.hpp
 gumbel_params/corelib/ncbidbg.hpp gumbel_params/sls_alp_regression.hpp

 utils.hh utils.cc Alphabet.cc Alphabet.hh MultiSequence.hh MultiSequence.cc
 io.hh io.cc

  semaphores.hh
  SubsetSuffixArrayUser.cc SubsetSuffixArrayUser.hh
  MultiSequenceUser.cc MultiSequenceUser.hh
)

#SPOBJ = Alphabet.o MultiSequence.o fileMap.o split/cbrc_linalg.o	\
#split/last-split.o split/cbrc_split_aligner.o split/last-split-main.o	\
#split/cbrc_unsplit_alignment.o
#
#MBOBJ = last-merge-batches.o

add_executable(lastdb ${LASTDB_SOURCE_FILES})
add_executable(lastal ${LASTAL_SOURCE_FILES})