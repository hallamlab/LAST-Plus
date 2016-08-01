cmake_minimum_required(VERSION 3.2)
project(ThreadedLAST)

set(DIRECTORIES gumbel_params 
	gumbel_params/common
	gumbel_params/corelib)

include_directories(${DIRECTORIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -w -m64 -pthread")

set(LASTDB_SOURCE_FILES
	src/Alphabet.cc
	src/Alphabet.hh
	src/MultiSequence.cc
	src/ScoreMatrixRow.hh
	src/CyclicSubsetSeed.cc
	src/SubsetSuffixArray.cc
	src/LastdbArguments.cc
	src/LastdbArguments.hh
	src/SequenceFormat.hh
	src/io.cc
	src/io.hh
	src/fileMap.cc
	src/SubsetSuffixArraySort.cc
	src/SubsetSuffixArray.hh
	src/CyclicSubsetSeed.hh
	src/Mmap.hh
	src/MultiSequenceQual.cc
	src/MultiSequence.hh
	src/stringify.hh
	src/lastdb.cc
	src/VectorOrMmap.hh
	src/fileMap.hh
	src/qualityScoreUtil.hh
	src/version.hh
	)

set(LASTAL_SOURCE_FILES
	src/externalsort.hh
	src/externalsort.cc
	src/utilities.hh
	src/utilities.cc
	src/heapsort.cc
	src/heapsort.hh
	src/linereader.hh
	src/linereader.cc
	src/Alphabet.cc
	src/Alphabet.hh
	src/MultiSequence.cc
	src/MultiSequence.hh
	src/fileMap.hh
	src/io.hh
	src/CyclicSubsetSeed.cc
	src/SubsetSuffixArray.cc
	src/SubsetSuffixArray.hh
	src/VectorOrMmap.hh
	src/LastalArguments.cc
	src/LastalArguments.hh
	src/SequenceFormat.hh
	src/io.cc
	src/fileMap.cc
	src/ScoreMatrix.cc
	src/ScoreMatrix.hh
	src/DiagonalTable.cc
	src/DiagonalTable.hh
	src/SegmentPair.cc
	src/Alignment.cc
	src/Alphabet.hh
	src/GeneralizedAffineGapCosts.hh
	src/GeneticCode.hh
	src/GappedXdropAligner.cc
	src/ScoreMatrixRow.hh
	src/SegmentPairPot.cc
	src/SegmentPairPot.hh
	src/SubsetSuffixArray.cc
	src/SubsetSuffixArray.hh
	src/CyclicSubsetSeed.hh
	src/VectorOrMmap.hh
	src/Mmap.hh
	src/fileMap.hh
	src/AlignmentPot.cc
	src/AlignmentPot.hh
	src/Alignment.hh
	src/GeneralizedAffineGapCosts.cc
	src/Centroid.cc 
	src/Centroid.hh 
	src/GappedXdropAligner.hh
	src/GeneralizedAffineGapCosts.hh 
	src/SegmentPair.hh
	src/GappedXdropAlignerInl.hh
	src/LambdaCalculator.cc 
	src/LambdaCalculator.hh
	src/lambda_calculator.hh 
	src/lambda_calculator.cc
	src/TwoQualityScoreMatrix.cc
	src/TwoQualityScoreMatrix.hh
	src/OneQualityScoreMatrix.cc
	src/OneQualityScoreMatrix.hh
	src/QualityPssmMaker.cc 
	src/QualityPssmMaker.hh
	src/qualityScoreUtil.hh
	src/GeneticCode.cc
	src/gaplessXdrop.cc 
	src/gaplessXdrop.hh
	src/io.cc io.hh
	src/gaplessPssmXdrop.cc 
	src/gaplessPssmXdrop.hh
	src/gaplessTwoQualityXdrop.cc
	src/AlignmentWrite.cc
	src/stringify.hh
	src/MultiSequenceQual.cc 
	src/MultiSequence.hh
	src/GappedXdropAlignerPssm.cc
	src/GappedXdropAligner2qual.cc
	src/GappedXdropAligner3frame.cc
	src/lastal.cc 
	src/LastalArguments.hh 
	src/SequenceFormat.hh
	src/SubsetSuffixArray.hh
	src/CyclicSubsetSeed.hh 
	src/Mmap.hh
	src/AlignmentPot.hh 
	src/SegmentPairPot.hh 
	src/ScoreMatrix.hh
	src/DiagonalTable.hh
	src/gaplessTwoQualityXdrop.hh 
	src/version.hh
	src/nrutil.hh
	src/nrutil.cc 
	src/ludcmp.cc 
	src/lubksb.cc
	src/LastexArguments.cc
	src/lastex.cc 
	src/LastexArguments.hh
	src/version.hh
	src/utils.hh 
	src/utils.cc 
	src/Alphabet.cc 
	src/Alphabet.hh 
	src/MultiSequence.hh 
	src/MultiSequence.cc
	src/io.hh 
	src/io.cc
	src/semaphores.hh
	src/SubsetSuffixArrayUser.cc 
	src/SubsetSuffixArrayUser.hh
	src/MultiSequenceUser.cc 
	src/MultiSequenceUser.hh

	gumbel_params/mcf_local_alignment_evaluer.hpp
	gumbel_params/sls_pvalues.hpp 
	gumbel_params/sls_normal_distr_array.hpp
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
	)

add_executable(lastdb ${LASTDB_SOURCE_FILES})
add_executable(lastal ${LASTAL_SOURCE_FILES})
