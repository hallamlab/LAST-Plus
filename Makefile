CXX=g++
CXXFLAGS= -O3 -w -pthread -m64

DBOBJ = src/Alphabet.o \
				src/MultiSequence.o \
				src/CyclicSubsetSeed.o \
				src/SubsetSuffixArray.o \
				src/LastdbArguments.o \
				src/io.o \
				src/fileMap.o  \
				src/SubsetSuffixArraySort.o \
				src/MultiSequenceQual.o \
				src/lastdb.o \
				src/utilities.o \

ALOBJ = src/lastal.o \
				src/Alphabet.o \
				src/MultiSequence.o \
				src/CyclicSubsetSeed.o \
				src/SubsetSuffixArray.o \
				src/LastalArguments.o \
				src/io.o \
				src/fileMap.o \
				src/ScoreMatrix.o  \
				src/DiagonalTable.o \
				src/SegmentPair.o \
				src/Alignment.o \
				src/GappedXdropAligner.o \
				src/SegmentPairPot.o \
				src/AlignmentPot.o \
				src/GeneralizedAffineGapCosts.o \
				src/Centroid.o \
				src/LambdaCalculator.o \
				src/TwoQualityScoreMatrix.o \
				src/OneQualityScoreMatrix.o \
				src/QualityPssmMaker.o \
				src/GeneticCode.o \
				src/gaplessXdrop.o \
				src/gaplessPssmXdrop.o \
				src/gaplessTwoQualityXdrop.o \
				src/AlignmentWrite.o \
				src/MultiSequenceQual.o \
				src/GappedXdropAlignerPssm.o \
				src/GappedXdropAligner2qual.o \
				src/GappedXdropAligner3frame.o \
				src/lambda_calculator.o \
				src/nrutil.o \
				src/LastexArguments.o \
				src/lastex.o  \
				src/utils.o  \
				src/externalsort.o \
				src/linereader.o \
				src/utilities.o \
				src/heapsort.o \
				src/tempfiles.o \
				src/ludcmp.o \
				src/lubksb.o \
				gumbel_params/mcf_local_alignment_evaluer.o \
				gumbel_params/njn_dynprogprob.o \
				gumbel_params/njn_dynprogproblim.o  \
				gumbel_params/njn_dynprogprobproto.o \
				gumbel_params/njn_ioutil.o  \
				gumbel_params/njn_localmaxstat.o \
				gumbel_params/njn_localmaxstatmatrix.o \
				gumbel_params/njn_localmaxstatutil.o \
				gumbel_params/njn_matrix.o \
				gumbel_params/random_gen.o \
				gumbel_params/sls_alp.o \
				gumbel_params/sls_alp_data.o \
				gumbel_params/sls_alp_regression.o \
				gumbel_params/sls_alp_sim.o \
				gumbel_params/sls_pvalues.o \

ALL=lastal+ lastdb+

VPATH=src:gumbel_params

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -Igumbel_params -Isrc -c -o $@ $<

%.o : %.cc
	$(CXX) $(CXXFLAGS) -Igumbel_params -Isrc -c -o $@ $<

all: $(ALL)

lastal+: $(ALOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(ALOBJ)

lastdb+: $(DBOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(DBOBJ) 

clean:
	rm -f */*.o $(ALL)
