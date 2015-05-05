#ifndef __LASTEX__
#define __LASTEX__
#include "LastexArguments.hh"
#include "ScoreMatrix.hh"
#include "Alphabet.hh"
#include "io.hh"
#include "gumbel_params/mcf_local_alignment_evaluer.hpp"
#include "utils.hh"


double evalueForSequences( int score, double seqLen1, double seqLen2 );

void initializeEvalueCalculator(cbrc::LastexArguments args, cbrc::ScoreMatrix socreMatrix, SequenceStatistics _stats1, SequenceStatistics __stats2);

void makeEvaluer();

void __makeScoreMatrix( const std::string& matrixFile );

SequenceStatistics readStats( const std::string& fileName );

double getLambda();

double getK();

#endif

