// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith
// BLAST-like pair-wise sequence alignment, using suffix arrays.
#include "lastal.hh"

threadData **threadDatas;
pthread_t *threads;
pthread_t writerThread;

std::queue<int> idInputQueue;
std::queue<MultiSequence*> inputQueue;

std::queue<int> idOutputQueue;
std::queue< std::vector< std::string >* > outputQueue;

SEM_T readerSema;
SEM_T writerSema;
SEM_T ioSema;
SEM_T terminationSema;
SEM_T roundSema;
SEM_T roundCheckSema;
SEM_T inputOutputQueueSema;

bool roundDone = 0;

unsigned volumes = 1;
unsigned volume = 0;
int readSequences = 0;
int doneSequences = 0;

countT refSequences = 0;
countT refLetters = 0;
countT maxRefSequences = 0;

void createStructures(std::string &matrixFile){

  alph = new Alphabet();
  geneticCode = new GeneticCode();

  scoreMatrix = new ScoreMatrix();

  if (args->maskLowercase > 0)
    oneQualityScoreMatrixMasked = new OneQualityScoreMatrix();
  if (args->maskLowercase < 3)
    oneQualityScoreMatrix = new OneQualityScoreMatrix();
  if (args->outputType > 3) {
    oneQualityExpMatrix = new OneQualityExpMatrix();
  }

  if (args->inputFormat == sequenceFormat::prb) {
    qualityPssmMaker = new QualityPssmMaker();
  }

  twoQualityScoreMatrixMasked = new TwoQualityScoreMatrix();
  twoQualityScoreMatrix = new TwoQualityScoreMatrix();

  readOuterPrj(args->lastdbName + ".prj", volumes, refSequences, refLetters);

  if (minSeedLimit > 1) {
    if (args->outputType == 0)
      ERR("can't use option -j 0: need to re-run lastdb with i <= 1");
    if (minSeedLimit > args->oneHitMultiplicity)
      ERR("can't use option -m < " + stringify(minSeedLimit) +
          ": need to re-run lastdb with i <= " +
          stringify(args->oneHitMultiplicity));
  }

  bool isMultiVolume = volumes > 1;
  args->setDefaultsFromAlphabet(alph->letters == alph->dna, alph->isProtein(),
      isCaseSensitiveSeeds, isMultiVolume);
  makeScoreMatrix(matrixFile);

  gapCosts.assign(args->gapExistCost, args->gapExtendCost,
      args->insExistCost, args->insExtendCost, args->gapPairCost);
  if (args->outputType > 0) {
    calculateScoreStatistics();
  }
  args->setDefaultsFromMatrix(lambdaCalculator->lambda());
  minScoreGapless = args->calcMinScoreGapless(refLetters, numOfIndexes);
  if (!isMultiVolume) {
    args->minScoreGapless = minScoreGapless;
  }
  if (args->outputType > 0) {
    makeQualityScorers();
  }

  if (args->isTranslated()) {
    if (alph->letters == alph->dna) { // allow user-defined alphabet
      ERR("expected protein database, but got DNA");
    }
    queryAlph->fromString(queryAlph->dna);
    if (args->geneticCodeFile.empty()) {
      geneticCode->fromString(geneticCode->standard);
    } else {
      geneticCode->fromFile(args->geneticCodeFile);
    }
    geneticCode->codeTableSet(*alph, *queryAlph);
  } else {
    queryAlph = alph;
  }
}

void threadData::prepareThreadData(int identifier){

  if (args->outputType == 0) {  // we just want match counts
    matchCounts = new std::vector< std::vector<countT> >();  
  }

  outputVectorQueue = new std::queue< std::vector<std::string>*>();
  for (int j=0; j<2; j++){
    outputVectorQueue->push(new std::vector<std::string>() );
  }

  gappedXdropAligner = new GappedXdropAligner();

  queryQueue = new std::queue<MultiSequence*>();
  MultiSequence *query1 = new MultiSequence();
  MultiSequence *query2 = new MultiSequence();

  if (args->isTranslated()) {
    query1->initForAppending(3);
    query2->initForAppending(3);
  } else {
    queryAlph = alph;
    query1->initForAppending(1);
    query2->initForAppending(1);
  }
  queryAlph->tr(query1->seqWriter(), query1->seqWriter() + query1->unfinishedSize());
  queryAlph->tr(query2->seqWriter(), query2->seqWriter() + query2->unfinishedSize());

  queryQueue->push(query1);
  queryQueue->push(query2);

  idInputQueue.push(identifier);

  centroid = new Centroid(*gappedXdropAligner);
  this->identifier = identifier;

  round = 0;

#ifdef __APPLE__
  char name[40];

  sprintf(name, "/readSema%d", identifier);
  sem_unlink(name);
  readSema = sem_open(name, O_CREAT, 0644, 0);

  sprintf(name, "/writeSema%d", identifier);
  sem_unlink(name);
  writeSema = sem_open(name, O_CREAT, 0644, 2);
#elif __linux
  sem_init(&readSema, 0, 0);
  sem_init(&writeSema, 0, 2);
#endif
}

void makeScoreMatrix(const std::string &matrixFile) {
  if (!matrixFile.empty()) {
    scoreMatrix->fromString(matrixFile);
  }
  else if (args->matchScore < 0 && args->mismatchCost < 0 && alph->isProtein()) {
    const char *b = args->isTranslated() ? "BLOSUM80" : "BLOSUM62";
    scoreMatrix->fromString(ScoreMatrix::stringFromName(b));
  }
  else {
    scoreMatrix->matchMismatch(args->matchScore, args->mismatchCost, alph->letters);
  }

  scoreMatrix->init(alph->encode);

  // If the input is a PSSM, the score matrix is not used, and its
  // maximum score should not be used.  Here, we try to set it to a
  // high enough value that it has no effect.  This is a kludge - it
  // would be nice to use the maximum PSSM score.
  if (args->inputFormat == sequenceFormat::pssm) scoreMatrix->maxScore = 10000;
  // This would work, except the maxDrops aren't finalized yet:
  // maxScore = std::max(args.maxDropGapped, args.maxDropFinal) + 1;
}

void makeQualityScorers() {
  const ScoreMatrixRow *m = scoreMatrix->caseSensitive;  // case isn't relevant
  double lambda = lambdaCalculator->lambda();
  const std::vector<double> &lp1 = lambdaCalculator->letterProbs1();
  bool isPhred1 = isPhred(referenceFormat);
  int offset1 = qualityOffset(referenceFormat);
  const std::vector<double> &lp2 = lambdaCalculator->letterProbs2();
  bool isPhred2 = isPhred(args->inputFormat);
  int offset2 = qualityOffset(args->inputFormat);

  if (referenceFormat == sequenceFormat::fasta) {
    if (isFastq(args->inputFormat)) {
      if (args->maskLowercase > 0)
        oneQualityScoreMatrixMasked->init(m, alph->size, lambda,
            &lp2[0], isPhred2, offset2,
            alph->canonical, true);
      if (args->maskLowercase < 3)
        oneQualityScoreMatrix->init(m, alph->size, lambda,
            &lp2[0], isPhred2, offset2,
            alph->canonical, false);
      if (args->outputType > 3) {
        const OneQualityScoreMatrix &m = (args->maskLowercase < 3) ?
          *oneQualityScoreMatrix : *oneQualityScoreMatrixMasked;
        oneQualityExpMatrix->init(m, args->temperature);
      }
    } else if (args->inputFormat == sequenceFormat::prb) {
      bool isMatchMismatch = (args->matrixFile.empty() && args->matchScore > 0);
      qualityPssmMaker->init(m, alph->size, lambda, isMatchMismatch,
          args->matchScore, -args->mismatchCost,
          offset2, alph->canonical);
    }
  } else {
    if (isFastq(args->inputFormat)) {
      if (args->maskLowercase > 0)
        twoQualityScoreMatrixMasked->init(m, lambda, &lp1[0], &lp2[0],
            isPhred1, offset1, isPhred2, offset2,
            alph->canonical, true);
      if (args->maskLowercase < 3)
        twoQualityScoreMatrix->init(m, lambda, &lp1[0], &lp2[0],
            isPhred1, offset1, isPhred2, offset2,
            alph->canonical, false);
      if (args->outputType > 3) {
        ERR("fastq-versus-fastq column probabilities not implemented");
      }
    } else {
      ERR("when the reference is fastq, the query must also be fastq");
    }
  }
}

void calculateScoreStatistics() {
  LOG("calculating matrix probabilities...");
  // the case-sensitivity of the matrix makes no difference here
  lambdaCalculator->calculate(scoreMatrix->caseSensitive, alph->size);
  if (lambdaCalculator->isBad()) {
    if (isQuality(args->inputFormat) ||
        (args->temperature < 0 && args->outputType > 3))
      ERR("can't get probabilities for this score matrix");
    else LOG("can't get probabilities for this score matrix");
  } else {
    LOG("lambda=" << lambdaCalculator->lambda());
  }
}

void readOuterPrj(const std::string &fileName, unsigned &volumes,
    countT &refSequences, countT &refLetters) {

  std::ifstream f(fileName.c_str());
  if (!f) ERR("can't open file: " + fileName);
  unsigned version = 0;
  unsigned tmp_maxRefSequences = 0;

  std::string line, word;
  while (getline(f, line)) {
    std::istringstream iss(line);
    getline(iss, word, '=');
    if (word == "version") iss >> version;
    if (word == "alphabet") iss >> (*alph);
    if (word == "numofsequences"){ 
      iss >> refSequences;
      if(refSequences > maxRefSequences){
        tmp_maxRefSequences = refSequences;
      }
    }
    if (word == "numofletters") iss >> refLetters;
    if (word == "maxunsortedinterval") iss >> minSeedLimit;
    if (word == "masklowercase") iss >> isCaseSensitiveSeeds;
    if (word == "sequenceformat") iss >> referenceFormat;
    if (word == "volumes") iss >> volumes;
    if (word == "numofindexes") iss >> numOfIndexes;
  }

  if (volumes == 1){
    maxRefSequences = tmp_maxRefSequences;
  }

  if (f.eof() && !f.bad()) f.clear();
  if (alph->letters.empty() || refSequences + 1 == 1 || refLetters + 1 == 1 ||
      isCaseSensitiveSeeds < 0 || referenceFormat >= sequenceFormat::prb){
    f.setstate(std::ios::failbit);
  }
  if (!f) ERR("can't read database prj file: " + fileName);
  if (version < 294 && version > 0)
    ERR("the lastdb files are old: please re-run lastdb");
}

void readInnerPrj(const std::string &fileName,
    indexT &seqCount, indexT &seqLen) {

  std::ifstream f(fileName.c_str());
  if (!f) ERR("can't open file: " + fileName);

  std::string line, word;
  while (getline(f, line)) {
    std::istringstream iss(line);
    getline(iss, word, '=');
    if (word == "numofsequences"){ 
      iss >> seqCount;
      if(seqCount > maxRefSequences){
        maxRefSequences = seqCount;
      }
    }
    if (word == "numofletters") iss >> seqLen;
    if (word == "numofindexes") iss >> numOfIndexes;
  }

  if (f.eof() && !f.bad()) f.clear();
  if (seqCount + 1 == 0 || seqLen + 1 == 0 ){
    f.setstate(std::ios::failbit);
  }
  if (!f) ERR("can't read volume prj file: " + fileName);
}

void threadData::writeCounts(std::ostream &out) {

  LOG("writing...");
  std::stringstream outstream;

  for (indexT i = 0; i < matchCounts->size(); ++i) {
    outstream << query->seqName(i) << "\n";

    for (indexT j = args->minHitDepth; j < matchCounts[i].size(); ++j) {
      outstream << j << "\t" << (*matchCounts)[i][j] << "\n";
    }

    outstream << "\n";
    out << outstream.str();
  }
}

void threadData::countMatches(char strand) {
  LOG("counting...");
  indexT seqNum = strand == '+' ? 0 : query->finishedSequences() - 1;

  for (indexT i = 0; i < query->finishedSize(); i += args->queryStep) {
    if (strand == '+') {
      for (; ;) {
        if (seqNum == query->finishedSequences()) return;
        if (query->seqEnd(seqNum) > i) break;
        ++seqNum;
      }
      // speed-up:
      if (args->minHitDepth > query->seqEnd(seqNum) - i) continue;
    }
    else {
      indexT j = query->finishedSize() - i;
      for (; ;) {
        if (seqNum + 1 == 0) return;
        if (query->seqBeg(seqNum) < j) break;
        --seqNum;
      }
      // speed-up:
      if (args->minHitDepth > j - query->seqBeg(seqNum)) continue;
    }

    for (unsigned x = 0; x < numOfIndexes; ++x)
      suffixArrays[x].countMatches((*matchCounts)[seqNum], query->seqReader() + i, text->seqReader());
  }
}

void threadData::alignGapless(SegmentPairPot &gaplessAlns, char strand) {

  Dispatcher dis(Phase::gapless, text, query, scoreMatrix, twoQualityScoreMatrix,
      twoQualityScoreMatrixMasked, referenceFormat, alph);
  DiagonalTable dt;  // record already-covered positions on each diagonal
  countT matchCount = 0, gaplessExtensionCount = 0, gaplessAlignmentCount = 0;

  for (indexT i = 0; i < query->finishedSize(); i += args->queryStep) {

    for (unsigned x = 0; x < numOfIndexes; ++x) {

      const indexT *beg;
      const indexT *end;

      suffixArrays[x].match(beg, end, dis.b + i, dis.a, args->oneHitMultiplicity, args->minHitDepth);
      matchCount += end - beg;

      // Tried: if we hit a delimiter when using contiguous seeds, then
      // increase "i" to the delimiter position.  This gave a speed-up
      // of only 3%, with 34-nt tags.

      indexT gaplessAlignmentsPerQueryPosition = 0;

      for ( /* no-op*/; beg < end; ++beg) { // loop over suffix-array matches

        if (gaplessAlignmentsPerQueryPosition == args->maxGaplessAlignmentsPerQueryPosition) break;

        indexT j = *beg;  // coordinate in the reference sequence

        if (dt.isCovered(i, j)) continue;

        int fs = dis.forwardGaplessScore(j, i);
        int rs = dis.reverseGaplessScore(j, i);
        int score = fs + rs;
        ++gaplessExtensionCount;

        // Tried checking the score after isOptimal & addEndpoint, but
        // the number of extensions decreased by < 10%, and it was
        // slower overall.
        if (score < minScoreGapless) continue;

        indexT tEnd = dis.forwardGaplessEnd(j, i, fs);
        indexT tBeg = dis.reverseGaplessEnd(j, i, rs);
        indexT qBeg = i - (j - tBeg);
        if (!dis.isOptimalGapless(tBeg, tEnd, qBeg)) continue;
        SegmentPair sp(tBeg, qBeg, tEnd - tBeg, score);

        if (args->outputType == 1) {  // we just want gapless alignments
          Alignment aln(identifier);
          aln.fromSegmentPair(sp);
          aln.write(args->scoreCutoff, args->evalueCutoff, *text, *query, strand, args->isTranslated(), *alph, args->outputFormat,
              outputVector);
        }
        else {
          gaplessAlns.add(sp);  // add the gapless alignment to the pot
        }

        ++gaplessAlignmentsPerQueryPosition;
        ++gaplessAlignmentCount;
        dt.addEndpoint(sp.end2(), sp.end1());
      }
    }
  }
  LOG("initial matches=" << matchCount);
  LOG("gapless extensions=" << gaplessExtensionCount);
  LOG("gapless alignments=" << gaplessAlignmentCount);
}

void threadData::alignGapped(AlignmentPot &gappedAlns, SegmentPairPot &gaplessAlns, Phase::Enum phase) {

  Dispatcher dis(phase, text, query, scoreMatrix, twoQualityScoreMatrix, twoQualityScoreMatrixMasked,
      referenceFormat, alph);
  indexT frameSize = args->isTranslated() ? (query->finishedSize() / 3) : 0;
  countT gappedExtensionCount = 0, gappedAlignmentCount = 0;

  // Redo the gapless extensions, using gapped score parameters.
  // Without this, if we self-compare a huge sequence, we risk getting
  // huge gapped extensions.
  for (size_t i = 0; i < gaplessAlns.size(); ++i) {
    SegmentPair &sp = gaplessAlns.items[i];

    int fs = dis.forwardGaplessScore(sp.beg1(), sp.beg2());
    int rs = dis.reverseGaplessScore(sp.beg1(), sp.beg2());
    indexT tEnd = dis.forwardGaplessEnd(sp.beg1(), sp.beg2(), fs);
    indexT tBeg = dis.reverseGaplessEnd(sp.beg1(), sp.beg2(), rs);
    indexT qBeg = sp.beg2() - (sp.beg1() - tBeg);
    sp = SegmentPair(tBeg, qBeg, tEnd - tBeg, fs + rs);

    if (!dis.isOptimalGapless(tBeg, tEnd, qBeg)) {
      SegmentPairPot::mark(sp);
    }
  }

  erase_if(gaplessAlns.items, SegmentPairPot::isMarked);

  gaplessAlns.cull(args->cullingLimitForGaplessAlignments);
  gaplessAlns.sort();  // sort by score descending, and remove duplicates

  LOG("redone gapless alignments=" << gaplessAlns.size());

  for (size_t i = 0; i < gaplessAlns.size(); ++i) {
    SegmentPair &sp = gaplessAlns.get(i);

    if (SegmentPairPot::isMarked(sp)) continue;

    Alignment aln(identifier);
    AlignmentExtras extras;  // not used
    aln.seed = sp;

    dis.shrinkToLongestIdenticalRun(aln.seed);

    // do gapped extension from each end of the seed:
    // third...
    aln.makeXdrop(*gappedXdropAligner, *centroid, dis.a, dis.b, args->globality,
        dis.m, scoreMatrix->maxScore, gapCosts, dis.d,
        args->frameshiftCost, frameSize, dis.p,
        dis.t, dis.i, dis.j, *alph, extras);
    ++gappedExtensionCount;

    if (aln.score < args->minScoreGapped) continue;

    if (!aln.isOptimal(dis.a, dis.b, args->globality, dis.m, dis.d, gapCosts,
          args->frameshiftCost, frameSize, dis.p,
          dis.t, dis.i, dis.j)) {
      // If retained, non-"optimal" alignments can hide "optimal"
      // alignments, e.g. during non-redundantization.
      continue;
    }

    gaplessAlns.markAllOverlaps(aln.blocks);
    gaplessAlns.markTandemRepeats(aln.seed, args->maxRepeatDistance);

    if (phase == Phase::final) gappedAlns.add(aln);
    else SegmentPairPot::markAsGood(sp);

    ++gappedAlignmentCount;
  }

  LOG("gapped extensions=" << gappedExtensionCount);
  LOG("gapped alignments=" << gappedAlignmentCount);
}

void threadData::alignFinish(const AlignmentPot &gappedAlns, char strand) {

  Dispatcher dis(Phase::final, text, query, scoreMatrix, twoQualityScoreMatrix,
      twoQualityScoreMatrixMasked, referenceFormat, alph);
  indexT frameSize = args->isTranslated() ? (query->finishedSize() / 3) : 0;

  if (args->outputType > 3) {
    if (dis.p) {
      LOG("exponentiating PSSM...");
      centroid->setPssm(dis.p, query->finishedSize(), args->temperature,
          *oneQualityExpMatrix, dis.b, dis.j);
    }
    else {
      centroid->setScoreMatrix(dis.m, args->temperature);
    }
    centroid->setOutputType(args->outputType);
  }

  LOG("finishing...");

  for (size_t i = 0; i < gappedAlns.size(); ++i) {
    const Alignment &aln = gappedAlns.items[i];
    if (args->outputType < 4) {
      aln.write(args->scoreCutoff, args->evalueCutoff, *text, *query, strand, args->isTranslated(), *alph, args->outputFormat, outputVector);
    }else{  // calculate match probabilities:
      Alignment probAln(identifier);
      AlignmentExtras extras;
      probAln.seed = aln.seed;
      probAln.makeXdrop(*gappedXdropAligner, *centroid,
          dis.a, dis.b, args->globality,
          dis.m, scoreMatrix->maxScore, gapCosts, dis.d,
          args->frameshiftCost, frameSize, dis.p, dis.t,
          dis.i, dis.j, *alph, extras,
          args->gamma, args->outputType);
      probAln.write(args->scoreCutoff, args->evalueCutoff, *text, *query, strand, args->isTranslated(),
          *alph, args->outputFormat, outputVector, extras);
    }
  }
}

void threadData::makeQualityPssm(bool isApplyMasking) {
  if (!isQuality(args->inputFormat) || isQuality(referenceFormat)) return;

  LOG("making PSSM...");
  query->resizePssm();

  const uchar *seqBeg = query->seqReader();
  const uchar *seqEnd = seqBeg + query->finishedSize();
  const uchar *q = query->qualityReader();
  int *pssm = *query->pssmWriter();

  if (args->inputFormat == sequenceFormat::prb) {
    qualityPssmMaker->make(seqBeg, seqEnd, q, pssm, isApplyMasking);
  }
  else {
    const OneQualityScoreMatrix &m =
      isApplyMasking ? *oneQualityScoreMatrixMasked : *oneQualityScoreMatrix;
    makePositionSpecificScoreMatrix(m, seqBeg, seqEnd, q, pssm);
  }
}

void threadData::scan(char strand) {

  if (args->outputType == 0) {  // we just want match counts
    countMatches(strand);
    return;
  }

  bool isApplyMasking = (args->maskLowercase > 0);
  makeQualityPssm(isApplyMasking);

  LOG("scanning...");

  SegmentPairPot gaplessAlns;
  alignGapless(gaplessAlns, strand);
  if (args->outputType == 1) return;  // we just want gapless alignments

  if (args->maskLowercase == 1) makeQualityPssm(false);

  AlignmentPot gappedAlns;

  if (args->maskLowercase == 2 || args->maxDropFinal != args->maxDropGapped) {
    alignGapped(gappedAlns, gaplessAlns, Phase::gapped);
    erase_if(gaplessAlns.items, SegmentPairPot::isNotMarkedAsGood);
  }
  if (args->maskLowercase == 2) makeQualityPssm(false);

  alignGapped(gappedAlns, gaplessAlns, Phase::final);

  if (args->outputType > 2) {  // we want non-redundant alignments
    gappedAlns.eraseSuboptimal();
    LOG("nonredundant gapped alignments=" << gappedAlns.size());
  }

  gappedAlns.sort();  // sort by score
  alignFinish(gappedAlns, strand);
}

void threadData::translateAndScan(char strand) {

  if (args->isTranslated()) {
    LOG("translating...");
    std::vector<uchar> translation(query->finishedSize());
    geneticCode->translate(query->seqReader(),
        query->seqReader() + query->finishedSize(), &translation[0]);

    query->swapSeq(translation);

    scan(strand);
    query->swapSeq(translation);
  }
  else scan(strand);
}

void readIndex(const std::string &baseName, indexT seqCount) {

  LOG("reading " << baseName << "...");
  text->fromFiles(baseName, seqCount, isFastq(referenceFormat));
  LOG("reading suffix " << baseName << "...");
  for (unsigned x = 0; x < numOfIndexes; ++x) {
    if (numOfIndexes > 1) {
      suffixArrays[x].fromFiles(baseName + char('a' + x), isCaseSensitiveSeeds, alph->encode);
    } else {
      suffixArrays[x].fromFiles(baseName, isCaseSensitiveSeeds, alph->encode);
    }
  }
  LOG("finished reading " << baseName << "...");
}

void readVolume(unsigned volumeNumber) {

  std::string baseName = args->lastdbName + stringify(volumeNumber);
  indexT seqCount = indexT(-1);
  indexT seqLen = indexT(-1);
  readInnerPrj(baseName + ".prj", seqCount, seqLen);
  minScoreGapless = args->calcMinScoreGapless(seqLen, numOfIndexes);
  readIndex(baseName, seqCount);
}

void threadData::reverseComplementPssm() {

  ScoreMatrixRow *beg = query->pssmWriter();
  ScoreMatrixRow *end = beg + query->finishedSize();

  while (beg < end) {
    --end;
    for (unsigned i = 0; i < scoreMatrixRowSize; ++i) {
      unsigned j = queryAlph->complement[i];
      if (beg < end || i < j) std::swap((*beg)[i], (*end)[j]);
    }
    ++beg;
  }
}

void threadData::reverseComplementQuery() {
  LOG("reverse complementing...");
  queryAlph->rc(query->seqWriter(), query->seqWriter() + query->finishedSize());
  if (isQuality(args->inputFormat)) {
    std::reverse(query->qualityWriter(),
        query->qualityWriter() + query->finishedSize() * query->qualsPerLetter());
  } else if (args->inputFormat == sequenceFormat::pssm) {
    reverseComplementPssm();
  }
}

void threadData::scanAllVolumes(){

  if (args->strand == 2 && round > 0) {
    reverseComplementQuery();
  }

  if (args->strand != 0) {
    translateAndScan('+');
  }

  if (args->strand == 2 || (args->strand == 0 && round == 0)) {
    reverseComplementQuery();
  }

  if (args->strand != 1) {
    translateAndScan('-');
  }

  //if (args->outputType == 0){
  // writeCounts(out);
  //}

  LOG("query batch done!");
}

void writeHeader(std::ostream &out) {

  out << "# LAST version " <<

#include "version.hh"
    << "\n" << "#\n";
  args->writeCommented(out);
  out << "# Reference sequences=" << refSequences << " normal letters=" << refLetters << "\n" << "#\n";

  if (args->outputType == 0) {
    out << "# length\tcount\n";
  } else {
    if (args->outputFormat != 2) {
      if (args->inputFormat != sequenceFormat::pssm || !args->matrixFile.empty()) {
        scoreMatrix->writeCommented(out);
        out << "#\n";
      }
      out << "# Coordinates are 0-based.  For - strand matches, coordinates\n";
      out << "# in the reverse complement of the 2nd sequence are used.\n";
      out << "#\n";
    }
    if (args->outputFormat == 0) {  // tabular format
      out << "# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\t"
        << "name2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n";
    } else if (args->outputFormat == 2) { //blast-like format
      out <<
        "# Q_ID\tS_ID\tIDENT\tALIGN_LEN\tMISMATCHES\tGAPS\tQ_BEG\tQ_END\tS_BEG\tS_END\tE_VAL\tBIT_SCORE\n";
    } else {  // MAF format
      out << "# name start alnSize strand seqSize alignment\n";
    }
  }
  out << "#\n";
}

std::istream& appendFromFasta(std::istream &in, MultiSequence *query) {

  indexT maxSeqLen = args->batchSize;
  if (maxSeqLen < args->batchSize) maxSeqLen = indexT(-1);
  if (query->finishedSequences() == 0) maxSeqLen = indexT(-1);

  int count = 0;

  while(count < INPUT_SIZE && !in.eof() ){

    size_t oldUnfinishedSize = query->unfinishedSize();

    if (args->inputFormat == sequenceFormat::fasta) {
      query->appendFromFasta(in, maxSeqLen);
    } else if (args->inputFormat == sequenceFormat::prb) {
      query->appendFromPrb(in, maxSeqLen, queryAlph->size, queryAlph->decode);
    } else if (args->inputFormat == sequenceFormat::pssm) {
      query->appendFromPssm(in, maxSeqLen, queryAlph->encode,
          args->maskLowercase > 1);
    } else {
      query->appendFromFastq(in, maxSeqLen);
    }
    if (!query->isFinished() && query->finishedSequences() == 0) {
      ERR("encountered a sequence that's too long");
    }
    // encode the newly-read sequence
    queryAlph->tr(query->seqWriter() + oldUnfinishedSize,
        query->seqWriter() + query->unfinishedSize());

    if (isPhred(args->inputFormat))  // assumes one quality code per letter:
      checkQualityCodes(query->qualityReader() + oldUnfinishedSize,
          query->qualityReader() + query->unfinishedSize(),
          qualityOffset(args->inputFormat));
    count++;
  }

  return in;
}

void initializeEvalueCalulator(const std::string dbPrjFile, ScoreMatrix *scoreMatrix,
    std::string dbfilePrj) {
  SequenceStatistics _stats1, _stats2;

  fastaFileSequenceStats(dbfilePrj, &_stats1);
  _stats2 = readStats(dbPrjFile);

  cbrc::LastexArguments args1(args->gapExistCost, args->gapExtendCost);
  initializeEvalueCalculator(args1, *scoreMatrix, _stats1, _stats2);

  makeEvaluer();
}

void initializeThreads() {

  threadDatas = new threadData*[args->threadNum];
  threads = new pthread_t[args->threadNum];

  for (int i = 0; i < args->threadNum; i++) {
    threadData *data = new threadData();
    threadDatas[i] = data;
  }
  pthread_t thread;
}

void initializeSemaphores() {

#ifdef __APPLE__
  sem_unlink("/inputOutputQueueSema");
  if (( inputOutputQueueSema = sem_open("/inputOutputQueueSema", O_CREAT, 0644, 1))== SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/ioSema");
  if (( ioSema = sem_open("/ioSema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/readerSema");
  if (( readerSema = sem_open("/readerSema", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/writerSema");
  if (( writerSema = sem_open("/writerSema", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/terminationSema");
  if (( terminationSema = sem_open("/terminationSema", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/roundSema");
  if (( roundSema = sem_open("/roundSema", O_CREAT, 0644, 0)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/roundCheckSema");
  if (( roundCheckSema = sem_open("/roundCheckSema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
#elif __linux
  sem_init(&readerSema, 0, 0);
  sem_init(&writerSema, 0, 0);
  sem_init(&ioSema, 0, 1);
  sem_init(&terminationSema, 0, 0);
  sem_init(&roundSema, 0, 0);
  sem_init(&roundCheckSema, 0, 1);
  sem_init(&inputOutputQueueSema, 0, 1);
#endif
}

void *writerFunction(void *arguments){

  int id;
  bool state;
  std::vector< std::string >* current;
  std::ofstream outFileStream;
  std::ostream &out = openOut(args->outFile, outFileStream);

  while (1) {
    SEM_WAIT(writerSema);

    SEM_WAIT(inputOutputQueueSema);
    current = outputQueue.front();
    outputQueue.pop();
    id = idOutputQueue.front();
    idOutputQueue.pop();
    SEM_POST(inputOutputQueueSema);

    threadData *data = threadDatas[id];

    SEM_WAIT(ioSema);
    for (int j=0; j < current->size(); j++) {
      out << current->at(j);
    }
    SEM_POST(ioSema);

    current->clear();

    SEM_WAIT(inputOutputQueueSema);
    data->outputVectorQueue->push(current);
    doneSequences++;
    SEM_POST(inputOutputQueueSema);
    SEM_POST(data->writeSema);

    SEM_WAIT(roundCheckSema);
    if (roundDone && readSequences == doneSequences && readSequences){
      SEM_POST(roundSema);
      if (volume+1 == volumes){
        SEM_POST(terminationSema);
        break;
      }
    }
    SEM_POST(roundCheckSema);
  }
}

void readerFunction( std::istream& in ){

  int id;
  MultiSequence *current;

  if (volumes == 1) {
    readIndex(args->lastdbName, refSequences);
  }

  for (unsigned i = 0; i < volumes; ++i){

    volume = i;
    LOG(i+1 << " out of " << volumes);

    SEM_WAIT(roundCheckSema);
    roundDone = 0;
    SEM_POST(roundCheckSema);

    if (text->unfinishedSize() == 0 || volumes > 1){

      SEM_WAIT(ioSema);
      readVolume(i);
      SEM_POST(ioSema);

      for(int j=0; j<args->threadNum; j++){
        threadDatas[j]->round = i;
      }
    }

    while(!in.eof() ){

      SEM_WAIT(readerSema);
      SEM_WAIT(inputOutputQueueSema);
      id = idInputQueue.front();
      idInputQueue.pop();
      threadData *data = threadDatas[id];
      current = inputQueue.front();
      inputQueue.pop();
      SEM_POST(inputOutputQueueSema);

      SEM_WAIT(ioSema);
      appendFromFasta(in, current);
      SEM_POST(ioSema);

      data->queryQueue->push(current);
      readSequences++;

      SEM_POST(data->readSema);
    }
    in.clear();
    in.seekg(0);

    SEM_WAIT(roundCheckSema);
    roundDone = 1;
    SEM_POST(roundCheckSema);

    SEM_WAIT(roundSema);
    if (i+1 == volumes){
      SEM_WAIT(terminationSema);
      continue;
    }
  }
}

void *threadFunction(void *__threadData){

  struct threadData *data = (struct threadData *) __threadData;

  if (args->outputType == 0) {
    data->matchCounts->clear();
    data->matchCounts->resize(data->query->finishedSequences());
  }

  while (1) {
    SEM_WAIT(data->readSema);
    SEM_WAIT(data->writeSema);

    data->outputVector = data->outputVectorQueue->front();
    data->outputVectorQueue->pop();
    data->query = data->queryQueue->front();
    data->queryQueue->pop();

    data->scanAllVolumes();

    data->query->reinitForAppending();

    SEM_WAIT(inputOutputQueueSema);
    idInputQueue.push(data->identifier);
    inputQueue.push(data->query);
    idOutputQueue.push(data->identifier);
    outputQueue.push(data->outputVector);
    SEM_POST(inputOutputQueueSema);

    SEM_POST(readerSema);
    SEM_POST(writerSema);
  }
  return (void *) 0;
}

void lastal(int argc, char **argv) {


  lambdaCalculator = new LambdaCalculator();
  args = new LastalArguments();

  args->fromArgs(argc, argv);
  std::string matrixFile;

  if (!args->matrixFile.empty()) {
    matrixFile = ScoreMatrix::stringFromName(args->matrixFile);
    args->fromString(matrixFile);  // read options from the matrix file
    args->fromArgs(argc, argv);  // command line overrides matrix file
  }

  initializeThreads();
  initializeSemaphores();

  createStructures(matrixFile);

  suffixArrays = new SubsetSuffixArray[numOfIndexes];
  text = new MultiSequence();
  doneSequences -= (args->threadNum)*2;

  for (int i=0; i<args->threadNum; i++){
    threadDatas[i]->prepareThreadData(i);
  }

  char defaultInputName[] = "-";
  char *defaultInput[] = {defaultInputName, 0};
  char **inputBegin = argv + args->inputStart;

  initializeEvalueCalulator(args->lastdbName + ".prj", scoreMatrix, *inputBegin);

  for (char **i = *inputBegin ? inputBegin : defaultInput; *i; ++i) {
    std::ifstream inFileStream;
    std::istream &in = openIn(*i, inFileStream);

    pthread_create(&writerThread, 0, writerFunction, 0);
    for (int j = 0; j < args->threadNum; j++) {
      pthread_create(&(threads[j]), 0, threadFunction, (void *) threadDatas[j]);
    }

    for (int j = 0; j < args->threadNum; j++) {
      for(int k=0; k<2; k++) {
        SEM_POST(threadDatas[j]->readSema);
      }
    }
    readerFunction(in);
  }

  text->closeFiles();
  for (size_t x = 0; x < numOfIndexes; x++) {
    suffixArrays[x].closeFiles();
  }

  //now sort the LAST output on the disk
  disk_sort_file(std::string("/tmp"), args->outFile, std::string(args->outFile) + std::string("sort"),
      maxRefSequences, orf_extractor_from_blast);

  // parse the top k hits from the file
  if (args->topHits < 1000 ){
    topHits(args->outFile, args->topHits);
  }
}


int main(int argc, char **argv) {

  try {
    lastal(argc, argv);
    return EXIT_SUCCESS;
  } catch (const std::bad_alloc &e) { 
    std::cerr << "lastal: memory exception\n";
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception &e) {
    std::cerr << "lastal: " << e.what() << '\n';
    return EXIT_FAILURE;
  } catch (int i) {
    return i;
  }
}
