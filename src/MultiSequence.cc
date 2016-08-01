// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

#include "MultiSequence.hh"
#include "io.hh"
#include <sstream>
#include <algorithm>  // upper_bound
#include <cassert>
#include <cctype>  // isspace
#include <iterator>  // istreambuf_iterator
#include <iostream>
#include <algorithm> //reverse

using namespace cbrc;

std::string line;

void MultiSequence::initForAppending( indexT padSizeIn ){
  padSize = padSizeIn;
  seq.v.assign( padSize, ' ' );
  ends.v.assign( 1, padSize );
  names.v.clear();
  nameEnds.v.assign( 1, 0 );
}

void MultiSequence::reinitForAppending(){
  seq.v.erase( seq.v.begin(), seq.v.begin() + ends.v.back() - padSize );
  names.v.erase( names.v.begin(),
      names.v.begin() + nameEnds.v[ finishedSequences() ] );
  ends.v.resize(1);
  nameEnds.v.resize(1);
  if( !names.v.empty() ) nameEnds.v.push_back( names.v.size() );
}

void MultiSequence::fromFiles( const std::string& baseName, indexT seqCount,
    std::size_t qualitiesPerLetter ){
  ends.m.open( baseName + ".ssp", seqCount + 1 );
  seq.m.open( baseName + ".tis", ends.m.back() );
  nameEnds.m.open( baseName + ".sds", seqCount + 1 );
  names.m.open( baseName + ".des", nameEnds.m.back() );

  padSize = ends.m[0];

  qualityScores.m.open( baseName + ".qua",
      ends.m.back() * qualitiesPerLetter );
}

void MultiSequence::closeFiles(){
  ends.m.close();
  seq.m.close();
  nameEnds.m.close();
  names.m.close();
  qualityScores.m.close();
}


void MultiSequence::toFiles( const std::string& baseName ) const{
  memoryToBinaryFile( ends.begin(), ends.end(), baseName + ".ssp" );

  memoryToBinaryFile( seq.begin(), seq.begin() + ends.back(),
      baseName + ".tis" );

  memoryToBinaryFile( nameEnds.begin(), nameEnds.begin() + ends.size(),
      baseName + ".sds" );

  memoryToBinaryFile( names.begin(),
      names.begin() + nameEnds[ finishedSequences() ],
      baseName + ".des" );

  memoryToBinaryFile( qualityScores.begin(),
      qualityScores.begin() + ends.back() * qualsPerLetter(),
      baseName + ".qua" );
}

void MultiSequence::addName( std::string& name ){
  names.v.insert( names.v.end(), name.begin(), name.end() );
  nameEnds.v.push_back( names.v.size() );
  if( nameEnds.v.back() < names.v.size() )
    throw std::runtime_error("the sequence names are too long");
}

std::istream& MultiSequence::readFastaName( std::istream& stream ){

  std::string tmp;
  size_t pos;

  line = "";
  getline( stream, line );

  pos = line.find_first_of(" \t");
  if (pos != std::string::npos){
    tmp = line.substr(0, pos);
  }else{
    tmp = line;
  }

  if( !stream ) return stream;
  addName(tmp);
  return stream;
}

std::istream& MultiSequence::appendFromFasta( std::istream& stream, indexT maxSeqLen ){

  if( isFinished() ){
    char c = '>';
    stream >> c;
    if( c != '>' )
      throw std::runtime_error("bad FASTA sequence data: missing '>'");
    readFastaName(stream);
    if( !stream ) return stream;
  }

  std::istreambuf_iterator<char> inpos(stream);
  std::istreambuf_iterator<char> endpos;
  while( inpos != endpos ){
    uchar c = *inpos;
    if( c == '>' ) break;  // we have hit the next FASTA sequence
    if( !std::isspace(c) ){
      seq.v.push_back(c);
    }
    ++inpos;
  }

  finish();
  return stream;
}
// Keep the -1 is infinity scenario.
std::istream&
MultiSequence::appendFromFastaLASTDB( std::istream& stream, indexT maxSeqLen, bool unlimited ){
  if( isFinished() ){
    char c = '>';
    stream >> c;
    if( c != '>' )
      throw std::runtime_error("bad FASTA sequence data: missing '>'");
    readFastaName(stream);
    if( !stream ) return stream;
  }

  std::istreambuf_iterator<char> inpos(stream);
  std::istreambuf_iterator<char> endpos;
  while( inpos != endpos ){
    uchar c = *inpos;
    if( c == '>' ) break;  // we have hit the next FASTA sequence
    if( !std::isspace(c) ){
      if(!unlimited){
        if( seq.v.size() >= maxSeqLen ) break;
      }
      seq.v.push_back(c);
    }
    ++inpos;
  }

  if(!unlimited){
    if( isFinishable(maxSeqLen) ) finish();
  }else{
    finish();
  }
  return stream;
}

void MultiSequence::finish(){
  assert( !isFinished() );
  seq.v.insert( seq.v.end(), padSize, ' ' );
  ends.v.push_back( seq.v.size() );
  assert( ends.v.back() == seq.v.size() );
}

void MultiSequence::unfinish(){
  assert( isFinished() );
  ends.v.pop_back();
  seq.v.erase( seq.v.end() - padSize, seq.v.end() );
}

bool MultiSequence::isFinishable( indexT maxSeqLen ) const{
  return seq.v.size() + padSize <= maxSeqLen;
}

MultiSequence::indexT MultiSequence::whichSequence( indexT coordinate ) const{

  const indexT* u = std::upper_bound( ends.begin(), ends.end(), coordinate );
  assert( u != ends.begin() && u != ends.end() );
  return u - ends.begin() - 1;
}

std::string MultiSequence::seqName( indexT seqNum ) const{
  return std::string( names.begin() + nameEnds[ seqNum ],
      names.begin() + nameEnds[ seqNum + 1 ] );
}

void MultiSequence::removeName(){
    // Remove the offensive name
    std::string offensive_name = "";
    std::size_t seq_name_length = nameEnds.v.back();
    if(nameEnds.v.size() > 1){
        seq_name_length -= nameEnds.v[nameEnds.v.size()-2];
    }
    for(std::size_t i=0; i<seq_name_length; i++){
        offensive_name += names.v.back();
        names.v.pop_back();
    }
    nameEnds.v.pop_back();

    // Print out the name of the offending sequence
    std::reverse(offensive_name.begin(), offensive_name.end());
    std::cerr << "Erroneous sequence : " << offensive_name << std::endl;
}


void MultiSequence::printOffensiveName(){
    // Remove the offensive name
    std::string offensive_name = "";
    std::size_t seq_name_length = nameEnds.v.back();
    if(nameEnds.v.size() > 1){
        seq_name_length -= nameEnds.v[nameEnds.v.size()-2];
    }
    for(std::size_t i=0; i<seq_name_length; i++){
        offensive_name += names.v[names.v.size()-1-i];
    }

    // Print out the name of the offending sequence
    std::reverse(offensive_name.begin(), offensive_name.end());
    std::cerr << "Erroneous sequence : " << offensive_name << std::endl;
}



//!! LAST+ error handling
void MultiSequence::removeLatest(){

    removeName();

    // Remove the offensive sequence from the MultiSequence object
    std::size_t seq_length = ends.v.back() - padSize;
    if(ends.v.size() > 1){
        seq_length -= ends.v[ends.v.size()-2];
    }
    for(std::size_t j=0; j<seq_length; j++){
        seq.v.pop_back();
    }
    ends.v.pop_back();
}
