#ifndef __OUTPUTSTRUCT_HH
#define __OUTPUTSTRUCT_HH

struct outputStruct {

  std::vector < std::string > *outputVector;
  bool done;

  outputStruct(){
    outputVector = new std::vector< std::string >(); 
    done = false;
  }
};

#endif
