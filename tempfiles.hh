#include <iostream>
#include <stack>
#include <stdlib.h>
#include <dirent.h>
#include <fstream>
#include <sys/stat.h>
#include <vector>
#include <cstring>

using namespace std;


typedef struct TempFiles {

      TempFiles(std::string _tempdir, std::string _basedir) : tempdir(_tempdir), basedir(_basedir), count(0), S(10) {} 
      void setFanOut(unsigned int i ) {  S = i; }



public:
     string nextFileName() ;
     void clear();
     unsigned int size();
     vector<string> getFileNames();

private:
     string filename(unsigned int i);
     string dirname(unsigned int i);
     string toString(unsigned int i);
     void remove_dir(char *path);
     vector<string> filenames;

     vector<string> names;
     std::string tempdir;
     std::string basedir ;
     unsigned int count; 
     unsigned int S; // fanout


} TEMPFILES;

