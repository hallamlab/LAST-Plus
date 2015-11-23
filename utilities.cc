#include "utilities.hh"
#include "lastal.hh"

char *split_n_pick(const string  &strn,  char *buf, char d, unsigned int n) {
  strcpy(buf, strn.c_str());

  char *v=buf;
  char *s1 = buf;
  v=s1;

  unsigned int i =0;

  while(*s1 != '\0') {
    if(*s1==d) {
      *s1 = '\0';
      i++;
      if(i > n) return  v ;
      v=s1+1;
    }
    s1++;
  }
  return v;
}

void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.push_back(s1);
  while(*s1 != '\0') {
    if(*s1==d) {
      *s1 = '\0';
      v.push_back(s1+1);
    }
    s1++;
  }
}

string random_str(const int len) {
  static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";

  string str;

  for (int i = 0; i < len; ++i) {
    str += alphanum[rand() % (sizeof(alphanum) - 1)];
  }
  return str;
}

string orf_extractor_from_blast(const string & line){
  char buf[10000];
  string orfid  = split_n_pick(line, buf, '\t', 0);
  return orfid;
}

double evalue_extractor_from_blast(const string &line){
  char buf[10000];
  string evaluestr  = split_n_pick(line, buf, '\t', 10);
  double evalue ;

  try{
    evalue = atof(evaluestr.c_str());
  }
  catch(...) {
    return 100;
  }
  return evalue;
}

double bit_score_extractor_from_blast(const string &line) {
  char buf[10000];
  string value  = split_n_pick(line, buf, '\t', 11);
  double bitscore ;

  try{
    bitscore = atof(value.c_str());
  }
  catch(...) {
    return 0;
  }
  return bitscore;
}


void topHits(std::string filename, int maxHits){
  int count=0;
  int location;
  std::ifstream input(filename.c_str());
  std::string tmp = filename + "_tmp";
  std::string current;
  std::string prevorfid;
  std::string currorfid;
  std::ofstream output(tmp.c_str());

  prevorfid = "";

  while(getline(input, current)) {
    location = current.find_first_of("\t");
    currorfid = current.substr(0,location);

    if(!(currorfid.compare(prevorfid) == 0 || prevorfid.size()==0 ) )
      count=0;

    if(count <  maxHits) {
      output << current << endl;
    }

    count++;
    prevorfid = currorfid;
  }
  std::rename(tmp.c_str(), filename.c_str());
}

// Niels' debug functions
void testFunction() {
  cout << "Hello from testFunction()" << endl;
}

void printOutSequence(const unsigned char* field, const unsigned char decode[]) {
  cout << "In printOutSequence" << endl;
  const unsigned char* s = field;
  int len = 0;
  while(*s) {
    len++;
    cout << decode[*s];
    s++;
  }
  cout << endl;
  cout << len << endl;
}
