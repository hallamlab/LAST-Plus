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

void split_seq_name(const string  &strn, std::vector<char *> &v, char *buf) {
    strcpy(buf, strn.c_str());

    if(buf[0]!='>') {
        v.push_back(buf);
        return;
    }

    char *s1 = buf+1;
    v.clear();
    v.push_back(s1);

    while(*s1 != '\0') {
        if(*s1==' ') {
            *s1 = '\0';
            v.push_back(s1+1);
            break;
        }
        s1++;
    }
}

std::string get_orf_name(std::string  &strn, std::vector<char *> &v, char *buf) {
    split(strn, v, buf, ';');
    if(v.size() == 0)  return std::string("");
    split(std::string(v[0]), v, buf, '=');
    if(v.size() < 2)  return std::string("");
    return std::string(v[1]);
}

bool matchString(const string &str, const string & stringtomatch, bool fromstart) {

    unsigned long pos = str.find(stringtomatch);
    if(fromstart &&  pos ==0 ) return true;

    if( !fromstart && pos > 0) return true;
    return false;

}

string to_string(unsigned long i) {
    char  c[100];
    char *p = c;
    int j = 0;
    char z = '0';

    while( i > 0 ) {
        if(i< 10) {
            *p='0' + i;
            p++;
            break;
        }
        else {
            j = i%10;
            i = (i - j)/10;
            *p =  z + j;
            p++;
            //*p++;
        }
    }
    *p = '\0';
    p--;

    return string(c);
}


char BUFFER[1000];

string ShortenORFId(const string &s, regex_t *r) {
    const char * p = s.c_str();
    char *buf = BUFFER;
    regmatch_t m[1];
    int nomatch = regexec(r, p, 1, m, 0);
    if (nomatch) {
        return s;
    }
    p += m[0].rm_so;
    int d = m[0].rm_eo - m[0].rm_so;
    while( d > 0) {
        *buf = *p;
        d--;
        buf++;  p++;
    }
    *buf = '\0';
    return string(BUFFER);
}

string ShortenORFId(const string &s) {

    const char * p = s.c_str();
    char BUFFER[200];
    strcpy(BUFFER, p);

    char *c;

    c = BUFFER ;

    c = c + strlen(BUFFER)-1;

    unsigned int S =0;

    while( c >= BUFFER ) {
        if(S==0) {
            if( isdigit(*c) || *c=='_') {
                if( *c=='_') S++;
            }
            else{
                return string(BUFFER);
            }
        }
        else if(S==1) {
            if( isdigit(*c))
                S++;
            else
                return string(BUFFER);
        }
        else {
            if( !isdigit(*c) ) {
                break;
            }
        }
        c--;
    }

    return string(c+1);
}


int compile_regex(regex_t * r, const char * regex_text) {
    int status = regcomp(r, regex_text, REG_EXTENDED|REG_NEWLINE);
    if (status != 0) {
        char error_message[MAX_ERROR_MSG];
        regerror (status, r, error_message, MAX_ERROR_MSG);
        printf ("Regex error compiling '%s': %s\n",
                regex_text, error_message);
        return 1;
    }
    return 0;
}

string getpattern(regex_t *r , const char *to_match, unsigned int no ) {
    /* "P" is a pointer into the string which points to the end of the previous match. */
    const char * p = to_match;
    /* "N_matches" is the maximum number of matches allowed. */ /* "M" contains the matches found. */
    regmatch_t m[100];
    char buf[1000];
    int i = no;
    int nomatch = regexec(r, p, no+1, m, 0);
    if (nomatch)  return string();
    int start;
    int finish;
    if (m[i].rm_so == -1) return string();
    start = m[i].rm_so + (p - to_match);
    finish = m[i].rm_eo + (p - to_match);
    memcpy(buf, to_match+ start, finish-start);
    buf[finish-start]='\0';
    return string(buf);
}



/*
 * Process the product field of the BLAST/LASTout.parsed file to extract taxonomy
 * contained in square brackets '[' ']'. Returns the expected taxonomy if found. Returns
 * 'no-taxonomy' otherwise.
 */
string getTaxonomyFromProduct(const char *str) {
    char buf[10000];
    const char *c, *b, *e;
    bool found = false;
    c = str;

    bool front = false;
// TODO: Basic implementaiton. Occassionally will run into problems with double taxonomies and '[[' ']]'
    while( *c!='\0') {
        // continue iterating through string until end
        if(*c=='[') {
            front = true;
            b = c+1;
            c++;
            continue;
        }
        if(*c==']' && front) {
            e = c;
            found = true;
        }
        c++;
    }

    if( found ) {
        unsigned int len = e - b;
        memcpy(buf,b, len);
        buf[len] ='\0';
        return string(buf);
    }

    return "no-taxonomy";
}

string getSEEDID(const char *str) {
    char buf[100];
    unsigned int S=0, i=0;
    const char *c;

    c = str;
    while( *c!='\0') {
        if(S==0)  
           if(*c=='[') 
              S++;
           else { 
             if( *c!=' ' || i==0 || buf[i-1]!=' ')  {
                 buf[i] =*c;
                 i++;
              }
             }
        else {
           if(*c=='[')  S++;
           if(*c==']')  S--;
        }
        c++;
    }

    
    if(i> 0 && buf[i-1]==' ') buf[i-1]='\0';
     else buf[i]='\0';
    return string(buf);
}


string getKEGGID(const char *str) {
    char buf[100];
    unsigned int S=0;
    const char *c;

    c = str;
    while( *c!='\0') {
        if(S==0)  
           if(*c=='K') buf[S++]=*c;        else  S=0;
        else if(1<=S && S<=5)  
            if( isdigit(*c)) buf[S++]=*c;  else S=0;
        else if( isspace(*c) or *c=='\n') { 
           S=7; break;
        }
        else 
            S=0;
        c++;
    }

    if( S==7 ) {
         buf[S] ='\0'; 
        return string(buf);
    }
    else {
       if(S==6 &&  *c=='\0') { 
          buf[S] ='\0'; 
          return string(buf);
       }

    }

    buf[0]='\0';
    return string(buf);
}

string getCOGID(const char *str) {
    char buf[100];
    unsigned int S=0;
    const char *c;

    c = str;
    while( *c!='\0') {
        if(S==0)  
           if(*c=='C') buf[S++]=*c;        else  S=0;
        else if(S==1)  
           if(*c=='O') buf[S++]=*c;        else S= 0;
        else if(S==2)  
           if(*c=='G') buf[S++]=*c;        else S = 0;
        else if(3<=S && S<=6)  
            if( isdigit(*c)) buf[S++]=*c;  else S=0;
        else if( isspace(*c) or *c=='\n') { 
           S= 8; break;
        }
        else 
            S=0;
        c++;
    }

    if( S==8 ) {
         buf[S] ='\0'; 
        return string(buf);
    }
    else {
       if(S==7 &&  *c=='\0') { 
          buf[S] ='\0'; 
          return string(buf);
       }

    }

    buf[0]='\0';
    return string(buf);
}

string getECNo(const char *str, unsigned int d) {
    char buf[100];
    unsigned int S=0;
    const char *c, *b, *e;
    bool found = false;

    c = str;
    while( *c!='\0') {

        if(S%2==0)  {
            if( isdigit(*c)) {
                if(S==0) b = c;
                S++;
                if(S==2*d+1) {
                    found = true;
                    e = c+1;
                }
            }
            else
                S=0;
        }
        else {// S%2 ==1
            if(S < 2*d +1 ) {
                if( isdigit(*c))
                    ;
                else if(*c=='.')
                    S++;
                else
                    S=0;
            }
            else  {//(S = 2*d +1 )
                if(!isdigit(*c)) {
                    e = c;
                    break;
                }
                else{
                    e = c+1;
                }

            }
        }
        c++;
    }

    if( found ) {
        unsigned int len = e - b;
        memcpy(buf,b, len);
        buf[len] ='\0';
        return string(buf);
    }

    buf[0]='\0';
    return string(buf);
}

string function_extractor_from_list(const string & line){
    return line;
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


unsigned int hashIntoBucket(const char *str, unsigned int index) {
    int hashValue = 0;
    char buffer[100000];

    char *first, *second;
    first  = buffer;

    char *x = buffer;
    const char *p =str;

    // Extract contig id and orf_id from string
    while( *p != '\0') {
        if( *p=='_' ) {
            *x = '\0';
            second = x +1;
        }
        else
            *x = *p;

        p++;
        x++;
    }
    *x = '\0';

    //std::cout << str <<  "  " << first << "  " << second << std::endl;

    // Find the longer of the two ids
    char *dest, *destfixed,  *src;
    dest = first;
    src = second;
    int lenextra = strlen(second) - strlen(first);
    if( lenextra > 0)  {
        dest = second;
        src = first;
    }
    else
        lenextra = -lenextra;

    destfixed = dest; // always remember the initial point of longer

    dest = dest + lenextra; // move to aligned position

    // XOR aligned bits of source and destination starting at aligned position

    while(*src != '\0') {
        *dest = (*dest ) ^ (*src);
        dest++; src++;
    }

    while( *destfixed != '\0') {
        //hashValue = (hashValue << 4) + (unsigned int)(*destfixed);
        hashValue = hashValue + (unsigned int)(*destfixed);
/*
       int hiBits = hashValue  & 0xF0000000;
       if(hiBits!=0)
          hashValue ^= hiBits>> 24;
       hashValue &= ~hiBits;
*/
        destfixed++;
    }
    return hashValue%index;
}


unsigned long long hashStringIntoBucket(const char *str, unsigned long long index) {
    int hashValue = 0;
    
    while( *str != '\0') {
      //hashValue = (hashValue << 4) + (unsigned int)(*destfixed);
      hashValue = hashValue + (unsigned int)(*str);
      int hiBits = hashValue  & 0xF0000000;
       if(hiBits!=0)
          hashValue ^= hiBits>> 24;
       hashValue &= ~hiBits;
       str++;
    }
    return hashValue%index;
}


string to_upper(const string &str) {

    char tempbuf[1000];
    
    char *buf = tempbuf;

    if( str.size() > 1000) 
       buf = new char[str.size() + 1];

    unsigned int i;
    for(i =0; i < str.size(); i++ )
       buf[i] = std::toupper(str[i]);
    buf[i]='\0';
    
    if( str.size() > 1000) 
        delete [] buf;
     

    return string(buf);
}

void topHits(std::string filename, int maxHits){

  //std::cout << "Parsing the output for only the top k : " << maxHits << std::endl; 

	int count=0;
	int location;
	std::ifstream input(filename.c_str());
	std::string tmp = filename + "_tmp";
	std::string current;
	std::string prevorfid;
	std::string currorfid;
	std::ofstream output(tmp.c_str());

	prevorfid = "";

	while(getline(input, current)){
		location = current.find_first_of("\t");
		currorfid = current.substr(0,location);

		if(!(currorfid.compare(prevorfid) == 0 || prevorfid.size()==0 ) )
      count=0;
             
	  if(count <  maxHits) 
	    output << current << "\n";

    count++;
		prevorfid = currorfid;
	}
	std::rename(tmp.c_str(), filename.c_str());

  //std::cout << "Finished parsing the output for top hits" << std::endl; 
}
