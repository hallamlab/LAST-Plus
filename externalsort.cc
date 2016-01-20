#include "externalsort.hh"

// Generate random string for output in order to allow mutiple LAST+ binaries to run simultaneously on a single machine.
// Check if the directory structure already exists. If it does we need to generate a new randstr

string generate_directory_name(){
  string random_string;
  int err = -1;
  do{
    random_string = random_str(30);
    string potential_directory = "/tmp" + random_string + "LASTtemp0";
    struct stat potential_directory_stat;
    err = stat(potential_directory.c_str(), &potential_directory_stat);
  } while(err != -1);
  return random_string;
}

bool comp_lines(const LINE &lhs, const LINE &rhs) {
  if (lhs->orfid < rhs->orfid) return true;

  if (lhs->orfid == rhs->orfid) {
    if (lhs->evalue < rhs->evalue) {
      return true;
    } else if (lhs->evalue == rhs->evalue) {
      if (lhs->bitscore > rhs->bitscore) {
        return true;
      }
    }
  }
  return false;
}

void free_lines(vector<Line *> &v) {
  vector<Line *>::iterator it;
  for (it = v.begin(); it != v.end(); it++)
    delete *it;
}

/* Sort the input sequences and divide them into blocks; return the number of blocks created */
int disk_sort_file(string outputdir, 
  string tobe_sorted_file_name, 
  string sorted_file_name,
  countT chunk_size, 
  string(*key_extractor)(const string &)) {

  // Create iterator for input fasta file
  std::ifstream inputfile;
  inputfile.open(tobe_sorted_file_name.c_str(), std::ifstream::in);

  countT curr_size = 0;
  countT batch = 0;


  // The current list of sequences to sort
  vector<Line *> lines;
  string line;
  Line *lineptr;

  string randstr = generate_directory_name();

  TEMPFILES *listptr = new  TEMPFILES( "/tmp", randstr + "LASTtemp0");
  listptr->clear();
  TEMPFILES *newlistptr = new  TEMPFILES( "/tmp", randstr + "LASTtemp1");
  newlistptr->clear();

  // Split input fasta into chunks to sort individually
  while (std::getline(inputfile, line).good()) {
    string orfid = key_extractor(line);
    double evalue = evalue_extractor_from_blast(line);
    double bitscore = bit_score_extractor_from_blast(line);
    lineptr = new Line;
    lineptr->setOrfId(orfid);
    lineptr->setLine(line);
    lineptr->setEvalue(evalue);
    lineptr->setBitscore(bitscore);
    lines.push_back(lineptr);

    if (curr_size > chunk_size) {
      // Sort the vector of sequence ids/lengths
      sort(lines.begin(), lines.end(), comp_lines);
      // Write the sequences to a file
      string fname =  listptr->nextFileName() ;

      write_sorted_sequences(lines, fname);
      free_lines(lines);
      batch++;
      // Clear the variables
      curr_size = 0;
      lines.clear();
    }
    curr_size++;
  }

  // Sort remaining sequences and write to last file
  if (lines.size() > 0) {
    sort(lines.begin(), lines.end(), comp_lines);
    string fname = listptr->nextFileName() ;
    write_sorted_sequences(lines, fname);
    free_lines(lines);
    lines.clear();
  }
  std::cout << std::endl;
}


std::vector<std::string> merge_some_files(const std::vector<std::string> &mergelist, 
                                          std::vector<TEMPFILES*> &directories){

  // Recursively merge in batches of 200 until all the remaining files can be fit into one batch
  std::size_t rounds = mergelist.size() / 200 + 1;

  std::string randstr = generate_directory_name();
  TEMPFILES *fileptr = new TEMPFILES( "/tmp", randstr + "LASTtemp0");
  directories.push_back(fileptr);
  //fileptr.clear();

  for(int i=0; i<=rounds-1; i++){
    std::vector<std::string>::const_iterator it = mergelist.begin() + i*200;
    std::string next_name = fileptr->nextFileName();
    if(i == rounds-1){
      std::vector<std::string> batch(it, mergelist.end()); 
      //print_vector(batch);
      merge_sorted_files(batch, next_name);
    }else{
      std::vector<std::string> batch(it, mergelist.begin() + (i+1)*200);
      //print_vector(batch);
      merge_sorted_files(batch, next_name);
    }
  }
  return fileptr->getFileNames();
}

/* Sort the input sequences and divide them into blocks */
int disk_sort_file(string outputdir, 
    string tobe_sorted_file_name, 
    string sorted_file_name, 
    countT chunk_size, 
    string(*key_extractor)(const string &), 
    const std::vector<std::string> &mergelist) {

  std::size_t num_files = mergelist.size();
  std::vector<std::string> files;
  std::vector<TEMPFILES*> directories;

/*
  std::cout << "NUM FILES : " << mergelist.size() << std::endl;
  std::cout << "NUM ROUNDS : " << mergelist.size()/200+1 << std::endl;
*/

  if(num_files > 200){
    while(num_files > 200){
      files = merge_some_files(mergelist, directories);
      num_files = files.size();
    }
    merge_sorted_files( files, tobe_sorted_file_name );
  }else{
    merge_sorted_files( mergelist, tobe_sorted_file_name );
  }
  
  for(int i=0; i<directories.size(); i++){
    directories[i]->clear();
    delete directories[i];
  }

  return 1;
}

int merge_sorted_files(const vector<string> &filenames, string sorted_file_name) {

  vector<istream_iterator<Line> > f_its;
  istream_iterator<Line> empty_it;

  vector<ifstream *> ifstream_for_filenames;
  // Open an istream_iterator for each fasta file
  int i;
  int S = filenames.size();
  Line *curr_lines = new Line[S];

  for (i = 0; i < S; i++) {
    ifstream *f_str = new ifstream(filenames[i].c_str()); // Make sure to keep the file stream "alive" outside this loop
    if ( f_str->peek() != std::ifstream::traits_type::eof() ){
      f_its.push_back(istream_iterator<Line>(*f_str));
      ifstream_for_filenames.push_back(f_str);
    } else {
      delete f_str;
    }
  }

  vector<pair<int, Line *> > values;
  pair<int, Line *> mod_line;
  S = f_its.size();

//!!
/*
  std::cout << "SIZE : " << S << std::endl;
  std::cout << "sorted_file_name : " << sorted_file_name << std::endl;
*/

  for (i = 0; i < S; i++) {
    if (f_its[i] != empty_it) {
      curr_lines[i] = *(f_its[i]);
      mod_line.first = i;
      mod_line.second = curr_lines + i;
      values.push_back(mod_line);
    }else {
      exit(0);
    }
  }

  try{
    build_heap(S, values);
  }
  catch(...) {
    std::cout << "exception \n";
  }

  // Open the output file
  ofstream outputfile;
  outputfile.open((sorted_file_name).c_str());
  if (!outputfile.is_open()) {
    cout << "Failed to  open output  file " << sorted_file_name << " ... exiting.\n";
    TEMPFILES *listptr = new  TEMPFILES( "/tmp", "LASTtemp0");
    listptr->clear();
    TEMPFILES *newlistptr = new  TEMPFILES( "/tmp", "LASTtemp1");
    newlistptr->clear();
    exit(0);
  }

  Line line;
  int iter_id;

  while (S > 0) {
    iter_id = values[0].first;
    // Get the top item off the heap (the longest remaining sequence in any file)
    line = *(f_its[iter_id]);
    line.print(outputfile);
    f_its[iter_id]++;

    // Add the next sequence to the top of the heap
    if (f_its[iter_id] != empty_it) {
      curr_lines[iter_id] = *(f_its[iter_id]);
      mod_line.first = iter_id;
      mod_line.second = curr_lines + iter_id;

      values[0] = mod_line;
    } else {
      values[0] = values[S - 1];
      S = S - 1;
    }

    // Re-heapify to make the top value percolate down if it isn't the longest
    if (S > 0) {
      heapify(values, 0, S);
    }
  }
  // Close last block file
  outputfile.close();

  f_its.clear();

  delete[] curr_lines;

  for(vector<ifstream *>::iterator it = ifstream_for_filenames.begin(); it!= ifstream_for_filenames.end(); ++it){
    (*it)->close();
    delete *it;
  }

  /*
     for (i = 0; i < S; i++) {
     remove(filenames[i].c_str());
     }
     */
  return 1;
}

/* Write the given sequences to a file in the order given by ids_lengths */
void write_sorted_sequences(vector<Line *> &lines, string filename) {
  countT i;
  ofstream output;
  output.open(filename.c_str(), std::ifstream::out);

  for (i = 0; i < lines.size(); i++) { ;
    output << lines[i]->line << std::endl;
  }
  output.close();
}

/* Remove the given file */
void remove_file(string filename) {
  if (remove(filename.c_str()) != 0) {
    cout << "Error deleting file " << filename << "\n";
  }
}
