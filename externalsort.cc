#include "externalsort.hh"
#include <algorithm>

#define STR_BUFFER_SIZE 10

typedef Line *LINE;

bool comp_lines(const LINE &lhs, const LINE &rhs) {
	if (lhs->orfid < rhs->orfid) return true;
	if (lhs->orfid == rhs->orfid) {
		return lhs->evalue < rhs->evalue;
	}
	return false;
}

void free_lines(vector<Line *> &v) {
	vector<Line *>::iterator it;
	for (it = v.begin(); it != v.end(); it++)
		delete *it;
}

/* Sort the input sequences and divide them into blocks; return the number of blocks created */
int disk_sort_file(string outputdir, string tobe_sorted_file_name, string sorted_file_name, \
    int chunk_size, string(*key_extractor)(const string &)) {

	string sorted_fname = outputdir + "/sorted.fasta";
	string sorted_fname_tmp = sorted_fname + ".tmp.";

	// Create iterator for input fasta file
	std::ifstream inputfile;
	inputfile.open(tobe_sorted_file_name.c_str(), std::ifstream::in);

	int curr_size = 0;
	int batch = 0;
	char buffer[STR_BUFFER_SIZE];

	vector<string> filenames;
	
	// The current list of sequences to sort
	vector<Line *> lines;
	pair<string, int> orfid_num;
	string line;
	int count = 0;
	Line *lineptr;

	// Split input fasta into chunks to sort individually
	while (std::getline(inputfile, line).good()) {

		string orfid = key_extractor(line);
		float evalue = evalue_extractor_from_blast(line);
		lineptr = new Line;
		lineptr->setOrfId(orfid);
		lineptr->setLine(line);
		lineptr->setEvalue(evalue);
		lines.push_back(lineptr);

		if (curr_size > chunk_size) {
			// Sort the vector of sequence ids/lengths
			sort(lines.begin(), lines.end(), comp_lines);
			// Write the sequences to a file
			sprintf(buffer, "%d", batch);
			string fname = sorted_fname_tmp + string(buffer);
			filenames.push_back(fname);
			write_sorted_sequences(lines, fname);

			free_lines(lines);
			batch++;
			// Clear the variables
			curr_size = 0;
			lines.clear();
		}
		curr_size++;
		count++;
	}

	// Sort remaining sequences and write to last file
	if (lines.size() > 0) {

		sort(lines.begin(), lines.end(), comp_lines);

		sprintf(buffer, "%d", batch);
		string fname = sorted_fname_tmp + string(buffer);
		filenames.push_back(fname);

		write_sorted_sequences(lines, fname);
		free_lines(lines);
		lines.clear();
	}
	inputfile.close();

	// Merge the sorted files and write into blocks
	int num_blocks = merge_sorted_files_create_blocks(filenames, chunk_size, outputdir, sorted_file_name);
	
	// Remove the individual sorted files
	unsigned int i;
	for (i = 0; i < filenames.size(); i++) {
		remove_file(filenames.at(i));
	}
	filenames.clear();
	lines.clear();
	rename(sorted_file_name.c_str(), tobe_sorted_file_name.c_str());

	return num_blocks;
}

/* Merge the individual sorted files while writing them to blocks; return the number of blocks created */
int merge_sorted_files_create_blocks(vector<string> &filenames, float block_mb, string outputdir,
                                     string sorted_file_name) {
	vector<istream_iterator<Line> > f_its;
	istream_iterator<Line> empty_it;

	// Open an istream_iterator for each fasta file
	int i;
	int S = filenames.size();
	Line *curr_lines = new Line[S];;

	for (i = 0; i < S; i++) {
		ifstream *f_str = new ifstream(
				filenames[i].c_str()); // Make sure to keep the file stream "alive" outside this loop
		f_its.push_back(istream_iterator<Line>(*f_str));
	}
	
	vector<pair<int, Line *> > values;
	pair<int, Line *> mod_line;

	for (i = 0; i < S; i++) {
		if (f_its[i] != empty_it) {
			curr_lines[i] = *(f_its[i]);
			mod_line.first = i;
			mod_line.second = curr_lines + i;
			values.push_back(mod_line);
		}
	}

	build_heap(S, values);

	// Set up variables for keeping track of blocks
	int block_num = 0;
	int block_curr_bytes = 0;
	
	// Open the output file
	ofstream outputfile;

	outputfile.open((sorted_file_name).c_str());
	if (!outputfile.is_open()) {
		cout << "Could not open sorted.block.0 ... exiting.\n";
		exit(0);
	}

	Line line;
	int iter_id;

	while (S > 0) {
		iter_id = values[0].first;
		// Get the top item off the heap (the longest remaining sequence in any file)
		line = *(f_its[iter_id]);

		line.print(outputfile);
		int seqlength = line.line.size();
		block_curr_bytes += seqlength;

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

	for (i = 0; i < S; i++) {
		remove(filenames[i].c_str());
	}
	return block_num + 1;
}

/* Write the given sequences to a file in the order given by ids_lengths */
void write_sorted_sequences(vector<Line *> &lines, string filename) {
	unsigned int i;
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
