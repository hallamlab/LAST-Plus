//
// Created by david on 20/05/15.
//

#include "inputOutput.h"
#include "lastal.hh"

std::queue<int> inputQueue;
std::queue<int> outputQueue;
bool finishedReadingFlag;

void* writerFunction(void* arguments) {

	std::queue<int> currentOutputQueue;
	int id;
	int readerCounter;
	int writerCounter;
	std::ofstream outFileStream;
	std::ostream &out = openOut(args.outFile, outFileStream);
	out.precision(3);

	while(1) {
		SEM_WAIT(writerSema);

		SEM_WAIT(inputOutputQueueSema);
		for (int j=0; j<outputQueue.size(); j++) {
			currentOutputQueue.push(outputQueue.front());
			outputQueue.pop();
		}
		SEM_POST(inputOutputQueueSema);

		while( currentOutputQueue.size() != 0 ) {
			SEM_WAIT(ioSema);

			id = currentOutputQueue.front();
			currentOutputQueue.pop();
			threadData *data = threadDatas->at(id);
			for (int j=0; j<data->outputVector->size(); j++){
				out << data->outputVector->at(j);
			}

			SEM_POST(ioSema);
			SEM_POST(data->writeSema);
		}

		SEM_WAIT(inputOutputQueueSema);
		readerCounter = inputQueue.size();
		writerCounter = outputQueue.size();
		SEM_POST(inputOutputQueueSema);

		if(finishedReadingFlag == 1 && writerCounter == 0 && readerCounter == 0)   {
			SEM_POST(terminationSema);
			break;
		}
	}
}

void readerFunction(char** argv){

	std::queue<int> currentInputQueue;
	int id;
	int count;
	char defaultInputName[] = "-";
	char *defaultInput[] = {defaultInputName, 0};
	char **inputBegin = argv + args.inputStart;

	for (char **i = *inputBegin ? inputBegin : defaultInput; *i; ++i) {
		std::ifstream inFileStream;
		std::istream &in = openIn(*i, inFileStream);

		if (args.outputType == 0) {
			matchCounts.clear();
			matchCounts.resize(query.finishedSequences());
		}
		if (volumes + 1 == 0) volumes = 1;
		for (unsigned i = 0; i < volumes; ++i) {
			if (text.unfinishedSize() == 0 || volumes > 1) readVolume(i);

			while (in) {
				SEM_WAIT(readerSema);

				SEM_WAIT(inputOutputQueueSema);
				for (int j = 0; j < inputQueue.size(); j++) {
					currentInputQueue.push(inputQueue.front());
					inputQueue.pop();
				}
				SEM_POST(inputOutputQueueSema);

				while (currentInputQueue.size() != 0) {
					SEM_WAIT(ioSema);

					count = 0;
					id = currentInputQueue.front();
					currentInputQueue.pop();
					threadData *data = threadDatas->at(id);
					// read in the data
					while (in || count > 10000) {
						data->appendFromFasta(in);
					}
					SEM_POST(ioSema);

					SEM_POST(data->readSema);
				}
			}
		}
		finishedReadingFlag = 1;
		SEM_WAIT(terminationSema);
	}
}


void* threadFunction(void *args) {

	struct threadData *data = (struct threadData*)args;

	while(1) {
		SEM_WAIT(data->readSema);
		SEM_WAIT(data->writeSema);

		data->scanAllVolumes(volumes);
		data->query.reinitForAppending();

		SEM_WAIT(inputOutputQueueSema);
		inputQueue.push( data->identifier );
		outputQueue.push( data->identifier );
		SEM_POST(inputOutputQueueSema);

		SEM_POST(readerSema);
		SEM_POST(writerSema);
	}
	return (void*) 0;
}
