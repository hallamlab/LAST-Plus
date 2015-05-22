//
// Created by david on 20/05/15.
//
#include <pthread.h>
#include <queue>
#include <string>
#include "semaphores.hh"

#ifndef THREADEDLAST_INPUTOUTPUT_H
#define THREADEDLAST_INPUTOUTPUT_H

void* writerFunction(void* args);
void  readerFunction(char** argv);
void* threadFunction(void *_thread_datas);

#endif //THREADEDLAST_INPUTOUTPUT_H
