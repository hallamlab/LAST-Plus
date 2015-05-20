//
// Created by david on 20/05/15.
//
#include <pthread.h>
#include <queue>
#include <string>
#include "semaphores.hh"

#ifndef THREADEDLAST_INPUTOUTPUT_H
#define THREADEDLAST_INPUTOUTPUT_H

void* writer_func(void* args);
void conductWork();
void* thread_func(void *_thread_datas);

#endif //THREADEDLAST_INPUTOUTPUT_H
