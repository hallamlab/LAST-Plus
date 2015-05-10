#ifndef __SEMAPHORES_HH
#define __SEMAPHORES_HH

#include <semaphore.h>
#include <vector>

#ifdef MAC_SEM
  typedef *sem_t SEM_T;
  #define SEM_POST(x) sem_post(x)
  #define SEM_WAIT(x) sem_wait(X)
#else
  typedef sem_t SEM_T;
  #define SEM_POST(x) sem_post(&x)
  #define SEM_WAIT(x) sem_wait(&x)
#endif

extern std::vector<SEM_T> *outputSemaphores;
extern SEM_T ioSema;

#endif
