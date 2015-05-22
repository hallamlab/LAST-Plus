#ifndef __SEMAPHORES_HH
#define __SEMAPHORES_HH

#include <semaphore.h>

#ifdef __APPLE__
  typedef sem_t* SEM_T;
  #define SEM_POST(x) sem_post(x)
  #define SEM_WAIT(x) sem_wait(x)
#elif __linux 
  typedef sem_t SEM_T;
  #define SEM_POST(x) sem_post(&x)
  #define SEM_WAIT(x) sem_wait(&x)
#endif


#endif
