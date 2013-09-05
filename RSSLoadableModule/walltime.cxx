/* wall clock */
/* see "Software Optimization for High Performance Computing" p. 135 */

#include <sys/time.h>

#include "walltime.h"

#include <iostream>

double GLOBAL_GTH818N_START_TIME;
long GLOBAL_GTH818N_base_sec = 0;
long GLOBAL_GTH818N_base_usec = 0;


double walltime( double *t0 )
{
  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  struct timezone tzp;
  static long base_sec = 0;
  static long base_usec = 0;

  (void) gettimeofday(&tp,&tzp);
  if (base_sec == 0)
    {
      base_sec = tp.tv_sec;
      base_usec = tp.tv_usec;
    }

  time = (double) (tp.tv_sec - base_sec);
  mic = (double) (tp.tv_usec - base_usec);
  time = (time + mic * mega) - *t0;

  //std::cout<<"tt = "<<time<<std::endl<<std::flush;

  return(time);
}


void gth818nTIC()
{
  double mic;
  double mega = 0.000001;
  struct timeval tp;
  struct timezone tzp;
  // static long GLOBAL_GTH818N_base_sec = 0;
  // static long base_usec = 0;

  (void) gettimeofday(&tp,&tzp);
  if (GLOBAL_GTH818N_base_sec == 0)
    {
      GLOBAL_GTH818N_base_sec = tp.tv_sec;
      GLOBAL_GTH818N_base_usec = tp.tv_usec;
    }

  GLOBAL_GTH818N_START_TIME = (double) (tp.tv_sec - GLOBAL_GTH818N_base_sec);
  mic = (double) (tp.tv_usec - GLOBAL_GTH818N_base_usec);
  GLOBAL_GTH818N_START_TIME += mic*mega;

  //std::cout<<GLOBAL_GTH818N_START_TIME<<std::endl<<std::flush;

  return;
}


double gth818nTOC()
{
  double mic, time;
  double mega = 0.000001;
  struct timeval tp;
  struct timezone tzp;
  // static long base_sec = 0;
  // static long base_usec = 0;

  (void) gettimeofday(&tp,&tzp);
  if (GLOBAL_GTH818N_base_sec == 0)
    {
      GLOBAL_GTH818N_base_sec = tp.tv_sec;
      GLOBAL_GTH818N_base_usec = tp.tv_usec;
    }

  time = (double) (tp.tv_sec - GLOBAL_GTH818N_base_sec);
  mic = (double) (tp.tv_usec - GLOBAL_GTH818N_base_usec);
  time = (time + mic * mega) - GLOBAL_GTH818N_START_TIME;

  // std::cout<<GLOBAL_GTH818N_START_TIME<<std::endl<<std::flush;
  // std::cout<<time<<std::endl<<std::flush;

  return(time);
}

/* Here give an example for usage:

int main(int argc, char** argv)
{
  double timeZero = 0.0;
  double startTime = 0.0;
  double ellapseTime = 0.0;

  startTime = walltime(&timeZero); // this is necessary
  
  // do the job want to timing

  ellapseTime = walltime(&startTime); 

  std::cout<<"wall time = "<<ellapseTime<<std::endl;

  return 0;
}
*/
