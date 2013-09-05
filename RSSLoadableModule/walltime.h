/* wall clock */
/* see "Software Optimization for High Performance Computing" p. 135 */

#ifndef walltime_h_
#define walltime_h_


extern double GLOBAL_GTH818N_START_TIME;
extern long GLOBAL_GTH818N_base_sec;
extern long GLOBAL_GTH818N_base_usec;


double walltime( double *t0 );

void gth818nTIC();
double gth818nTOC(); // return ellapsed time


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


#endif
