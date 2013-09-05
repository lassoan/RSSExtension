#ifndef SFLSThresholdSegmentor3D_hpp_
#define SFLSThresholdSegmentor3D_hpp_

// std
#include <algorithm>

// local
#include "SFLSThresholdSegmentor3D.h"


/* ============================================================
   basicInit    */
void CSFLSThresholdSegmentor3D::basicInit()
{
  SuperClassType::basicInit();
    
  m_threshold = 0.5;
}  


/* ============================================================
   computeForce    */
void CSFLSThresholdSegmentor3D::computeForce()
{
  this->m_force.clear();

  long n = this->m_lz.size();
  double* kappaOnZeroLS = new double[ n ];

  {
    long i = 0;
    for (typename CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz, ++i)
      {
        long ix = itz->SFLSNodeComponent1;
        long iy = itz->SFLSNodeComponent2;
        long iz = itz->SFLSNodeComponent3;

        kappaOnZeroLS[i] = this->computeKappa(ix, iy, iz);
      }
  }


  for (long i = 0; i < n; ++i)
    {
      this->m_force.push_back((this->m_curvatureWeight)*kappaOnZeroLS[i]);
    }

    
  delete[] kappaOnZeroLS;
}


/* ============================================================
   doSegmenation    */
void CSFLSThresholdSegmentor3D::doSegmenation()
{
  /*============================================================
   * From the initial mask, generate: 1. SFLS, 2. mp_label and
   * 3. mp_phi.      
   */
  this->initializeSFLS();

  for (unsigned int it = 0; it < this->m_numIter; ++it)
    {
      computeForce();

      this->normalizeForce();

      this->oneStepLevelSetEvolution();
    }
}


#endif
