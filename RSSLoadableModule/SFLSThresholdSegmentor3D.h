#ifndef SFLSThresholdSegmentor3D_h_
#define SFLSThresholdSegmentor3D_h_

#include "SFLSSegmentor3D.h"


class CSFLSThresholdSegmentor3D : public CSFLSSegmentor3D
{
public:
  typedef CSFLSSegmentor3D SuperClassType;
  typedef typename SuperClassType::CSFLSLayer CSFLSLayer;

  typedef CSFLSThresholdSegmentor3D Self;


  /*================================================================================
    ctor */
  CSFLSThresholdSegmentor3D() : CSFLSSegmentor3D()
  {
    basicInit();
  }

  void basicInit();

  void setThreshold(double t) {m_threshold = t;}
  double getThreshold() { return m_threshold;}


  /* ============================================================
   * functions
   * ============================================================*/

  //void doThresholdSegmenation();
  void doSegmenation();

  /* ============================================================
     computeForce    */
  void computeForce();

private:
  double m_threshold;

};



#endif
