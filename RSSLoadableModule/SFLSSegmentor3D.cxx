#ifndef SFLSSegmentor3D_hpp_
#define SFLSSegmentor3D_hpp_

#include "SFLSSegmentor3D.h"

#include <algorithm>
#include <cmath>

#include <csignal>
#include <cassert>

#include <fstream>

#include <vtkVersion.h>
#include <vtkImageData.h>

#include "vtkImageCast.h"


// dbg
#include "vtkMetaImageWriter.h"


CSFLSSegmentor3D::CSFLSSegmentor3D() : CSFLS()
{
    basicInit();
}


/* ============================================================
   basicInit    */
void
CSFLSSegmentor3D::basicInit()
{
    mp_img = 0;
    mp_label = 0;
    mp_mask = 0; // 0, non-0 mask for object
    mp_phi = 0;

    m_numIter = 100;
    m_timeStep = 1.0;

    m_nx = 0;
    m_ny = 0;
    m_nz = 0;

    m_increment0 = 0;
    m_increment1 = 0;
    m_increment2 = 0;

    mp_img_ptr = 0;
//    mp_label_ptr = 0;
//    mp_mask_ptr = 0;
//    mp_phi_ptr = 0;

    m_dx = 1.0;
    m_dy = 1.0;
    m_dz = 1.0;

    m_curvatureWeight = 0.0;

    m_insideVoxelCount = 0;
    m_insideVolume = 0;

    m_maxVolume = 1e10; // in mm^3
    m_maxRunningTime = 3600; // in sec

    m_keepZeroLayerHistory = false;

    m_currentIteration = 0;


    m_done = false;
}

/* ============================================================
   setNumIter                                                */
void CSFLSSegmentor3D::setNumIter(unsigned long n)
{
    m_numIter = n;
}

/* ============================================================
   setImage    */

void CSFLSSegmentor3D::setImage(vtkImageData* img)
{
    vtkImageCast* castFilter = vtkImageCast::New();
    castFilter->SetInput(img);
    castFilter->SetOutputScalarTypeToShort();
    castFilter->Update();

    mp_img = castFilter->GetOutput();
    /* TODO This is to avoid the cases that the starting extent is not 0, 0, 0.
This should not happen if the images are freshly loaded from file. Buy may happen when the input images are cropped from another image in memory.*/

    int* size = img->GetDimensions();

    if (m_nx + m_ny + m_nz == 0)
    {
        std::cout<<"Extent = ("<<size[0]<<", "<<size[1]<<", "<<size[2]<<")\n";
        m_nx = size[0];
        m_ny = size[1];
        m_nz = size[2];

        img->GetSpacing(m_dx, m_dy, m_dz);
    }
    else if ( m_nx != static_cast<long>(size[0]) || m_ny != static_cast<long>(size[1]) || static_cast<long>(size[2]) )
    {
        std::cerr<<"image sizes do not match. abort\n";
        raise(SIGABRT);
    }

    mp_img->GetIncrements(m_increment0, m_increment1, m_increment2);

    mp_img_ptr = static_cast<TPixel*>(mp_img->GetScalarPointer(0, 0, 0));

    return;
}


/* ============================================================
   setMaxVolume    */
void CSFLSSegmentor3D::setMaxVolume(double v)
{
    if (v <= 0)
    {
        std::cerr<<"Error: max volume >= 0\n";
        raise(SIGABRT);
    }

    m_maxVolume = v*1000; // v is in mL, m_maxVolume is in mm^3

    return;
}

/* ============================================================ */
void CSFLSSegmentor3D::setMaxRunningTime(double t)
{
    if (t <= 0)
    {
        std::cerr<<"Error: t <= 0\n";
        raise(SIGABRT);
    }

    m_maxRunningTime = t*60; // t is in min, m_maxRunningTime is in second

    return;
}


/* ============================================================
   setCurvatureWeight    */
void CSFLSSegmentor3D::setCurvatureWeight(double a)
{
    if (a < 0)
    {
        std::cerr<<"Error: curvature weight < 0\n";
        raise(SIGABRT);
    }

    m_curvatureWeight = a;

    return;
}



/* ============================================================
   setMask    */
void CSFLSSegmentor3D::setMask(vtkImageData* mask)
{
    mp_mask = mask;
    int* size = mask->GetDimensions();

    if (m_nx + m_ny + m_nz == 0)
    {
        m_nx = size[0];
        m_ny = size[1];
        m_nz = size[2];
    }
    else if ( m_nx != (long)size[0] || m_ny != (long)size[1] || m_nz != (long)size[2] )
    {
        std::cerr<<"image sizes do not match. abort\n";
        raise(SIGABRT);
    }

    return;
}


bool CSFLSSegmentor3D::getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(long ix, long iy, long iz, LabelPixelType* mp_label_this_ptr, LevelSetPixelType* mp_phi_this_ptr, LevelSetPixelType& thePhi)
{
    /*--------------------------------------------------
   *
   * Look in all the neighbors, to find the phi value of the nbhd:
   * this nbhd should satisfy: 1. its layer is strictly closer to
   * the zero layer hence its value is thought to be updated. 2. If
   * there are several nbhd's belonging to the same layer, choose
   * the one whose phi value has the smallest abs value.  If (ix,
   * iy) is outside, go through all nbhd who is in the layer of
   * label = mylevel-1 pick the SMALLEST phi. If (ix, iy) is inside,
   * go through all nbhd who is in the layer of label = mylevel+1
   * pick the LARGEST phi.
   */

    LabelPixelType mylevel = mp_label_this_ptr[0];
    bool foundNbhd = false;

    LevelSetPixelType itsPhi = 0.0;


    /*
    TODO the if branch in this may be removed by some clever constraint on the level set not evolving to the boundary region.
*/

    if (mylevel > 0)
    {
        // find the SMALLEST phi
        thePhi = 10000.0;

        if ((ix+1 < m_nx) && mp_label_this_ptr[m_increment0] == mylevel-1 )
        {
            itsPhi = mp_phi_this_ptr[m_increment0];
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ( (ix-1 >= 0) && mp_label_this_ptr[-m_increment0] == mylevel-1)
        {
            itsPhi = mp_phi_this_ptr[-m_increment0];
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ( (iy+1 < m_ny) && mp_label_this_ptr[m_increment1] == mylevel-1 )
        {
            itsPhi = mp_phi_this_ptr[m_increment1];
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ((iy-1 >= 0) && mp_label_this_ptr[-m_increment1] == mylevel-1 )
        {

            itsPhi = mp_phi_this_ptr[-m_increment1];
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ((iz+1 < m_nz) && mp_label_this_ptr[m_increment2] == mylevel-1)
        {
            itsPhi = mp_phi_this_ptr[m_increment2];
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ( (iz-1 >= 0) && mp_label_this_ptr[m_increment2] == mylevel-1 )
        {
            itsPhi = mp_phi_this_ptr[-m_increment2];
            thePhi = thePhi<itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }
    }
    else
    {
        // find the LARGEST phi
        thePhi = -10000;

        if ( (ix+1 < m_nx) && mp_label_this_ptr[m_increment0] == mylevel+1 )
        {
            itsPhi = mp_phi_this_ptr[m_increment0];
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ((ix-1 >= 0) && mp_label_this_ptr[-m_increment0] == mylevel+1)
        {
            itsPhi = mp_phi_this_ptr[-m_increment0];
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }

        if ((iy+1 < m_ny) && mp_label_this_ptr[m_increment1] == mylevel+1)
        {
            itsPhi = mp_phi_this_ptr[m_increment1];
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }


        if ((iy-1 >= 0) && mp_label_this_ptr[-m_increment1] == mylevel+1)
        {
            itsPhi = mp_phi_this_ptr[-m_increment1];
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }


        if ((iz+1 < m_nz) && mp_label_this_ptr[m_increment2] ==  mylevel+1)
        {
            itsPhi = mp_phi_this_ptr[m_increment2];
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }


        if ((iz-1 >= 0) && mp_label_this_ptr[-m_increment2] == mylevel+1)
        {
            itsPhi = mp_phi_this_ptr[-m_increment2];
            thePhi = thePhi>itsPhi?thePhi:itsPhi;

            foundNbhd = true;
        }
    }

    return foundNbhd;
}


/* ============================================================
   normalizeForce
   Normalize m_force s.t. max(abs(m_force)) < 0.5 */
void CSFLSSegmentor3D::normalizeForce()
{
    unsigned long nLz = m_lz.size();

    if (m_force.size() != nLz )
    {
        std::cerr<<"m_force.size() = "<<m_force.size()<<std::endl;
        std::cerr<<"nLz = "<<nLz<<std::endl;

        std::cerr<<"m_force.size() != nLz, abort.\n";
        raise(SIGABRT);
    }

    LevelSetPixelType fMax = fabs( m_force.front() );

    {
        long nf = m_force.size();

        for (long itf = 0; itf < nf; ++itf)
        {
            float v = fabs(m_force[itf]);
            fMax = fMax>v?fMax:v;
        }
    }
    fMax /= 0.49;

    {
        long nf = m_force.size();

        for (long itf = 0; itf < nf; ++itf)
        {
            m_force[itf] /= (fMax + 1e-10);
        }
    }
}


/* ============================================================
   updateInsideVoxelCount    */
void CSFLSSegmentor3D::updateInsideVoxelCount()
{
    m_insideVoxelCount -= m_lIn2out.size();
    m_insideVoxelCount += m_lOut2in.size();

    m_insideVolume = m_insideVoxelCount*m_dx*m_dy*m_dz;

    //dbg
    std::cout<<"m_insideVolume = "<<m_insideVolume<<std::endl;
    //dbg, end

    return;
}



/* ============================================================
   oneStepLevelSetEvolution    */
void CSFLSSegmentor3D::oneStepLevelSetEvolution()
{
    // create 'changing status' lists
    CSFLSLayer Sz;
    CSFLSLayer Sn1;
    CSFLSLayer Sp1;
    CSFLSLayer Sn2;
    CSFLSLayer Sp2;

    m_lIn2out.clear();
    m_lOut2in.clear();

    /*--------------------------------------------------
    1. add F to phi(Lz), create Sn1 & Sp1
    scan Lz values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========                */
    {
        long nz = m_lz.size();
        std::vector<CSFLSLayer::iterator> m_lzIterVct( nz );
        {
            long iiizzz = 0;
            for (CSFLSLayer::iterator itz = m_lz.begin(); itz != m_lz.end(); ++itz)
                m_lzIterVct[iiizzz++] = itz;
        }

        for (long iiizzz = 0; iiizzz < nz; ++iiizzz)
        {
            long itf = iiizzz;

            CSFLSLayer::iterator itz = m_lzIterVct[iiizzz];

            long ix = itz->SFLSNodeComponent1;
            long iy = itz->SFLSNodeComponent2;
            long iz = itz->SFLSNodeComponent3;

            LevelSetPixelType* phi_old_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
            LevelSetPixelType phi_old = phi_old_ptr[0];

            LevelSetPixelType phi_new = phi_old + m_force[itf];

            /*----------------------------------------------------------------------
          Update the lists of pt who change the state, for faster
          energy fnal computation. */
            if ( phi_old <= 0.0 && phi_new > 0.0 )
            {
                m_lIn2out.push_back(NodeType(ix, iy, iz));
            }

            if( phi_old > 0.0  && phi_new <= 0.0)
            {
                m_lOut2in.push_back(NodeType(ix, iy, iz));
            }

            phi_old_ptr[0] = phi_new;

            if(phi_new > 0.5)
            {
                Sp1.push_back(*itz);
                itz = m_lz.erase(itz);
            }
            else if (phi_new < -0.5)
            {
                Sn1.push_back(*itz);
                itz = m_lz.erase(itz);
            }
            else
            {
                ++itz;
            }
            /*--------------------------------------------------
          NOTE, mp_label are (should) NOT update here. They should
          be updated with Sz, Sn/p's
          --------------------------------------------------*/
        }
    }

    /*--------------------------------------------------
    2. update Ln1,Lp1,Lp2,Lp2, ****in that order****

    2.1 scan Ln1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ==========                     */
    for (CSFLSLayer::iterator itn1 = m_ln1.begin(); itn1 != m_ln1.end(); )
    {
        long ix = itn1->SFLSNodeComponent1;
        long iy = itn1->SFLSNodeComponent2;
        long iz = itn1->SFLSNodeComponent3;

        LevelSetPixelType* phiPtr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        LabelPixelType* mp_label_this_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));

        LevelSetPixelType thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, mp_label_this_ptr, phiPtr, thePhi);


        if (found)
        {
            LevelSetPixelType phi_new = thePhi - 1.0;

            phiPtr[0] = phi_new;

            if (phi_new >= -0.5)
            {
                Sz.push_back(*itn1);
                itn1 = m_ln1.erase(itn1);
            }
            else if (phi_new < -1.5)
            {
                Sn2.push_back(*itn1);
                itn1 = m_ln1.erase(itn1);
            }
            else
            {
                ++itn1;
            }
        }
        else
        {
            /*--------------------------------------------------
            No nbhd in inner (closer to zero contour) layer, so
            should go to Sn2. And the phi shold be further -1
          */
            Sn2.push_back(*itn1);
            itn1 = m_ln1.erase(itn1);

            phiPtr[0] -= 1.0;
        }
    }


    /*--------------------------------------------------
    2.2 scan Lp1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========          */
    for (CSFLSLayer::iterator itp1 = m_lp1.begin(); itp1 != m_lp1.end();)
    {
        long ix = itp1->SFLSNodeComponent1;
        long iy = itp1->SFLSNodeComponent2;
        long iz = itp1->SFLSNodeComponent3;

        LevelSetPixelType* phiPtr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        LabelPixelType* mp_label_this_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));

        LevelSetPixelType thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, mp_label_this_ptr, phiPtr, thePhi);

        if (found)
        {
            LevelSetPixelType phi_new = thePhi + 1.0;
            phiPtr[0] = phi_new;

            if (phi_new <= 0.5)
            {
                Sz.push_back(*itp1);
                itp1 = m_lp1.erase(itp1);
            }
            else if (phi_new > 1.5)
            {
                Sp2.push_back(*itp1);
                itp1 = m_lp1.erase(itp1);
            }
            else
            {
                ++itp1;
            }
        }
        else
        {
            /*--------------------------------------------------
            No nbhd in inner (closer to zero contour) layer, so
            should go to Sp2. And the phi shold be further +1
          */

            Sp2.push_back(*itp1);
            itp1 = m_lp1.erase(itp1);

            phiPtr[0] += 1.0;
        }
    }



    /*--------------------------------------------------
    2.3 scan Ln2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ==========                                      */
    for (CSFLSLayer::iterator itn2 = m_ln2.begin(); itn2 != m_ln2.end(); )
    {
        long ix = itn2->SFLSNodeComponent1;
        long iy = itn2->SFLSNodeComponent2;
        long iz = itn2->SFLSNodeComponent3;

        LevelSetPixelType* phiPtr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        LabelPixelType* mp_label_this_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));
        LabelPixelType* mp_label_mask_pixel_ptr = static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz));

        LevelSetPixelType thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, mp_label_this_ptr, phiPtr, thePhi);

        if (found)
        {
            LevelSetPixelType phi_new = thePhi - 1.0;
            phiPtr[0] = phi_new;

            if (phi_new >= -1.5)
            {
                Sn1.push_back(*itn2);
                itn2 = m_ln2.erase(itn2);
            }
            else if (phi_new < -2.5)
            {
                itn2 = m_ln2.erase(itn2);
                phiPtr[0] = -3.0;
                mp_label_this_ptr[0] = -3;
                mp_label_mask_pixel_ptr[0] = 1;
            }
            else
            {
                ++itn2;
            }
        }
        else
        {
            itn2 = m_ln2.erase(itn2);
            phiPtr[0] = -3.0;
            mp_label_this_ptr[0] = -3;
            mp_label_mask_pixel_ptr[0] = 1;
        }
    }



    /*--------------------------------------------------
    2.4 scan Lp2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
    ========= */
    for (CSFLSLayer::iterator itp2 = m_lp2.begin(); itp2 != m_lp2.end(); )
    {
        long ix = itp2->SFLSNodeComponent1;
        long iy = itp2->SFLSNodeComponent2;
        long iz = itp2->SFLSNodeComponent3;

        LevelSetPixelType* phiPtr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        LabelPixelType* mp_label_this_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));
        LabelPixelType* mp_label_mask_pixel_ptr = static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz));


        LevelSetPixelType thePhi;
        bool found = getPhiOfTheNbhdWhoIsClosestToZeroLevelInLayerCloserToZeroLevel(ix, iy, iz, mp_label_this_ptr, phiPtr, thePhi);

        if (found)
        {
            LevelSetPixelType phi_new = thePhi + 1.0;
            phiPtr[0] = phi_new;

            if (phi_new <= 1.5)
            {
                Sp1.push_back(*itp2);
                itp2 = m_lp2.erase(itp2);
            }
            else if (phi_new > 2.5)
            {
                itp2 = m_lp2.erase(itp2);
                phiPtr[0] = 3.0;
                mp_label_this_ptr[0] = 3;
                mp_label_mask_pixel_ptr[0] = 0;

            }
            else
            {
                ++itp2;
            }
        }
        else
        {
            itp2 = m_lp2.erase(itp2);
            phiPtr[0] = 3.0;
            mp_label_this_ptr[0] = 3;
            mp_label_mask_pixel_ptr[0] = 0;

        }
    }


    /*--------------------------------------------------
    3. Deal with S-lists Sz,Sn1,Sp1,Sn2,Sp2
    3.1 Scan Sz */
    for (CSFLSLayer::iterator itSz = Sz.begin(); itSz != Sz.end(); ++itSz)
    {
        long ix = itSz->SFLSNodeComponent1;
        long iy = itSz->SFLSNodeComponent2;
        long iz = itSz->SFLSNodeComponent3;

        m_lz.push_back(*itSz);
        *(static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz))) = 0;
        *(static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz))) = 1;

    }


    /*--------------------------------------------------
    3.2 Scan Sn1     */
    for (CSFLSLayer::iterator itSn1 = Sn1.begin(); itSn1 != Sn1.end(); ++itSn1)
    {
        long ix = itSn1->SFLSNodeComponent1;
        long iy = itSn1->SFLSNodeComponent2;
        long iz = itSn1->SFLSNodeComponent3;

        LevelSetPixelType* mp_phi_at_idx_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        LevelSetPixelType mp_phi_at_idx = mp_phi_at_idx_ptr[0];

        m_ln1.push_back(*itSn1);

        *(static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz))) = -1;
        *(static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz))) = 1;


        if ( ix+1 < m_nx )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[m_increment0], -3.0))
            {
                Sn2.push_back(NodeType(ix+1, iy, iz));
                mp_phi_at_idx_ptr[m_increment0] = mp_phi_at_idx - 1.0;
            }
        }

        if ( ix-1 >= 0 )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[-m_increment0], -3.0))
            {
                Sn2.push_back(NodeType(ix-1, iy, iz));
                mp_phi_at_idx_ptr[-m_increment0] = mp_phi_at_idx - 1.0;
            }
        }


        if ( iy+1 < m_ny )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[m_increment1], -3.0))
            {
                Sn2.push_back(NodeType(ix, iy+1, iz));
                mp_phi_at_idx_ptr[m_increment1] = mp_phi_at_idx - 1.0;
            }
        }


        if ( iy-1 >= 0 )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[-m_increment1], -3.0))
            {
                Sn2.push_back(NodeType(ix, iy-1, iz));
                mp_phi_at_idx_ptr[-m_increment1] = mp_phi_at_idx - 1.0;
            }
        }

        if ( iz+1 < m_nz )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[m_increment2], -3.0))
            {
                Sn2.push_back(NodeType(ix, iy, iz+1));
                mp_phi_at_idx_ptr[m_increment2] = mp_phi_at_idx - 1.0;
            }
        }


        if ( iz-1 >= 0 )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[-m_increment2], -3.0))
            {
                Sn2.push_back(NodeType(ix, iy, iz-1));
                mp_phi_at_idx_ptr[-m_increment2] = mp_phi_at_idx - 1.0;
            }
        }
    }


    /*--------------------------------------------------
    3.3 Scan Sp1     */
    for (CSFLSLayer::iterator itSp1 = Sp1.begin(); itSp1 != Sp1.end(); ++itSp1)
    {
        long ix = itSp1->SFLSNodeComponent1;
        long iy = itSp1->SFLSNodeComponent2;
        long iz = itSp1->SFLSNodeComponent3;

        LevelSetPixelType* mp_phi_at_idx_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        LevelSetPixelType mp_phi_at_idx = mp_phi_at_idx_ptr[0];

        m_lp1.push_back(*itSp1);

        *(static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz))) = 1;
        *(static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz))) = 0;


        if ( ix+1 < m_nx )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[m_increment0], 3.0))
            {
                Sp2.push_back(NodeType(ix+1, iy, iz));
                mp_phi_at_idx_ptr[m_increment0] = mp_phi_at_idx + 1.0;
            }
        }

        if ( ix-1 >= 0 )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[-m_increment0], 3.0))
            {
                Sp2.push_back(NodeType(ix-1, iy, iz));
                mp_phi_at_idx_ptr[-m_increment0] = mp_phi_at_idx + 1.0;
            }
        }


        if ( iy+1 < m_ny )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[m_increment1], 3.0))
            {
                Sp2.push_back(NodeType(ix, iy+1, iz));
                mp_phi_at_idx_ptr[m_increment1] = mp_phi_at_idx + 1.0;
            }
        }

        if ( iy-1 >= 0 )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[-m_increment1], 3.0))
            {
                Sp2.push_back(NodeType(ix, iy-1, iz));
                mp_phi_at_idx_ptr[-m_increment1] = mp_phi_at_idx + 1.0;
            }
        }

        if ( iz+1 < m_nz )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[m_increment2], 3.0))
            {
                Sp2.push_back(NodeType(ix, iy, iz+1));
                mp_phi_at_idx_ptr[m_increment2] = mp_phi_at_idx + 1.0;
            }
        }


        if ( iz-1 >= 0 )
        {
            if (floatingEqual<LevelSetPixelType>(mp_phi_at_idx_ptr[-m_increment2], 3.0))
            {
                Sp2.push_back(NodeType(ix, iy, iz-1));
                mp_phi_at_idx_ptr[-m_increment2] = mp_phi_at_idx + 1.0;
            }
        }
    }


    /*--------------------------------------------------
    3.4 Scan Sn2     */
    {
        //debug
        int aaa = 0;
        for (CSFLSLayer::iterator itSn2 = Sn2.begin(); itSn2 != Sn2.end(); ++itSn2, ++aaa)
        {
            long ix = itSn2->SFLSNodeComponent1;
            long iy = itSn2->SFLSNodeComponent2;
            long iz = itSn2->SFLSNodeComponent3;

            m_ln2.push_back(*itSn2);

            *(static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz))) = -2;
            *(static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz))) = 1;
        }
    }



    /*--------------------------------------------------
    3.5 Scan Sp2     */
    for (CSFLSLayer::iterator itSp2 = Sp2.begin(); itSp2 != Sp2.end(); ++itSp2)
    {
        long ix = itSp2->SFLSNodeComponent1;
        long iy = itSp2->SFLSNodeComponent2;
        long iz = itSp2->SFLSNodeComponent3;

        m_lp2.push_back(*itSp2);

        *(static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz))) = 2;
        *(static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz))) = 0;
    }

}

/*================================================================================
  initializeLabel*/
void CSFLSSegmentor3D::initializeLabel()
{
    if (m_nx + m_ny + m_nz == 0 || m_increment0 + m_increment1 + m_increment2 == 0)
    {
        std::cerr<<"set mp_img first.\n";
        abort();
    }

    //find interface and mark as 0, create Lz
    LabelImagePixelType defaultLabel = 0;

    int* size = mp_img->GetDimensions();

    mp_label = vtkImageData::New();
    mp_label->SetDimensions(size);
    mp_label->SetOrigin(mp_img->GetOrigin());
    mp_label->SetSpacing(mp_img->GetSpacing());
    mp_label->SetInformation(mp_img->GetInformation());
    //mp_label->SetExtent(mp_img->GetExtent());
#if VTK_MAJOR_VERSION <= 5
    mp_label->SetNumberOfScalarComponents(1);
    mp_label->SetScalarTypeToShort();
    mp_label->AllocateScalars();
#else
    mp_label->AllocateScalars(VTK_SHORT, 1);
#endif

    LabelPixelType* mp_label_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(0, 0, 0));

    mp_label_mask = vtkImageData::New();
    mp_label_mask->SetDimensions(size);
    mp_label_mask->SetOrigin(mp_img->GetOrigin());
    mp_label_mask->SetSpacing(mp_img->GetSpacing());
    mp_label_mask->SetInformation(mp_img->GetInformation());
    //mp_label_mask->SetExtent(mp_img->GetExtent());
#if VTK_MAJOR_VERSION <= 5
    mp_label_mask->SetNumberOfScalarComponents(1);
    mp_label_mask->SetScalarTypeToShort();
    mp_label_mask->AllocateScalars();
#else
    mp_label_mask->AllocateScalars(VTK_SHORT, 1);
#endif

    LabelImagePixelType* mp_label_mask_pixel_ptr = static_cast<LabelImagePixelType*>(mp_label_mask->GetScalarPointer(0, 0, 0));


    long n = size[0]*size[1]*size[2];
    for (long i = 0; i < n; i += m_increment0)
    {
        mp_label_ptr[i] = defaultLabel;
        mp_label_mask_pixel_ptr[i] = defaultLabel;
    }

    std::cout<<"mp_label_mask is initied\n"<<std::flush;


    // dbg
    vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
    writer->SetFileName("/tmp/mp_label_mask-1.mhd");
    writer->SetRAWFileName("/tmp/mp_label_mask-1.raw");
    writer->SetInput(mp_label_mask);
    writer->Write();
    // dbg, end



    //    for (int z = 0; z < size[2]; z++)
    //    {
    //        for (int y = 0; y < size[1]; y++)
    //        {
    //            for (int x = 0; x < size[0]; x++)
    //            {
    //                *(static_cast<LabelPixelType*>(mp_label->GetScalarPointer(x,y,z))) = defaultLabel;
    //            }
    //        }
    //    }


    return;
}


/*================================================================================
  initializePhi*/
void CSFLSSegmentor3D::initializePhi()
{
    if (m_nx + m_ny + m_nz == 0)
    {
        std::cerr<<"set mp_img first.\n";
        abort();
    }

    LevelSetPixelType arbitraryInitPhi = 1000.0;

    int* size = mp_img->GetDimensions();

    mp_phi = vtkImageData::New();

    mp_phi->SetDimensions(size);
    mp_phi->SetOrigin(mp_img->GetOrigin());
    mp_phi->SetSpacing(mp_img->GetSpacing());
    mp_phi->SetInformation(mp_img->GetInformation());

#if VTK_MAJOR_VERSION <= 5
    mp_phi->SetNumberOfScalarComponents(1);
    mp_phi->SetScalarTypeToFloat();
    mp_phi->AllocateScalars();
#else
    mp_phi->AllocateScalars(VTK_FLOAT, 1);
#endif

    int startIdx[] = {0, 0, 0};
    LevelSetPixelType* mp_phi_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(startIdx));
    long n = size[0]*size[1]*size[2];
    for (long i = 0; i < n; i += m_increment0)
    {
        mp_phi_ptr[i] = arbitraryInitPhi;
    }


    //    // dbg
    //    vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
    //    writer->SetFileName("initPhi.mhd");
    //    writer->SetRAWFileName("initPhi.raw");
    //    writer->SetInput(mp_phi);
    //    writer->Write();
    //    // dbg, end


    return;
}


/* ============================================================
   initializeSFLSFromMask    */
void CSFLSSegmentor3D::initializeSFLSFromMask()
{
    if (!mp_mask)
    {
        std::cerr<<"set mp_mask first.\n";
        abort();
    }


    initializePhi();
    initializeLabel();

    MaskPixelType* mp_mask_pixel_ptr = 0;
    LevelSetPixelType* mp_phi_pixel_ptr = 0;
    LabelPixelType* mp_label_pixel_ptr = 0;
    LabelPixelType* mp_label_mask_pixel_ptr = 0;

    {
      // dbg
      char mhdName[1000];
      sprintf(mhdName, "/tmp/mp_label_mask-initSFLS-%d.mhd", 1);
      char rawName[1000];
      sprintf(rawName, "/tmp/mp_label_mask-initSFLS-%d.raw", 1);

      vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
      writer->SetFileName(mhdName);
      writer->SetRAWFileName(rawName);
      writer->SetInput(mp_label_mask);
      writer->Write();
      // dbg, end
    }

    for (long iz = 0; iz < m_nz; ++iz)
    {
        for (long iy = 0; iy < m_ny; ++iy)
        {
            for (long ix = 0; ix < m_nx; ++ix)
            {
                //mark the inside and outside of label and phi
                mp_mask_pixel_ptr = static_cast<MaskPixelType*>(mp_mask->GetScalarPointer(ix, iy, iz));
                mp_phi_pixel_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
                mp_label_pixel_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));
                mp_label_mask_pixel_ptr = static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz));

                if( mp_mask_pixel_ptr[0] == 0 )
                {
                    mp_label_pixel_ptr[0] = 3;
                    mp_phi_pixel_ptr[0] = 3.0;
                    mp_label_mask_pixel_ptr[0] = 0;
                }
                else
                {
                    mp_label_pixel_ptr[0] = -3;
                    mp_phi_pixel_ptr[0] = -3.0;
                    mp_label_mask_pixel_ptr[0] = 1;

                    ++m_insideVoxelCount;

                    if ( (iy+1 < m_ny && *(mp_mask_pixel_ptr + m_increment1) == 0)	\
                         || (iy-1 >= 0 && *(mp_mask_pixel_ptr - m_increment1) == 0)	\
                         || (ix+1 < m_nx && *(mp_mask_pixel_ptr + m_increment0) == 0) \
                         || (ix-1 >= 0 && *(mp_mask_pixel_ptr - m_increment0) == 0)	\
                         || (iz+1 < m_nz && *(mp_mask_pixel_ptr + m_increment2) == 0) \
                         || (iz-1 >= 0 && *(mp_mask_pixel_ptr - m_increment2) == 0) )
                    {
                        m_lz.push_back(NodeType(ix, iy, iz));
                        mp_label_pixel_ptr[0] = 0;
                        mp_phi_pixel_ptr[0] = 0.0;

                        mp_label_mask_pixel_ptr[0] = 1;
                    }
                }
            }
        }
    }



    {
      // dbg
      char mhdName[1000];
      sprintf(mhdName, "/tmp/mp_label_mask-initSFLS-%d.mhd", 2);
      char rawName[1000];
      sprintf(rawName, "/tmp/mp_label_mask-initSFLS-%d.raw", 2);

      vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
      writer->SetFileName(mhdName);
      writer->SetRAWFileName(rawName);
      writer->SetInput(mp_label_mask);
      writer->Write();
      // dbg, end
    }

    m_insideVolume = m_insideVoxelCount*m_dx*m_dy*m_dz;

    //scan Lz to create Ln1 and Lp1
    for (CSFLSLayer::const_iterator it = m_lz.begin(); it != m_lz.end(); ++it)
    {
        long ix = it->SFLSNodeComponent1;
        long iy = it->SFLSNodeComponent2;
        long iz = it->SFLSNodeComponent3;

        mp_phi_pixel_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        mp_label_pixel_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));
        mp_label_mask_pixel_ptr = static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz));

        if(ix+1 < m_nx)
        {
            int idx[] = {ix+1, iy, iz};

            if ( mp_label_pixel_ptr[m_increment0] == 3 )
            {
                mp_label_pixel_ptr[m_increment0] = 1;
                mp_phi_pixel_ptr[m_increment0] = 1.0;
                mp_label_mask_pixel_ptr[m_increment0] = 0;

                m_lp1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
            else if ( mp_label_pixel_ptr[m_increment0] == -3 )
            {
                mp_label_pixel_ptr[m_increment0] = -1;
                mp_phi_pixel_ptr[m_increment0] = -1.0;

                mp_label_mask_pixel_ptr[m_increment0] = 1;

                m_ln1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
        }

        if(ix-1 >= 0)
        {
            int idx[] = {ix-1, iy, iz};

            if ( mp_label_pixel_ptr[-m_increment0] == 3 )
            {
                mp_label_pixel_ptr[-m_increment0] = 1;
                mp_phi_pixel_ptr[-m_increment0] = 1.0;
                mp_label_mask_pixel_ptr[-m_increment0] = 0;

                m_lp1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
            else if ( mp_label_pixel_ptr[-m_increment0] == -3 )
            {
                mp_label_pixel_ptr[-m_increment0] = -1;
                mp_phi_pixel_ptr[-m_increment0] = -1.0;
                mp_label_mask_pixel_ptr[-m_increment0] = 1;

                m_ln1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
        }

        if(iy+1 < m_ny)
        {
            int idx[] = {ix, iy+1, iz};


            if ( mp_label_pixel_ptr[m_increment1] == 3 )
            {
                mp_label_pixel_ptr[m_increment1] = 1;
                mp_phi_pixel_ptr[m_increment1] = 1.0;
                mp_label_mask_pixel_ptr[m_increment1] = 0;

                m_lp1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
            else if (mp_label_pixel_ptr[m_increment1] == -3 )
            {
                mp_label_pixel_ptr[m_increment1] = -1;
                mp_phi_pixel_ptr[m_increment1] = -1.0;
                mp_label_mask_pixel_ptr[m_increment1] = 1;

                m_ln1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
        }

        if(iy-1 >= 0)
        {
            int idx[] = {ix, iy-1, iz};

            if ( mp_label_pixel_ptr[-m_increment1] == 3 )
            {
                mp_label_pixel_ptr[-m_increment1] = 1;
                mp_phi_pixel_ptr[-m_increment1] = 1.0;
                mp_label_mask_pixel_ptr[-m_increment1] = 0;

                m_lp1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
            else if (mp_label_pixel_ptr[-m_increment1] == -3 )
            {
                mp_label_pixel_ptr[-m_increment1] = -1;
                mp_phi_pixel_ptr[-m_increment1] = -1.0;
                mp_label_mask_pixel_ptr[-m_increment1] = 1;

                m_ln1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
        }

        if(iz+1 < m_nz)
        {
            int idx[] = {ix, iy, iz+1};

            if ( mp_label_pixel_ptr[m_increment2] == 3 )
            {
                mp_label_pixel_ptr[m_increment2] = 1;
                mp_phi_pixel_ptr[m_increment2]= 1.0;
                mp_label_mask_pixel_ptr[m_increment2] = 0;

                m_lp1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
            else if (mp_label_pixel_ptr[m_increment2] == -3 )
            {
                mp_label_pixel_ptr[m_increment2] = -1;
                mp_phi_pixel_ptr[m_increment2] = -1.0;
                mp_label_mask_pixel_ptr[m_increment2] = 1;

                m_ln1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
        }

        if(iz-1 >= 0)
        {
            int idx[] = {ix, iy, iz-1};

            if ( mp_label_pixel_ptr[-m_increment2] == 3 )
            {
                mp_label_pixel_ptr[-m_increment2] = 1;
                mp_phi_pixel_ptr[-m_increment2] = 1.0;
                mp_label_mask_pixel_ptr[-m_increment2] = 0;

                m_lp1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
            else if (mp_label_pixel_ptr[-m_increment2] == -3 )
            {
                mp_label_pixel_ptr[-m_increment2] = -1;
                mp_phi_pixel_ptr[-m_increment2] = -1.0;
                mp_label_mask_pixel_ptr[-m_increment2] = 1;

                m_ln1.push_back( NodeType(idx[0], idx[1], idx[2]) );
            }
        }
    }


    //scan Ln1 to create Ln2
    for (CSFLSLayer::const_iterator it = m_ln1.begin(); it != m_ln1.end(); ++it)
    {
        long ix = it->SFLSNodeComponent1;
        long iy = it->SFLSNodeComponent2;
        long iz = it->SFLSNodeComponent3;

        mp_phi_pixel_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        mp_label_pixel_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));
        mp_label_mask_pixel_ptr = static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz));

        if(ix+1 < m_nx && mp_label_pixel_ptr[m_increment0] == -3 )
        {
            mp_label_pixel_ptr[m_increment0] = -2;
            mp_phi_pixel_ptr[m_increment0] = -2.0;
            mp_label_mask_pixel_ptr[m_increment0] = 1;

            m_ln2.push_back( NodeType(ix+1, iy, iz) );
        }

        if(ix-1 >= 0 && mp_label_pixel_ptr[-m_increment0] == -3 )
        {
            mp_label_pixel_ptr[-m_increment0] = -2;
            mp_phi_pixel_ptr[-m_increment0] = -2.0;
            mp_label_mask_pixel_ptr[-m_increment0] = 1;

            m_ln2.push_back( NodeType(ix-1, iy, iz) );
        }

        if(iy+1 < m_ny && mp_label_pixel_ptr[m_increment1] == -3 )
        {
            mp_label_pixel_ptr[m_increment1] = -2;
            mp_phi_pixel_ptr[m_increment1] = -2.0;
            mp_label_mask_pixel_ptr[m_increment1] = 1;

            m_ln2.push_back( NodeType(ix, iy+1, iz) );
        }

        if(iy-1 >= 0 && mp_label_pixel_ptr[-m_increment1] == -3 )
        {
            mp_label_pixel_ptr[-m_increment1] = -2;
            mp_phi_pixel_ptr[-m_increment1] = -2.0;
            mp_label_mask_pixel_ptr[-m_increment1] = 1;

            m_ln2.push_back( NodeType(ix, iy-1, iz) );
        }

        if(iz+1 < m_nz && mp_label_pixel_ptr[m_increment2] == -3 )
        {
            mp_label_pixel_ptr[m_increment2] = -2;
            mp_phi_pixel_ptr[m_increment2] = -2.0;
            mp_label_mask_pixel_ptr[m_increment2] = 1;

            m_ln2.push_back( NodeType(ix, iy, iz+1) );
        }

        if(iz-1 >= 0 && mp_label_pixel_ptr[-m_increment2] == -3 )
        {
            mp_label_pixel_ptr[-m_increment2] = -2;
            mp_phi_pixel_ptr[-m_increment2] = -2.0;
            mp_label_mask_pixel_ptr[-m_increment2] = 1;

            m_ln2.push_back( NodeType(ix, iy, iz-1) );
        }
    }

    {
      // dbg
      char mhdName[1000];
      sprintf(mhdName, "/tmp/mp_label_mask-initSFLS-%d.mhd", 3);
      char rawName[1000];
      sprintf(rawName, "/tmp/mp_label_mask-initSFLS-%d.raw", 3);

      vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
      writer->SetFileName(mhdName);
      writer->SetRAWFileName(rawName);
      writer->SetInput(mp_label_mask);
      writer->Write();
      // dbg, end
    }

    //scan Lp1 to create Lp2
    for (CSFLSLayer::const_iterator it = m_lp1.begin(); it != m_lp1.end(); ++it)
    {
        long ix = it->SFLSNodeComponent1;
        long iy = it->SFLSNodeComponent2;
        long iz = it->SFLSNodeComponent3;

        mp_phi_pixel_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(ix, iy, iz));
        mp_label_pixel_ptr = static_cast<LabelPixelType*>(mp_label->GetScalarPointer(ix, iy, iz));
        mp_label_mask_pixel_ptr = static_cast<LabelPixelType*>(mp_label_mask->GetScalarPointer(ix, iy, iz));

        if(ix+1 < m_nx && mp_label_pixel_ptr[m_increment0] == 3 )
        {
            mp_label_pixel_ptr[m_increment0] = 2;
            mp_phi_pixel_ptr[m_increment0] = 2.0;
            mp_label_mask_pixel_ptr[m_increment0] = 0;

            m_lp2.push_back( NodeType(ix+1, iy, iz) );
        }

        if(ix-1 >= 0 && mp_label_pixel_ptr[-m_increment0] == 3 )
        {
            mp_label_pixel_ptr[-m_increment0] = 2;
            mp_phi_pixel_ptr[-m_increment0] = 2.0;
            mp_label_mask_pixel_ptr[-m_increment0] = 0;

            m_lp2.push_back( NodeType(ix-1, iy, iz) );
        }

        if(iy+1 < m_ny && mp_label_pixel_ptr[m_increment1] == 3 )
        {
            mp_label_pixel_ptr[m_increment1] = 2;
            mp_phi_pixel_ptr[m_increment1] = 2.0;
            mp_label_mask_pixel_ptr[m_increment1] = 0;

            m_lp2.push_back( NodeType(ix, iy+1, iz) );
        }

        if(iy-1 >= 0 && mp_label_pixel_ptr[-m_increment1] == 3 )
        {
            mp_label_pixel_ptr[-m_increment1] = 2;
            mp_phi_pixel_ptr[-m_increment1] = 2.0;
            mp_label_mask_pixel_ptr[-m_increment1] = 0;

            m_lp2.push_back( NodeType(ix, iy-1, iz) );
        }

        if(iz+1 < m_nz && mp_label_pixel_ptr[m_increment2] == 3 )
        {
            mp_label_pixel_ptr[m_increment2] = 2;
            mp_phi_pixel_ptr[m_increment2] = 2.0;
            mp_label_mask_pixel_ptr[m_increment2] = 0;

            m_lp2.push_back( NodeType(ix, iy, iz+1) );
        }

        if(iz-1 >= 0 && mp_label_pixel_ptr[-m_increment2] == 3 )
        {
            mp_label_pixel_ptr[-m_increment2] = 2;
            mp_phi_pixel_ptr[-m_increment2] = 2.0;
            mp_label_mask_pixel_ptr[-m_increment2] = 0;

            m_lp2.push_back( NodeType(ix, iy, iz-1) );
        }
    }

    {
      // dbg
      char mhdName[1000];
      sprintf(mhdName, "/tmp/mp_label_mask-initSFLS-%d.mhd", 4);
      char rawName[1000];
      sprintf(rawName, "/tmp/mp_label_mask-initSFLS-%d.raw", 4);

      vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
      writer->SetFileName(mhdName);
      writer->SetRAWFileName(rawName);
      writer->SetInput(mp_label_mask);
      writer->Write();
      // dbg, end
    }
}

///* getLevelSetFunction */
//vtkImageData* CSFLSSegmentor3D::getLevelSetFunction()
//{
//    //   if (!m_done)
//    //     {
//    //       std::cerr<<"Error: not done.\n";
//    //       raise(SIGABRT);
//    //     }

//    return mp_phi;
//}


/*============================================================
   computeKappa
   Compute kappa at a point in the zero level set           */
double CSFLSSegmentor3D::computeKappa(long ix, long iy, long iz)
{
    //double kappa = 0;

    double dx = 0;
    double dy = 0;
    double dz = 0;

    double dxx = 0;
    double dyy = 0;
    double dzz = 0;

    double dx2 = 0;
    double dy2 = 0;
    double dz2 = 0;

    double dxy = 0;
    double dxz = 0;
    double dyz = 0;

    char xok = 0;
    char yok = 0;
    char zok = 0;

    int idx[] = {ix, iy, iz};


    if( ix+1 < m_nx && ix-1 >=0 )
    {
        xok = 1;
    }

    if( iy+1 < m_ny && iy-1 >=0 )
    {
        yok = 1;
    }

    if( iz+1 < m_nz && iz-1 >=0 )
    {
        zok = 1;
    }

    LevelSetPixelType* mp_phi_idx_ptr = static_cast<LevelSetPixelType*>(mp_phi->GetScalarPointer(idx));
    LevelSetPixelType tmp = mp_phi_idx_ptr[0];

    //assert(!isnan(tmp));

    if (xok)
    {
        LevelSetPixelType tmp1 = mp_phi_idx_ptr[-m_increment0];
        LevelSetPixelType tmp2 = mp_phi_idx_ptr[m_increment0];

        dx  = (tmp2 - tmp1 )/(2.0*m_dx);
        dxx = (tmp2 - 2.0*tmp + tmp1)/(m_dx*m_dx);
        dx2 = dx*dx;
    }

    if (yok)
    {
        LevelSetPixelType tmp3 = mp_phi_idx_ptr[-m_increment1];
        LevelSetPixelType tmp4 = mp_phi_idx_ptr[m_increment1];

        dy  = ((tmp4 - tmp3 ))/(2.0*m_dy);
        dyy = (tmp4 - 2*tmp + tmp3)/(m_dy*m_dy);
        dy2 = dy*dy;
    }

    if (zok)
    {
        LevelSetPixelType tmp5 = mp_phi_idx_ptr[-m_increment2];
        LevelSetPixelType tmp6 = mp_phi_idx_ptr[m_increment2];

        dz  = ((tmp6 - tmp5 ))/(2.0*m_dz);
        dzz = (tmp6 - 2.0*tmp + tmp5)/(m_dz*m_dz);
        dz2 = dz*dz;
    }


    if(xok && yok)
    {
        dxy = 0.25*(mp_phi_idx_ptr[m_increment0 + m_increment1]\
                    + mp_phi_idx_ptr[-m_increment0 - m_increment1]\
                    - mp_phi_idx_ptr[m_increment0 - m_increment1]\
                    - mp_phi_idx_ptr[-m_increment0 + m_increment1])/(m_dx*m_dy);
    }

    if(xok && zok)
    {
        dxy = 0.25*(mp_phi_idx_ptr[m_increment0 + m_increment2]\
                    + mp_phi_idx_ptr[-m_increment0 - m_increment2]\
                    - mp_phi_idx_ptr[m_increment0 - m_increment2]\
                    - mp_phi_idx_ptr[-m_increment0 + m_increment2])/(m_dx*m_dz);
    }

    if(yok && zok)
    {
        dxy = 0.25*(mp_phi_idx_ptr[m_increment1 + m_increment2]\
                    + mp_phi_idx_ptr[-m_increment1 - m_increment2]\
                    - mp_phi_idx_ptr[m_increment1 - m_increment2]\
                    - mp_phi_idx_ptr[-m_increment1 + m_increment2])/(m_dy*m_dz);
    }

    double k = (dxx*(dy2 + dz2) + dyy*(dx2 + dz2) + dzz*(dx2 + dy2) - 2*dx*dy*dxy - 2*dx*dz*dxz - 2*dy*dz*dyz)/(dx2 + dy2 + dz2 + m_eps);

    //assert(!isnan(k));

    return k;
}



#endif
