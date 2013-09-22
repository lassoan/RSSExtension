#ifndef SFLSRobustStatSegmentor3DLabelMap_txx_
#define SFLSRobustStatSegmentor3DLabelMap_txx_

#include "SFLSRobustStatSegmentor3DLabelMap.h"

#include <algorithm>
#include <ctime>
#include <cstdlib>

#include <limits>

// //debug//
#include <fstream>
// //DEBUG//

#include "vtkImageData.h"

#include "vtkImageCast.h"

// dbg
#include "vtkMetaImageWriter.h"


/* ============================================================   */

void
CSFLSRobustStatSegmentor3DLabelMap::basicInit()
{
    SuperClassType::basicInit();

    m_statNeighborX = 1;
    m_statNeighborY = 1;
    m_statNeighborZ = 1;

    m_kernelWidthFactor = 10.0;

    m_inputImageIntensityMin = 0;
    m_inputImageIntensityMax = 0;

    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::setInputLabelImage(vtkImageData* l)
{
    vtkImageCast* castFilter = vtkImageCast::New();
    castFilter->SetInput(l);
    castFilter->SetOutputScalarTypeToShort();
    castFilter->Update();

    m_inputLabelImage = castFilter->GetOutput();
    m_inputLabelImage_buffer_ptr = static_cast<LabelImagePixelType*>(m_inputLabelImage->GetScalarPointer(0, 0, 0));

    int* size = m_inputLabelImage->GetDimensions();

    if (this->m_nx + this->m_ny + this->m_nz == 0)
    {
        this->m_nx = size[0];
        this->m_ny = size[1];
        this->m_nz = size[2];
    }
    else if ( this->m_nx != (long)size[0] || this->m_ny != (long)size[1] || this->m_nz != (long)size[2] )
    {
        std::cerr<<"Error: image sizes do not match with label image size.\n";
        abort();;
    }


//    //dbg
//    vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
//    writer->SetFileName("/tmp/m_inputLabelImage.mhd");
//    writer->SetRAWFileName("/tmp/m_inputLabelImage.raw");
//    writer->SetInput(m_inputLabelImage);
//    writer->Write();
//    //dbg, end

    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::computeForce()
{
//    int np = mp_img->GetNumberOfPoints();
//    std::cout<<"0 np = "<<np<<std::endl<<std::flush;

    double fmax = std::numeric_limits<double>::min();
    double kappaMax = std::numeric_limits<double>::min();

    long n = this->m_lz.size();
    double* kappaOnZeroLS = new double[ n ];
    double* cvForce = new double[ n ];


    std::vector<CSFLSLayer::iterator> m_lzIterVct( n );
    {
        long iiizzz = 0;
        for (CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz)
            m_lzIterVct[iiizzz++] = itz;
    }

//    std::cout<<"77777777777777777\n"<<std::flush;

    for (long i = 0; i < n; ++i)
    {
//        std::cout<<i<<std::endl<<std::flush;

        CSFLSLayer::iterator itz = m_lzIterVct[i];

        long ix = itz->SFLSNodeComponent1;
        long iy = itz->SFLSNodeComponent2;
        long iz = itz->SFLSNodeComponent3;

//        std::cout<<ix<<'\t'<<iy<<'\t'<<iz<<std::endl<<std::flush;

        int idx[] = {ix, iy, iz};

        kappaOnZeroLS[i] = this->computeKappa(ix, iy, iz);

//        std::cout<<"kappaOnZeroLS[i] = "<<kappaOnZeroLS[i]<<std::endl<<std::flush;

        std::vector<FeatureImagePixelType> f(m_numberOfFeature);

        computeFeatureAt(idx, f);

//        std::cout<<"feature = "<<f[0]<<'\t'<<f[1]<<'\t'<<f[2]<<std::endl<<std::flush;

        double a = -kernelEvaluationUsingPDF(f);

//        std::cout<<"a = "<<a<<std::endl<<std::flush;


        fmax = fmax>fabs(a)?fmax:fabs(a);
        kappaMax = kappaMax>fabs(kappaOnZeroLS[i])?kappaMax:fabs(kappaOnZeroLS[i]);

        cvForce[i] = a;
    }


    this->m_force.resize(n);
    for (long i = 0; i < n; ++i)
    {
        //this->m_force.push_back(cvForce[i]/(fmax + 1e-10) +  (this->m_curvatureWeight)*kappaOnZeroLS[i]);
        this->m_force[i] = (1 - (this->m_curvatureWeight))*cvForce[i]/(fmax + 1e-10) \
                +  (this->m_curvatureWeight)*kappaOnZeroLS[i]/(kappaMax + 1e-10);
    }

    
    delete[] kappaOnZeroLS;
    delete[] cvForce;
}

/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::inputLableImageToSeeds()
{
    std::vector<long> thisSeed(3);

    int* size = m_inputLabelImage->GetDimensions();
    for (int iz = 0; iz < size[2]; ++iz)
    {
        for (int iy = 0; iy < size[1]; ++iy)
        {
            for (int ix = 0; ix < size[0]; ++ix)
            {
                if (*(static_cast<LabelImagePixelType*>(m_inputLabelImage->GetScalarPointer(ix, iy, iz))) != 0)
                {
                    thisSeed[0] = ix;
                    thisSeed[1] = iy;
                    thisSeed[2] = iz;
                    m_seeds.push_back(thisSeed);

                    //                    sf<<thisSeed[0]<<", "<<thisSeed[1]<<", "<<thisSeed[2]<<std::endl;
                }
            }
        }
    }

    return;
}

/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::getThingsReady()
{
    /*
   1. Generate mp_mask from seeds
   2. Compute feature at each point
   3. Extract feature at/around the seeds
  */
    inputLableImageToSeeds();

    seedToMask();

    //dialteSeeds();

    initFeatureComputedImage();
    initFeatureImage();


    //computeFeature();
    getFeatureAroundSeeds();
    estimateFeatureStdDevs();

    estimatePDFs();

    return;
}


/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap::initFeatureImage()
{
    if (!(this->mp_img))
    {
        std::cerr<<"Error: set input image first.\n";
        abort();;
    }

    for (long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature)
    {
        vtkImageData* fimg = vtkImageData::New();

        fimg->SetDimensions(mp_img->GetDimensions());
        fimg->SetOrigin(mp_img->GetOrigin());
        fimg->SetSpacing(mp_img->GetSpacing());
        fimg->SetInformation(mp_img->GetInformation());
#if VTK_MAJOR_VERSION <= 5
        //    std::cout<<"VTK_MAJOR_VERSION <= 5\n";
        fimg->SetNumberOfScalarComponents(1);
        fimg->SetScalarTypeToFloat();
        fimg->AllocateScalars();
#else
        fimg->AllocateScalars(VTK_FLOAT, 1);
#endif

        m_featureImageList.push_back(fimg);
    }

    return;
}


/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap::initFeatureComputedImage()
{
    if (!(this->mp_img))
    {
        std::cerr<<"Error: set input image first.\n";
        abort();;
    }

    m_featureComputed = vtkImageData::New();
    int* size = mp_img->GetDimensions();
    m_featureComputed->SetDimensions(size);
    m_featureComputed->SetOrigin(mp_img->GetOrigin());
    m_featureComputed->SetSpacing(mp_img->GetSpacing());
    m_featureComputed->SetInformation(mp_img->GetInformation());
#if VTK_MAJOR_VERSION <= 5
    //    std::cout<<"VTK_MAJOR_VERSION <= 5\n";
    m_featureComputed->SetNumberOfScalarComponents(1);
    m_featureComputed->SetScalarTypeToUnsignedChar();
    m_featureComputed->AllocateScalars();
#else
    m_featureComputed->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
#endif

    int startIdx[] = {0, 0, 0};
    m_featureComputed_pixel_type* m_featureComputed_ptr = static_cast<m_featureComputed_pixel_type*>(m_featureComputed->GetScalarPointer(startIdx));
    long n = size[0]*size[1]*size[2];
    for (long i = 0; i < n; i += m_increment0)
    {
        m_featureComputed_ptr[i] = 0;
    }

    return;
}


/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap::computeFeatureAt(int idx[3], std::vector<FeatureImagePixelType>& f)
{
    f.resize(m_numberOfFeature);
    m_featureComputed_pixel_type* m_featureComputed_idx_ptr = static_cast<m_featureComputed_pixel_type*>(m_featureComputed->GetScalarPointer(idx));

//    std::cout<<"m_featureComputed_idx_ptr[0] = "<<m_featureComputed_idx_ptr[0]<<std::endl<<std::flush;

    if (m_featureComputed_idx_ptr[0])
    {
//        std::cout<<"alive here\n"<<std::flush;

        // the feature at this pixel is computed, just retrive
        for (long i = 0; i < m_numberOfFeature; ++i)
        {
            f[i] = *(static_cast<FeatureImagePixelType*>((m_featureImageList[i])->GetScalarPointer(idx)));
        }
    }
    else
    {
//        std::cout<<"alive else 0 \n"<<std::flush;

        // compute the feature
        std::vector< TPixel > neighborIntensities;

        long ix = idx[0];
        long iy = idx[1];
        long iz = idx[2];

//        std::cout<<ix<<'\t'<<iy<<'\t'<<iz<<std::endl<<std::flush;

//        std::cout<<"mp_img = "<<mp_img<<std::endl<<std::flush;
//        int np = mp_img->GetNumberOfPoints();
//        std::cout<<"np = "<<np<<std::endl<<std::flush;

//        std::cout<<mp_img->Get

//        TPixel* img_ptr = static_cast<TPixel*>(this->mp_img->GetScalarPointer(ix, iy, iz));
        TPixel* img_ptr = mp_img_ptr + iz*m_nx*m_ny + iy*m_nx + ix;

//        std::cout<<mp_img_ptr<<std::endl<<std::flush;
//        std::cout<<img_ptr<<std::endl<<std::flush;

//        std::cout<<"alive else 1 \n"<<std::flush;



        for (long iiz = - m_statNeighborZ; iiz <= m_statNeighborZ; ++iiz)
        {
            for (long iiy = - m_statNeighborY; iiy <= m_statNeighborY; ++iiy)
            {
                for (long iix = - m_statNeighborX; iix <= m_statNeighborX; ++iix)
                {
                    if (0 <= ix-iix && ix+iix < this->m_nx && 0 <= iy-iiy && iy+iiy < this->m_ny && 0 <= iz-iiz && iz+iiz < this->m_nz)
                    {
                        TPixel tmp = img_ptr[iiz*m_increment2 + iiy*m_increment1 + iix*m_increment0];
                        //                        TPixel tmp = *(static_cast<TPixel*>(mp_img->GetScalarPointer(iix, iiy, iiz)));
                        neighborIntensities.push_back(tmp);
                    }
                }
            }
        }

//        std::cout<<"alive else\n"<<std::flush;

        getRobustStatistics(neighborIntensities, f);

        for (long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature)
        {
            *(static_cast<FeatureImagePixelType*>((m_featureImageList[ifeature])->GetScalarPointer(idx))) = static_cast<FeatureImagePixelType>(f[ifeature]);
        }

        m_featureComputed_idx_ptr[0] = 1; // mark as computed
    }

    return;
}


void
CSFLSRobustStatSegmentor3DLabelMap::doSegmenationBeforeIteration()
{

    getThingsReady();

    //    std::ofstream f("/tmp/d.txt", std::ios_base::app);
    //    f<<"m_maxRunningTime = "<<this->m_maxRunningTime<<std::endl;
    //    f.close();



    /*============================================================
     * From the initial mask, generate: 1. SFLS, 2. mp_label and
     * 3. mp_phi.
     */
    this->initializeSFLS();
}

// void
// CSFLSRobustStatSegmentor3DLabelMap::doSegmenationIteration()
// {
//     //douher::saveAsImage2< double >(mp_phi, "initPhi.nrrd");
//     while(!m_done)        //for (unsigned int it = 0; ; ++it)
//     {
//         inOneSegmentationIteration();
//     }
// }



void CSFLSRobustStatSegmentor3DLabelMap::inOneSegmentationIteration()
{
    if (m_done)
    {
        return;
    }

//    //dbg//
//    std::cout<<"In iteration "<<m_currentIteration<<std::endl<<std::flush;
//    //DBG//


//    // dbg
//    char mhdName[1000];
//    sprintf(mhdName, "/tmp/mp_label_mask-before-iteration-%d.mhd", m_currentIteration);
//    char rawName[1000];
//    sprintf(rawName, "/tmp/mp_label_mask-before-iteration-%d.raw", m_currentIteration);

//    vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
//    writer->SetFileName(mhdName);
//    writer->SetRAWFileName(rawName);
//    writer->SetInput(mp_label_mask);
//    writer->Write();
//    // dbg, end


    float oldVoxelCount = this->m_insideVoxelCount;

    computeForce();

//    std::cout<<"444444444444444444444444444444\n"<<std::flush;

    this->normalizeForce();

//    std::cout<<"5555555555555555555"<<std::flush;

    this->oneStepLevelSetEvolution();

//    std::cout<<"66666666666666666666"<<std::flush;


    /*----------------------------------------------------------------------
      If the level set stops growing, stop */
    this->updateInsideVoxelCount();
    if (m_currentIteration > 2 && oldVoxelCount >= this->m_insideVoxelCount)
    {
        m_done = true;
    }
    /* If the level set stops growing, stop
       ----------------------------------------------------------------------*/


    /*----------------------------------------------------------------------
      If the inside physical volume exceed expected volume, stop */
    float volumeIn = (this->m_insideVoxelCount)*(this->m_dx)*(this->m_dy)*(this->m_dz);
    if (volumeIn > (this->m_maxVolume))
    {
        m_done = true;
    }
    /*If the inside physical volume exceed expected volume, stop
      ----------------------------------------------------------------------*/

    ++m_currentIteration;

    if (m_currentIteration >= m_numIter)
    {
        m_done = true;
    }

    mp_label_mask->Modified();
}

/* ============================================================  */

// void
// CSFLSRobustStatSegmentor3DLabelMap::doSegmenation()
// {
//     double startingTime = clock();

//     getThingsReady();

//     //    std::ofstream f("/tmp/d.txt", std::ios_base::app);
//     //    f<<"m_maxRunningTime = "<<this->m_maxRunningTime<<std::endl;
//     //    f.close();



//     /*============================================================
//    * From the initial mask, generate: 1. SFLS, 2. mp_label and
//    * 3. mp_phi.
//    */
//     this->initializeSFLS();

//     for (unsigned int it = 0; it < this->m_numIter; ++it)
//         //for (unsigned int it = 0; ; ++it)
//     {
//         //dbg//
//         std::cout<<"In iteration "<<it<<std::endl<<std::flush;
//         //DBG//

//         // keep current zero contour as history is required
//         if (this->m_keepZeroLayerHistory)
//         {
//             (this->m_zeroLayerHistory).push_back(this->m_lz);
//         }



//         double oldVoxelCount = this->m_insideVoxelCount;

//         computeForce();

//         this->normalizeForce();

//         this->oneStepLevelSetEvolution();



//         /*----------------------------------------------------------------------
//         If the level set stops growing, stop */
//         this->updateInsideVoxelCount();
//         if (it > 2 && oldVoxelCount >= this->m_insideVoxelCount)
//         {
//             std::ofstream f("/tmp/o.txt");
//             f<<"stop grow\n";
//             f.close();

//             break;
//         }

//         /* If the level set stops growing, stop
//          ----------------------------------------------------------------------*/


//         /*----------------------------------------------------------------------
//         If the inside physical volume exceed expected volume, stop */
//         double volumeIn = (this->m_insideVoxelCount)*(this->m_dx)*(this->m_dy)*(this->m_dz);
//         if (volumeIn > (this->m_maxVolume))
//         {
//             //          std::fstream f("/tmp/o.txt", std::ios_base::app);
//             std::ofstream f("/tmp/o.txt");
//             f<<"m_maxVolume = "<<this->m_maxVolume<<std::endl;
//             f<<"volumeIn = "<<volumeIn<<std::endl;

//             f<<"reach max volume\n";
//             f.close();


//             break;
//         }
//         /*If the inside physical volume exceed expected volume, stop
//         ----------------------------------------------------------------------*/


//         double ellapsedTime = (clock() - startingTime)/static_cast<double>(CLOCKS_PER_SEC);
//         if (ellapsedTime > (this->m_maxRunningTime))
//         {
//             std::ofstream f("/tmp/o.txt");
//             f<<"running time = "<<ellapsedTime<<std::endl;
//             f<<"m_maxRunningTime = "<<this->m_maxRunningTime<<std::endl;
//             f.close();

//             break;
//         }


//     }

//     return;
// }


/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap::getRobustStatistics(std::vector<TPixel>& samples, std::vector<float>& robustStat)
{
    /* note, sample is sorted, so the order is changed */
    robustStat.resize(m_numberOfFeature);

    std::sort(samples.begin(), samples.end() );

    double n = samples.size();

    double q1 = n/4.0;
    double q1_floor;
    double l1 = modf(q1, &q1_floor);

    double q2 = n/2.0;
    double q2_floor;
    double l2 = modf(q2, &q2_floor);

    double q3 = 3.0*n/4.0;
    double q3_floor;
    double l3 = modf(q3, &q3_floor);

    double median = (1 - l2)*samples[static_cast<long>(q2_floor)] + l2*samples[static_cast<long>(q2_floor) + 1];

    double iqr = ( (1 - l3)*samples[static_cast<long>(q3_floor)] + l3*samples[static_cast<long>(q3_floor) + 1] ) \
            - ( (1 - l1)*samples[static_cast<long>(q1_floor)] + l1*samples[static_cast<long>(q1_floor) + 1] );

    robustStat[0] = median;
    robustStat[1] = iqr;

    /* next compute MAD */
    long nn = samples.size();
    std::vector<double> samplesDeMedian(nn);
    for (long i = 0; i < nn; ++i)
    {
        samplesDeMedian[i] = fabs(samples[i] - median);
    }

    std::sort(samplesDeMedian.begin(), samplesDeMedian.end() );

    double mad = (1 - l2)*samplesDeMedian[static_cast<long>(q2_floor)] + l2*samplesDeMedian[static_cast<long>(q2_floor) + 1];
    robustStat[2] = mad;


    return;
}


/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap
::seedToMask()
{
    if (!(this->mp_img))
    {
        std::cerr<<"Error: set input image first.\n";
        abort();
    }

    if (this->mp_mask)
    {
        /* Sometimes, the mask is not corresponding to the seed, like
         using the mean shape as the mask and just some other seed as
         the samples for feature. In such cases, do not touch mask
         and just return. */

        return;
    }


    long n = m_seeds.size();
    if (n == 0)
    {
        std::cerr << "Error: No seeds specified." << std::endl;
        abort();;
    }

    mp_mask = vtkImageData::New();
    mp_mask->SetDimensions(mp_img->GetDimensions());
    mp_mask->SetOrigin(mp_img->GetOrigin());
    mp_mask->SetSpacing(mp_img->GetSpacing());
    mp_mask->SetInformation(mp_img->GetInformation());
#if VTK_MAJOR_VERSION <= 5
    //    std::cout<<"VTK_MAJOR_VERSION <= 5\n";
    mp_mask->SetNumberOfScalarComponents(1);
    mp_mask->SetScalarTypeToUnsignedChar(); // MaskPixelType is uchar
    mp_mask->AllocateScalars();
#else
    mp_mask->AllocateScalars(VTK_UNSIGNED_CHAR, 1); // MaskPixelType is uchar
#endif

    for (long i = 0; i < n; ++i)
    {
        if (3 != m_seeds[i].size())
        {
            std::cerr<<"Error: 3 != m_seeds[i].size()\n";
            abort();
        }

        long ix = m_seeds[i][0];
        long iy = m_seeds[i][1];
        long iz = m_seeds[i][2];

        for (long iiz = iz - 1; iiz <= iz + 1; ++iiz)
        {
            for (long iiy = iy - 1; iiy <= iy + 1; ++iiy)
            {
                for (long iix = ix - 1; iix <= ix + 1; ++iix)
                {
                    if (0 <= iix && iix < this->m_nx && 0 <= iiy && iiy < this->m_ny && 0 <= iiz && iiz < this->m_nz)
                    {
                        //                      TIndex idx = {iix, iiy, iiz};
                        //                      this->mp_mask->SetPixel(idx, 1);
                        *(static_cast<MaskPixelType*>(mp_mask->GetScalarPointer(iix, iiy, iiz))) = 1;
                    }
                }
            }
        }
    }

    return;
}


/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap
::dialteSeeds()
{
    /* For each seed, add its 26 neighbors into the seed list. */

    if (!(this->mp_img))
    {
        std::cerr<<"Error: set input image first.\n";
        abort();
    }



    long n = m_seeds.size();
    std::vector<std::vector<long> > newSeeds;

    if (n == 0)
    {
        std::cerr << "Error: No seeds specified." << std::endl;
        abort();
    }


    for (long i = 0; i < n; ++i)
    {
        if (3 != m_seeds[i].size())
        {
            std::cerr<<"Error: 3 != m_seeds[i].size()\n";
            abort();
        }

        long ix = m_seeds[i][0];
        long iy = m_seeds[i][1];
        long iz = m_seeds[i][2];

        for (long iiz = iz - 1; iiz <= iz + 1; ++iiz)
        {
            for (long iiy = iy - 1; iiy <= iy + 1; ++iiy)
            {
                for (long iix = ix - 1; iix <= ix + 1; ++iix)
                {
                    if (0 <= iix && iix < this->m_nx && 0 <= iiy && iiy < this->m_ny && 0 <= iiz && iiz < this->m_nz)
                    {
                        /* Some locations may be added multiple times,
                         if the original seeds are close. But I think
                         this is fine */

                        std::vector<long> s(3);
                        s[0] = iix;
                        s[1] = iiy;
                        s[2] = iiz;

                        newSeeds.push_back(s);
                    }
                }
            }
        }
    }

    m_seeds.assign(newSeeds.begin(), newSeeds.end() );

    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::getFeatureAroundSeeds()
{
    if (!m_featureImageList[m_numberOfFeature-1])
    {
        // last feature image is not constructed
        std::cerr<<"Error: construct feature images first.\n";
        abort();;
    }

    long n = m_seeds.size();
    if (n == 0)
    {
        std::cerr << "Error: No seeds specified." << std::endl;
        abort();;
    }

    short ax = 0;
    short ay = 0;
    short az = 0;

    for (long i = 0; i < n; ++i)
    {
        if (3 != m_seeds[i].size())
        {
            std::cerr<<"Error: 3 != m_seeds[i].size()\n";
            abort();
        }

        long ix = m_seeds[i][0];
        long iy = m_seeds[i][1];
        long iz = m_seeds[i][2];

        for (long iiz = iz - az; iiz <= iz + az; ++iiz)
        {
            for (long iiy = iy - ay; iiy <= iy + ay; ++iiy)
            {
                for (long iix = ix - ax; iix <= ix + ax; ++iix)
                {
                    if (0 <= iix && iix < this->m_nx && 0 <= iiy && iiy < this->m_ny && 0 <= iiz && iiz < this->m_nz)
                    {
                        int idx[] = {iix, iiy, iiz};

                        std::vector<FeatureImagePixelType> featureHere(m_numberOfFeature);
                        computeFeatureAt(idx, featureHere);

                        m_featureAtTheSeeds.push_back(featureHere);
                    }
                }
            }
        }
    }


    return;
}

/* ============================================================ */

void
CSFLSRobustStatSegmentor3DLabelMap::estimateFeatureStdDevs()
{
    m_kernelStddev.assign(m_numberOfFeature, 0.0);

    long n = m_seeds.size(); // == m_featureAtTheSeeds.size()

    for (long i = 0; i < m_numberOfFeature; ++i)
    {
        double m = 0;
        for (long ii = 0; ii < n; ++ii)
        {
            m += m_featureAtTheSeeds[ii][i];
        }
        m /= n;

        for (long ii = 0; ii < n; ++ii)
        {
            m_kernelStddev[i] += (m_featureAtTheSeeds[ii][i] - m)*(m_featureAtTheSeeds[ii][i] - m);
        }

        m_kernelStddev[i] /= (n-1);
        m_kernelStddev[i] = sqrt(m_kernelStddev[i]);
    }

    return;
}


double
CSFLSRobustStatSegmentor3DLabelMap::kernelEvaluationUsingPDF(const std::vector<FeatureImagePixelType>& newFeature)
{
    double p = 1;

    for (long i = 0; i < m_numberOfFeature; ++i)
    {
        long idx = static_cast<long>(newFeature[i] - m_inputImageIntensityMin);

        double probOfThisFeature = m_PDFlearnedFromSeeds[i][idx];

        p *= probOfThisFeature;
    }

    return p;
}

/* ============================================================  */

double
CSFLSRobustStatSegmentor3DLabelMap::kernelEvaluation(const std::vector<FeatureImagePixelType>& newFeature)
{
    long n = m_seeds.size(); // == m_featureAtTheSeeds.size()

    double p = 1;
    //double p = 0;

    for (long i = 0; i < m_numberOfFeature; ++i)
    {
        double pp = 0.0;

        double stdDev = m_kernelStddev[i]/m_kernelWidthFactor; // /10 as in Eric's appendix

        double var2 = -1.0/(2*stdDev*stdDev);
        double c = 1.0/sqrt(2*m_pi)/stdDev;

        for (long ii = 0; ii < n; ++ii)
        {
            pp += exp(var2*(newFeature[i] - m_featureAtTheSeeds[ii][i])*(newFeature[i] - m_featureAtTheSeeds[ii][i]));
        }

        pp *= c;
        pp /= n;

        p *= pp;
        //p = p>pp?p:pp;
    }

    return p;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::setKernelWidthFactor(double f)
{
    if (f < 0.3)
    {
        m_kernelWidthFactor = 0.3;
    }

    if (f > 30.0)
    {
        m_kernelWidthFactor = 30.0;
    }

    m_kernelWidthFactor = f;


    //   std::ofstream fil("/tmp/d.txt", std::ios_base::app);
    //   fil<<"m_kernelWidthFactor = "<<m_kernelWidthFactor<<std::endl;
    //   fil.close();


    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::setIntensityHomogeneity(double h)
{
    //   std::ofstream fil("/tmp/d.txt", std::ios_base::app);
    //   fil<<"intensity homogeneity = "<<h<<std::endl;
    //   fil.close();


    double f = h*(30.0 - 0.3) + 0.3;

    setKernelWidthFactor(f);

    return;
}

/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::estimatePDFs()
{
    m_PDFlearnedFromSeeds.clear();

    computeMinMax(); // so we have the range of all pdfs

    long n = m_seeds.size();

    for (long ifeature = 0; ifeature < m_numberOfFeature; ++ifeature)
    {
        std::vector<double> thisPDF(m_inputImageIntensityMax - m_inputImageIntensityMin + 1);
        // assumption: TPixel are of integer types.

        double stdDev = m_kernelStddev[ifeature]/m_kernelWidthFactor; // /10 as in Eric's appendix
        double var2 = -1.0/(2*stdDev*stdDev);
        double c = 1.0/sqrt(2*m_pi)/stdDev;

        //#pragma omp parallel for
        for (TPixel a = m_inputImageIntensityMin; a <= m_inputImageIntensityMax; ++a)
        {
            long ia = static_cast<long>(a - m_inputImageIntensityMin);

            double pp = 0.0;
            for (long ii = 0; ii < n; ++ii)
            {
                pp += exp(var2*(a - m_featureAtTheSeeds[ii][ifeature])*(a - m_featureAtTheSeeds[ii][ifeature]));
            }

            pp *= c;
            pp /= n;

            thisPDF[ia] = pp;
        }


        m_PDFlearnedFromSeeds.push_back(thisPDF);
    }

    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::computeMinMax()
{
    if (!(this->mp_img))
    {
        std::cerr<<"Error: set input image first.\n";
        abort();;
    }

    m_inputImageIntensityMin = std::numeric_limits<TPixel>::max(); // yes, it's twisted so easity to compute.
    m_inputImageIntensityMax = std::numeric_limits<TPixel>::min();

    int* size = mp_img->GetDimensions();
    long n = size[0]*size[1]*size[2];

    for (long i = 0; i < n; i += m_increment0)
    {
        TPixel v = mp_img_ptr[i];
        m_inputImageIntensityMin = m_inputImageIntensityMin<v?m_inputImageIntensityMin:v;
        m_inputImageIntensityMax = m_inputImageIntensityMax>v?m_inputImageIntensityMax:v;
    }


    //    for (int iz = 0; iz < size[2]; ++iz)
    //    {
    //        for (int iy = 0; iy < size[1]; ++iy)
    //        {
    //            for (int ix = 0; ix < size[0]; ++ix)
    //            {
    //                 = *(static_cast<TPixel*>(mp_img->GetScalarPointer(ix, iy, iz)));


    //            }
    //        }
    //    }



    return;
}


#endif
