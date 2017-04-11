#include "SFLSRobustStatSegmentor3DLabelMap.h"

#include <algorithm>
#include <ctime>
#include <cstdlib>

#include <limits>

#include "vtkImageData.h"
#include "vtkImageCast.h"


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
    castFilter->SetInputData(l);
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

    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::computeForce()
{
    double fmax = std::numeric_limits<double>::min();
    double kappaMax = std::numeric_limits<double>::min();

    long n = this->m_lz.size();
    double* kappaOnZeroLS = new double[ n ];
    double* cvForce = new double[ n ];

    {
        long i = 0;

        for (CSFLSLayer::iterator itz = this->m_lz.begin(); itz != this->m_lz.end(); ++itz, ++i)
        {
            long ix = itz->SFLSNodeComponent1;
            long iy = itz->SFLSNodeComponent2;
            long iz = itz->SFLSNodeComponent3;

            int idx[] = {ix, iy, iz};

            kappaOnZeroLS[i] = this->computeKappa(ix, iy, iz);

            std::vector<FeatureImagePixelType> f(m_numberOfFeature);

            computeFeatureAt(idx, f);

            double a = -kernelEvaluationUsingPDF(f);

            fmax = fmax>fabs(a)?fmax:fabs(a);
            kappaMax = kappaMax>fabs(kappaOnZeroLS[i])?kappaMax:fabs(kappaOnZeroLS[i]);

            cvForce[i] = a;
        }
    }


    this->m_force.resize(n);
    for (long i = 0; i < n; ++i)
    {
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

    std::cout<<"in-inputLableImageToSeeds "<<size[0]<<'\t'<<size[1]<<'\t'<<size[2]<<"\n";

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
        fimg->AllocateScalars(VTK_FLOAT, 1);

        {
            /// init fimg to 0
            FeatureImagePixelType* fimg_ptr = static_cast<FeatureImagePixelType*>(fimg->GetScalarPointer(0, 0, 0));
            long n = (this->m_nx)*(this->m_ny)*(this->m_nz);
            for (long i = 0; i < n; i += m_increment0)
            {
                fimg_ptr[i] = 0.0;
            }
        }

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
    m_featureComputed->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

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

    if (m_featureComputed_idx_ptr[0])
    {
        // the feature at this pixel is computed, just retrive
        for (long i = 0; i < m_numberOfFeature; ++i)
        {
            f[i] = *(static_cast<FeatureImagePixelType*>((m_featureImageList[i])->GetScalarPointer(idx)));
        }
    }
    else
    {
        // compute the feature
        std::vector< TPixel > neighborIntensities;

        long ix = idx[0];
        long iy = idx[1];
        long iz = idx[2];

        TPixel* img_ptr = mp_img_ptr + iz*m_increment2 + iy*m_increment1 + ix;

        for (long iiz = - m_statNeighborZ; iiz <= m_statNeighborZ; ++iiz)
        {
            for (long iiy = - m_statNeighborY; iiy <= m_statNeighborY; ++iiy)
            {
                for (long iix = - m_statNeighborX; iix <= m_statNeighborX; ++iix)
                {
                    if (0 <= ix-iix && ix+iix < this->m_nx && 0 <= iy-iiy && iy+iiy < this->m_ny && 0 <= iz-iiz && iz+iiz < this->m_nz)
                    {
                        TPixel tmp = img_ptr[iiz*m_increment2 + iiy*m_increment1 + iix*m_increment0];
                        neighborIntensities.push_back(tmp);
                    }
                }
            }
        }

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

    /*============================================================
     * From the initial mask, generate: 1. SFLS, 2. mp_label and
     * 3. mp_phi.
     */
    this->initializeSFLS();
}


void CSFLSRobustStatSegmentor3DLabelMap::inOneSegmentationIteration()
{
    if (m_done)
    {
        return;
    }

    float oldVoxelCount = this->m_insideVoxelCount;

    computeForce();

    this->normalizeForce();

    this->oneStepLevelSetEvolution();


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
    mp_mask->AllocateScalars(VTK_UNSIGNED_CHAR, 1); // MaskPixelType is uchar

    {
        /// init mp_mask to 0
        MaskPixelType* mp_mask_pixel_ptr = static_cast<MaskPixelType*>(mp_mask->GetScalarPointer(0, 0, 0));

        long nn = (this->m_nx)*(this->m_ny)*(this->m_nz);
        for (long i = 0; i < nn; i += m_increment0)
        {
            mp_mask_pixel_ptr[i] = 0;
        }
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
                        *(static_cast<MaskPixelType*>(mp_mask->GetScalarPointer(iix, iiy, iiz))) = 1;
                    }
                }
            }
        }
    }

    return;
}


///* ============================================================ */

//void
//CSFLSRobustStatSegmentor3DLabelMap
//::dialteSeeds()
//{
//    /* For each seed, add its 26 neighbors into the seed list. */

//    if (!(this->mp_img))
//    {
//        std::cerr<<"Error: set input image first.\n";
//        abort();
//    }



//    long n = m_seeds.size();
//    std::vector<std::vector<long> > newSeeds;

//    if (n == 0)
//    {
//        std::cerr << "Error: No seeds specified." << std::endl;
//        abort();
//    }


//    for (long i = 0; i < n; ++i)
//    {
//        if (3 != m_seeds[i].size())
//        {
//            std::cerr<<"Error: 3 != m_seeds[i].size()\n";
//            abort();
//        }

//        long ix = m_seeds[i][0];
//        long iy = m_seeds[i][1];
//        long iz = m_seeds[i][2];

//        for (long iiz = iz - 1; iiz <= iz + 1; ++iiz)
//        {
//            for (long iiy = iy - 1; iiy <= iy + 1; ++iiy)
//            {
//                for (long iix = ix - 1; iix <= ix + 1; ++iix)
//                {
//                    if (0 <= iix && iix < this->m_nx && 0 <= iiy && iiy < this->m_ny && 0 <= iiz && iiz < this->m_nz)
//                    {
//                        /* Some locations may be added multiple times,
//                         if the original seeds are close. But I think
//                         this is fine */

//                        std::vector<long> s(3);
//                        s[0] = iix;
//                        s[1] = iiy;
//                        s[2] = iiz;

//                        newSeeds.push_back(s);
//                    }
//                }
//            }
//        }
//    }

//    m_seeds.assign(newSeeds.begin(), newSeeds.end() );

//    return;
//}


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

    return;
}


/* ============================================================  */

void
CSFLSRobustStatSegmentor3DLabelMap::setIntensityHomogeneity(double h)
{
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

    return;
}
