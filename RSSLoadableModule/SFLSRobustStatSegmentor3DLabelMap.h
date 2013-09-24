#ifndef SFLSRobustStatSegmentor3DLabelMap_h_
#define SFLSRobustStatSegmentor3DLabelMap_h_

#include "SFLSSegmentor3D.h"

#include <list>
#include <vector>



class CSFLSRobustStatSegmentor3DLabelMap : public CSFLSSegmentor3D
{
    /*----------------------------------------------------------------------
    just copy, not logic change */

public:
    typedef CSFLSSegmentor3D SuperClassType;

    typedef CSFLSRobustStatSegmentor3DLabelMap Self;

    typedef SuperClassType::TPixel TPixel;
    typedef SuperClassType::NodeType NodeType;
    typedef SuperClassType::CSFLSLayer CSFLSLayer;

    /*================================================================================
    ctor */
    CSFLSRobustStatSegmentor3DLabelMap() : CSFLSSegmentor3D()
    {
        basicInit();
    }

    void basicInit();

    /* just copy, not logic change
     ----------------------------------------------------------------------
     ----------------------------------------------------------------------
     ----------------------------------------------------------------------
     ---------------------------------------------------------------------- */
    typedef float FeatureImagePixelType;

    /* ============================================================
   * functions
   * ============================================================*/

    void setInputLabelImage(vtkImageData* l);

  //void doSegmenation();

    void doSegmenationBeforeIteration();
    void doSegmenationIteration();
    void inOneSegmentationIteration();

    virtual void computeForce();

    void setKernelWidthFactor(double f);
    void setIntensityHomogeneity(double h);

protected:
    /* data */
    vtkImageData* m_inputLabelImage;
    LabelImagePixelType* m_inputLabelImage_buffer_ptr;
    std::vector<std::vector<long> > m_seeds; // in IJK

    std::vector< std::vector<FeatureImagePixelType> > m_featureAtTheSeeds;

    long m_statNeighborX;
    long m_statNeighborY;
    long m_statNeighborZ;

    const static long m_numberOfFeature = 3;
    /* Store the robust stat as the feature at each point
     0: Meadian
     1: interquartile range (IRQ)
     2. median absolute deviation (MAD)
  */
    typedef unsigned char m_featureComputed_pixel_type;
    vtkImageData* m_featureComputed; // if feature at this point is computed, then is 1
    std::vector<vtkImageData*> m_featureImageList;


    double m_kernelWidthFactor; // kernel_width = empirical_std/m_kernelWidthFactor, Eric has it at 10.0


    /* fn */
    void initFeatureComputedImage();
    void initFeatureImage();

    //void computeFeature();
    void computeFeatureAt(int idx[], std::vector<FeatureImagePixelType>& f);

    void getRobustStatistics(std::vector<TPixel>& samples, std::vector<float>& robustStat);
    void inputLableImageToSeeds();
    void seedToMask();
//    void dialteSeeds();
    void getFeatureAroundSeeds();
    void estimateFeatureStdDevs();


    TPixel m_inputImageIntensityMin;
    TPixel m_inputImageIntensityMax;
    void computeMinMax();

    std::vector< std::vector<double> > m_PDFlearnedFromSeeds; // each feature corresponds to a inner std::vector<double>
    void estimatePDFs();

    virtual void getThingsReady();

    // kernel
    std::vector<double> m_kernelStddev;
    double kernelEvaluation(const std::vector<FeatureImagePixelType>& newFeature);
    double kernelEvaluationUsingPDF(const std::vector<FeatureImagePixelType>& newFeature);

};


#endif
