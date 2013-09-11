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

    //  typedef boost::shared_ptr< Self > Pointer;

    typedef SuperClassType::TPixel TPixel;
    typedef typename SuperClassType::NodeType NodeType;
    typedef typename SuperClassType::CSFLSLayer CSFLSLayer;

    /*================================================================================
    ctor */
    CSFLSRobustStatSegmentor3DLabelMap() : CSFLSSegmentor3D()
    {
        basicInit();
    }

    //  /* New */
    //  static Pointer New() { return Pointer(new Self); }

    void basicInit();

    /* just copy, not logic change
     ----------------------------------------------------------------------
     ----------------------------------------------------------------------
     ----------------------------------------------------------------------
     ---------------------------------------------------------------------- */

//    typedef SuperClassType::vtkImageDataPointer vtkImageDataPointer;
    typedef float FeatureImagePixelType;
    //  typedef typename SuperClassType::TCharImage TLabelImage;
    //  typedef typename TLabelImage::Pointer TLabelImagePointer;


    //  typedef typename SuperClassType::TDoubleImage TDoubleImage;
    //  typedef typename TDoubleImage::Pointer TDoubleImagePointer;

    //  typedef typename SuperClassType::TImage TImage;
    //  typedef typename SuperClassType::TFloatImage TFloatImage;


    //  typedef typename SuperClassType::MaskImageType TMaskImage;

    //  typedef typename SuperClassType::TIndex TIndex;
    //  typedef typename SuperClassType::TSize TSize;
    //  typedef typename SuperClassType::TRegion TRegion;

    /* ============================================================
   * functions
   * ============================================================*/

    void setInputLabelImage(vtkImageData* l);

  //void doSegmenation();

  void doSegmenationBeforeIteration();
  //void doSegmenationIteration();
  void inOneSegmentationIteration();

    virtual void computeForce();

    void setKernelWidthFactor(double f);
    void setIntensityHomogeneity(double h);

protected:
    /* data */
    //  TLabelImagePointer m_inputLabelImage;
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
    //  TLabelImagePointer m_featureComputed; // if feature at this point is computed, then is 1
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
    void dialteSeeds();
    void getFeatureAroundSeeds();
    void estimateFeatureStdDevs();


    TPixel m_inputImageIntensityMin;
    TPixel m_inputImageIntensityMax;
    void computeMinMax();

    std::vector< std::vector<double> > m_PDFlearnedFromSeeds; // each feature corresponds to a inner std::vector<double>
    void estimatePDFs();

    //void getFeatureAt(TDoubleImage::IndexType idx, std::vector<double>& f);

    virtual void getThingsReady();


    // kernel
    std::vector<double> m_kernelStddev;
    double kernelEvaluation(const std::vector<FeatureImagePixelType>& newFeature);
    double kernelEvaluationUsingPDF(const std::vector<FeatureImagePixelType>& newFeature);

};


#endif
