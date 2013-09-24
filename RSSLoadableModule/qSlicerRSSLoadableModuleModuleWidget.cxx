/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QDebug>

// SlicerQt includes
#include "qSlicerRSSLoadableModuleModuleWidget.h"
#include "ui_qSlicerRSSLoadableModuleModuleWidget.h"

// add by YG
#include "vtkMRMLNode.h"
#include "vtkMRMLScalarVolumeNode.h"
#include "vtkImageData.h"
//#include "vtkITKImageWriter.h"
//#include "vtkSmartPointer.h"
#include "vtkImageCast.h"

#include <iostream>

#include <cstdlib>


#include "SFLSRobustStatSegmentor3DLabelMap.h"

//#include "walltime.h"

#include "vtkVersion.h"
//#include "vtkSmartPointer.h"
#include <vtkImageData.h>

#include "vtkSlicerApplicationLogic.h"
#include "vtkMRMLSelectionNode.h"
#include "qSlicerApplication.h"

#include "vtkITKImageWriter.h"

#include "QTimer"

#include "vtkMRMLScene.h"

#include "qSlicerLayoutManager.h"
#include "vtkMRMLSliceNode.h"
#include "qMRMLSliceWidget.h"
#include "qMRMLSliceView.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

// dbg
#include "vtkMetaImageWriter.h"


//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerRSSLoadableModuleModuleWidgetPrivate: public Ui_qSlicerRSSLoadableModuleModuleWidget
{
public:
    qSlicerRSSLoadableModuleModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerRSSLoadableModuleModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerRSSLoadableModuleModuleWidgetPrivate::qSlicerRSSLoadableModuleModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerRSSLoadableModuleModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerRSSLoadableModuleModuleWidget::qSlicerRSSLoadableModuleModuleWidget(QWidget* _parent)
    : Superclass( _parent )
    , d_ptr( new qSlicerRSSLoadableModuleModuleWidgetPrivate )
{
    m_evolutionPaused = 0;

    //    m_UIRefreshInterval = 201;
}

//-----------------------------------------------------------------------------
qSlicerRSSLoadableModuleModuleWidget::~qSlicerRSSLoadableModuleModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerRSSLoadableModuleModuleWidget::setup()
{
    Q_D(qSlicerRSSLoadableModuleModuleWidget);
    d->setupUi(this);
    this->Superclass::setup();

    connect(d->InputVolumeMRMLNodeComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(onInputVolumeChanged(vtkMRMLNode*)));
    connect(d->InputLabelVolumeMRMLNodeComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(onInputLabelVolumeChanged(vtkMRMLNode*)));
    connect(d->OutputLabelVolumeMRMLNodeComboBox, SIGNAL(currentNodeChanged(vtkMRMLNode*)), this, SLOT(onOutputVolumeChanged(vtkMRMLNode*)));

    connect(d->applyButton, SIGNAL(clicked()), this, SLOT(applyPushButtonClicked()));
    connect(d->pauseButton, SIGNAL(clicked()), this, SLOT(pausePushButtonClicked()));


    d->SmoothnessSliderWidget->setSynchronizeSiblings(ctkSliderWidget::NoSynchronize);
    d->labelValueSliderWidget->setSynchronizeSiblings(ctkSliderWidget::NoSynchronize);

    d->SmoothnessSliderWidget->setDecimals(1);
    d->labelValueSliderWidget->setDecimals(0);

    d->SmoothnessSliderWidget->setValue(0.3);
}

void qSlicerRSSLoadableModuleModuleWidget::onInputVolumeChanged(vtkMRMLNode* node)
{
    if (node)
    {
        qDebug("onInputVolumeChanged");
        qDebug(node->GetName(), 100);

        vtkMRMLScalarVolumeNode* scalarVolNode = vtkMRMLScalarVolumeNode::SafeDownCast(node);
        vtkImageData* vtkImg = scalarVolNode->GetImageData();

        double spacing[3];
        vtkImg->GetSpacing(spacing);
        std::cout<<spacing[0]<<'\t'<<spacing[1]<<'\t'<<spacing[2]<<std::endl;

        int extent[6];
        vtkImg->GetExtent(extent);
        std::cout<<extent[0]<<'\t'<<extent[1]<<'\t'<<extent[2]<<'\t'<<extent[3]<<'\t'<<extent[4]<<'\t'<<extent[5]<<std::endl;
    }
}


void qSlicerRSSLoadableModuleModuleWidget::onInputLabelVolumeChanged(vtkMRMLNode* node)
{
    if (node)
    {
        qDebug("onInputLabelVolumeChanged");
        qDebug(node->GetName(), 100);

        vtkMRMLScalarVolumeNode* scalarVolNode = vtkMRMLScalarVolumeNode::SafeDownCast(node);
        vtkImageData* vtkImg = scalarVolNode->GetImageData();

        double spacing[3];
        vtkImg->GetSpacing(spacing);
        std::cout<<spacing[0]<<'\t'<<spacing[1]<<'\t'<<spacing[2]<<std::endl;

        int extent[6];
        vtkImg->GetExtent(extent);
        std::cout<<extent[0]<<'\t'<<extent[1]<<'\t'<<extent[2]<<'\t'<<extent[3]<<'\t'<<extent[4]<<'\t'<<extent[5]<<std::endl;
    }
}


void qSlicerRSSLoadableModuleModuleWidget::onOutputVolumeChanged(vtkMRMLNode* node)
{
    if (node)
    {
        qDebug("onOutputVolumeChanged");
        qDebug(node->GetName(), 100);
        qDebug(node->GetClassName(), 100);

    }
}


void qSlicerRSSLoadableModuleModuleWidget::pausePushButtonClicked()
{
    Q_D(qSlicerRSSLoadableModuleModuleWidget);
    if (!m_evolutionPaused)
    {
        d->pauseButton->setText("Continue");
        m_evolutionPaused = true;
    }
    else
    {
        d->pauseButton->setText("Pause");
        m_evolutionPaused = false;

        QTimer::singleShot(m_UIRefreshInterval, this, SLOT(oneLevelSetIteration()));
    }

    return;
}



void qSlicerRSSLoadableModuleModuleWidget::applyPushButtonClicked()
{
    std::cout<<"PushButton clicked.\n";

    Q_D(qSlicerRSSLoadableModuleModuleWidget);

    double curvatureWeight = d->SmoothnessSliderWidget->value();
    double intensityHomogeneity = 0.6;
    short labelValue = d->labelValueSliderWidget->value();

    std::cout<<"curvatureWeight = "<<curvatureWeight<<std::endl;
    std::cout<<"intensityHomogeneity = "<<intensityHomogeneity<<std::endl;
    std::cout<<"labelValue = "<<labelValue<<std::endl;

    vtkMRMLNode* inputNode = d->InputVolumeMRMLNodeComboBox->currentNode();
    vtkMRMLVolumeNode* inputVolumeNode = vtkMRMLVolumeNode::SafeDownCast(inputNode);
    vtkImageCast* castFilter = vtkImageCast::New();
    castFilter->SetInput(inputVolumeNode->GetImageData());
    castFilter->SetOutputScalarTypeToShort();// coz RSS input label is short pixel type
    castFilter->Update();
    vtkImageData* inputImageVTK = castFilter->GetOutput();
    inputImageVTK->SetSpacing(inputVolumeNode->GetSpacing() );
    inputImageVTK->SetOrigin(inputVolumeNode->GetOrigin());


//    double dx, dy, dz;
//    inputVolumeNode->GetSpacing(dx, dy, dz);
//    std::cout<<"dx, dy, dz = "<<dx<<'\t'<<dy<<'\t'<<dz<<'\n';


    vtkMRMLNode* inputLabelNode = d->InputLabelVolumeMRMLNodeComboBox->currentNode();    
    vtkMRMLVolumeNode* inputLabelVolumeNode = vtkMRMLVolumeNode::SafeDownCast(inputLabelNode);

    vtkImageData* inputLableImageVtkNoOrientationInfo = inputLabelVolumeNode->GetImageData();
    inputLableImageVtkNoOrientationInfo->SetOrigin(inputImageVTK->GetOrigin());
    inputLableImageVtkNoOrientationInfo->SetSpacing(inputImageVTK->GetSpacing());
    inputLableImageVtkNoOrientationInfo->SetInformation(inputImageVTK->GetInformation());

//    inputLabelVolumeNode->SetOrigin(inputVolumeNode->GetOrigin()); % these two lines does not add spatial/spacing info to the newly drawn label images.
//    inputLabelVolumeNode->CopyOrientation(inputVolumeNode);
    vtkImageCast* castFilter1 = vtkImageCast::New();
    castFilter1->SetInput(inputLableImageVtkNoOrientationInfo);
    castFilter1->SetOutputScalarTypeToShort(); // coz RSS input label is short pixel type
    castFilter1->Update();
    //    vtkImageData* inputLabelImageVtk = inputLabelVolumeNode->GetImageData();
    vtkImageData* newInputLabelImageVtk = preprocessLabelMap(castFilter1->GetOutput(), labelValue);
    newInputLabelImageVtk->SetSpacing(inputVolumeNode->GetSpacing() );
    newInputLabelImageVtk->SetOrigin(inputVolumeNode->GetOrigin());

    qSlicerApplication * app = qSlicerApplication::application();

    d->pauseButton->setEnabled(false);
    //d->applyButton->setEnabled(false);

    // do seg
    m_rssPointer = new CSFLSRobustStatSegmentor3DLabelMap();
    m_rssPointer->setImage(inputImageVTK);
    m_rssPointer->setInputLabelImage(newInputLabelImageVtk);

    m_rssPointer->setNumIter(10000); // a large enough number, s.t. will not be stoped by this creteria. because we have the pause button in loadable module
    double expectedVolume = 1000000; // a large enough number, s.t. will not be stoped by this creteria. because we have the pause button in loadable module
    m_rssPointer->setMaxVolume(expectedVolume);

    double maxRunningTime = 10000; // a large enough number, s.t. will not be stoped by this creteria. because we have the pause button in loadable module
    m_rssPointer->setMaxRunningTime(maxRunningTime);

    m_rssPointer->setIntensityHomogeneity(intensityHomogeneity);
    m_rssPointer->setCurvatureWeight(curvatureWeight / 1.5);

    m_rssPointer->doSegmenationBeforeIteration();

    d->pauseButton->setEnabled(true);
    d->applyButton->setEnabled(false);


//    // dbg
//    {
//      char mhdName[1000];
//      sprintf(mhdName, "/tmp/mp_label_mask-here-%d.mhd", 1);
//      char rawName[1000];
//      sprintf(rawName, "/tmp/mp_label_mask-here-%d.raw", 1);

//      vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
//      writer->SetFileName(mhdName);
//      writer->SetRAWFileName(rawName);
//      writer->SetInput(m_rssPointer->mp_label_mask);
//      writer->Write();
//    }
//    // dbg, end

    vtkMRMLNode* outputNode = d->OutputLabelVolumeMRMLNodeComboBox->currentNode();
    vtkMRMLScalarVolumeNode* outputVolumeNode = vtkMRMLScalarVolumeNode::SafeDownCast(outputNode);
    outputVolumeNode->SetOrigin(inputVolumeNode->GetOrigin());
    outputVolumeNode->CopyOrientation(inputVolumeNode);

    outputVolumeNode->SetAndObserveImageData(m_rssPointer->mp_label_mask);
    outputVolumeNode->SetLabelMap(1);

    //    outputVolumeNode->SetAndObserveImageData(m_rssPointer->mp_label_mask);
    //    outputVolumeNode->SetAndObserveImageData(m_rssPointer->mp_label);

    outputVolumeNode->SetDisplayVisibility(1);


//    std::cout<<"11111111111111111111111111111 np = "<<m_rssPointer->mp_img->GetNumberOfPoints()<<std::endl<<std::flush;

    //    vtkSlicerApplicationLogic *appLogic = this->module()->appLogic();
    vtkSlicerApplicationLogic *appLogic = app->applicationLogic();
    vtkMRMLSelectionNode *selectionNode = appLogic->GetSelectionNode();
    selectionNode->SetActiveLabelVolumeID(outputVolumeNode->GetID());
    appLogic->PropagateVolumeSelection();

//    std::cout<<"22222222222222222222222222 np = "<<m_rssPointer->mp_img->GetNumberOfPoints()<<std::endl<<std::flush;

//    // dbg
//    {
//      char mhdName[1000];
//      sprintf(mhdName, "/tmp/mp_label_mask-here-%d.mhd", 3);
//      char rawName[1000];
//      sprintf(rawName, "/tmp/mp_label_mask-here-%d.raw", 3);

//      vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
//      writer->SetFileName(mhdName);
//      writer->SetRAWFileName(rawName);
//      writer->SetInput(m_rssPointer->mp_label_mask);
//      writer->Write();
//    }
//    // dbg, end

    oneLevelSetIteration();

    QTimer::singleShot(m_UIRefreshInterval, this, SLOT(oneLevelSetIteration()));
}


vtkImageData* qSlicerRSSLoadableModuleModuleWidget::preprocessLabelMap(vtkImageData* originalLabelMap, short desiredLabel)
{
    /*
    If there is a single non-zero label in the originalLabelMap, then
    use that as the label map to feed the segmentor.

    If originalLabelMap contains multiple labels, extract
    desiredLabel, forms a new label map containing only desiredLabel
    to feed the segmentor.

    If originalLabelMap contains multiple labels, but do NOT contain
    desiredLabel, then use all the non-zero labels as a single label
    and output label value is desiredLabel.

    1. count number of different labels in the originalLabelMap

    2. if #= 1, return it

    3. if #= 2, if desiredLabel is not there, then return originalLableMap itself

    4. if #= 2, go thru originalLabelMap and fill in a new lable map
       with 1 where the label matches.


   */

    // 1.
    int* size = originalLabelMap->GetDimensions();
    long n = size[0]*size[1]*size[2];

    int startIdx[] = {0, 0, 0};
    short* originalLabelMap_buffer_ptr = static_cast<short*>(originalLabelMap->GetScalarPointer(startIdx));

    std::vector<short> uniqueLabels(n);
    for (long i = 0; i < n; ++i)
    {
        uniqueLabels[i] = originalLabelMap_buffer_ptr[i];
    }

    std::sort(uniqueLabels.begin(), uniqueLabels.end() );
    std::vector<short>::iterator itl = std::unique(uniqueLabels.begin(), uniqueLabels.end() );
    uniqueLabels.resize( itl - uniqueLabels.begin() );

    if( uniqueLabels[0] != 0 )
    {
        std::cerr << "Error: least label is not 0? no background?\n";
        abort();
    }

    short numOfLabels = uniqueLabels.size() - 1; // 0 not count

    // 2.
    if( 1 == numOfLabels )
    {
        return originalLabelMap;
    }

    // 3.
    if( !std::binary_search(uniqueLabels.begin(), uniqueLabels.end(), desiredLabel) )
    {
        return originalLabelMap;
    }

    // 4.
    vtkImageData* newLabelMap = vtkImageData::New();
    newLabelMap->SetDimensions(size);
    newLabelMap->SetOrigin(originalLabelMap->GetOrigin());
    newLabelMap->SetSpacing(originalLabelMap->GetSpacing());
    newLabelMap->SetInformation(originalLabelMap->GetInformation());
    //    newLabelMap->SetExtent(originalLabelMap->GetExtent());
//#if VTK_MAJOR_VERSION <= 5
    //    std::cout<<"VTK_MAJOR_VERSION <= 5\n";
    newLabelMap->SetNumberOfScalarComponents(1);
    newLabelMap->SetScalarTypeToShort();
    newLabelMap->AllocateScalars();
//#else
//    newLabelMap->AllocateScalars(VTK_SHORT, 1);
//#endif


    short* newLabelMap_buffer_ptr = static_cast<short*>(newLabelMap->GetScalarPointer(startIdx));

    for (long i = 0; i < n; ++i)
    {
        newLabelMap_buffer_ptr[i] = originalLabelMap_buffer_ptr[i] == desiredLabel?1:0;
    }

    return newLabelMap;
}

void qSlicerRSSLoadableModuleModuleWidget::oneLevelSetIteration()
{
//    std::cout<<"3333333333333333333333333 np = "<<m_rssPointer->mp_label_mask->GetNumberOfPoints()<<std::endl<<std::flush;

    if (m_rssPointer->m_done || m_evolutionPaused)
    {
        return;
    }
    else
    {
//        std::cout<<"44444444444444444444 np = "<<m_rssPointer->mp_img->GetNumberOfPoints()<<std::endl<<std::flush;


        m_rssPointer->inOneSegmentationIteration();

        Q_D(qSlicerRSSLoadableModuleModuleWidget);
        vtkMRMLNode* outputNode = d->OutputLabelVolumeMRMLNodeComboBox->currentNode();
        outputNode->Modified();

                QTimer::singleShot(m_UIRefreshInterval, this, SLOT(oneLevelSetIteration()));
    }
}
