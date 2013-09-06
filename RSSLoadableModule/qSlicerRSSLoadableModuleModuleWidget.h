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

#ifndef __qSlicerRSSLoadableModuleModuleWidget_h
#define __qSlicerRSSLoadableModuleModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerRSSLoadableModuleModuleExport.h"

// added by YG
#include "itkImage.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "SFLSRobustStatSegmentor3DLabelMap.h"

class qSlicerRSSLoadableModuleModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_RSSLOADABLEMODULE_EXPORT qSlicerRSSLoadableModuleModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerRSSLoadableModuleModuleWidget(QWidget *parent=0);
  virtual ~qSlicerRSSLoadableModuleModuleWidget();

public slots:

  void onInputVolumeChanged(vtkMRMLNode* node);
  void onInputLabelVolumeChanged(vtkMRMLNode* node);
  void onOutputVolumeChanged(vtkMRMLNode* node);

  void applyPushButtonClicked();
  void pausePushButtonClicked();

  void oneLevelSetIteration();

protected:
  QScopedPointer<qSlicerRSSLoadableModuleModuleWidgetPrivate> d_ptr;
  
  virtual void setup();

static const int m_UIRefreshInterval = 101;

private:
  Q_DECLARE_PRIVATE(qSlicerRSSLoadableModuleModuleWidget);
  Q_DISABLE_COPY(qSlicerRSSLoadableModuleModuleWidget);

//  itk::Image<short, 3>::Pointer getFinalMask(itk::Image<float, 3>::Pointer img, unsigned char l, float thod);

//  vtkImageData* thresholdLevelSetFunctionToGetFinalLabelImage(vtkImageData* levelSetFn, float threshold, short targetLabel);
  vtkImageData* preprocessLabelMap(vtkImageData* originalLabelMap, short desiredLabel);

  CSFLSRobustStatSegmentor3DLabelMap* m_rssPointer;
  bool m_evolutionPaused;


};

#endif
