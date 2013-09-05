/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerRSSLoadableModuleFooBarWidget.h"
#include "ui_qSlicerRSSLoadableModuleFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_RSSLoadableModule
class qSlicerRSSLoadableModuleFooBarWidgetPrivate
  : public Ui_qSlicerRSSLoadableModuleFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerRSSLoadableModuleFooBarWidget);
protected:
  qSlicerRSSLoadableModuleFooBarWidget* const q_ptr;

public:
  qSlicerRSSLoadableModuleFooBarWidgetPrivate(
    qSlicerRSSLoadableModuleFooBarWidget& object);
  virtual void setupUi(qSlicerRSSLoadableModuleFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerRSSLoadableModuleFooBarWidgetPrivate
::qSlicerRSSLoadableModuleFooBarWidgetPrivate(
  qSlicerRSSLoadableModuleFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerRSSLoadableModuleFooBarWidgetPrivate
::setupUi(qSlicerRSSLoadableModuleFooBarWidget* widget)
{
  this->Ui_qSlicerRSSLoadableModuleFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerRSSLoadableModuleFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerRSSLoadableModuleFooBarWidget
::qSlicerRSSLoadableModuleFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerRSSLoadableModuleFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerRSSLoadableModuleFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerRSSLoadableModuleFooBarWidget
::~qSlicerRSSLoadableModuleFooBarWidget()
{
}
