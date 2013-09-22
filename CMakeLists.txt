cmake_minimum_required(VERSION 2.8.9)

project(RSSExtension)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/RobustStatisticsSegmenter")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Yi Gao (UAB), Ron Kikinis (BWH), Sylvain Bouix (BWH), Martha Shenton (BWH), Allen Tannenbaum (SBU)")
set(EXTENSION_DESCRIPTION "User draw some initial seeds, this module perform an interative segmentation in 3D volume.")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/RSSExtension/RSSExtension.png?revision=21746&view=co")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/f/f1/MeningiomaCase1RSSLabel.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(RSSLoadableModule)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
