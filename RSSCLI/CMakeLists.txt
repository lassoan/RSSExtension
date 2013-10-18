
#-----------------------------------------------------------------------------
set(MODULE_NAME RSSCLI)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  )

set(MODULE_SRCS
    labelMapPreprocessor.h
    SFLS.h
    SFLSSegmentor3D.h
    SFLSSegmentor3D.txx
    SFLSRobustStatSegmentor3DLabelMap_single.h
    SFLSRobustStatSegmentor3DLabelMap_single.txx
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

