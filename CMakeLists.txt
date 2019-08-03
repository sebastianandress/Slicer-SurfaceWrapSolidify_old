cmake_minimum_required(VERSION 3.5)

project(SegmentEditorSRSFilter)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/sebastianandress/Slicer-SegmentEditorSRSFilter")
set(EXTENSION_CATEGORY "Segmentations")
set(EXTENSION_CONTRIBUTORS "Sebastian Andress (LMU Munich)")
set(EXTENSION_DESCRIPTION "This is an example of a simple extension")
set(EXTENSION_ICONURL "https://github.com/sebastianandress/Slicer-SegmentEditorSRSFilter/raw/master/SegmentEditorSRSFilter.png")
set(EXTENSION_SCREENSHOTURLS "https://github.com/sebastianandress/Slicer-SegmentEditorSRSFilter/raw/master/Resources/Screenshots/screenshot1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(SRSFilter)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
