if(NOT MSVC)
    return()
endif()

# use registry to find Qt https://stackoverflow.com/questions/15639781/how-to-find-the-qt5-cmake-module-on-windows
get_filename_component(_QT_CREATOR_BIN "[HKEY_CURRENT_USER\\Software\\Classes\\Applications\\QtProject.QtCreator.cpp\\shell\\Open\\Command]" PATH)
get_filename_component(_QT_ROOT "${_QT_CREATOR_BIN}\\..\\..\\.." ABSOLUTE)

file(GLOB _QT_VERSIONS "${_QT_ROOT}/6.*")
list(SORT _QT_VERSIONS COMPARE NATURAL ORDER DESCENDING)

set(_PLATFORM)
if(CMAKE_CL_64)
    set(_PLATFORM _64)
endif()

foreach(_QT_VERSION IN LISTS _QT_VERSIONS)
    set(_MSVC_YEAR)
    if((MSVC_TOOLSET_VERSION GREATER_EQUAL 142) AND (IS_DIRECTORY "${_QT_VERSION}/msvc2019${_PLATFORM}"))
        set(_MSVC_YEAR 2019)
    elseif((MSVC_TOOLSET_VERSION GREATER_EQUAL 141) AND (IS_DIRECTORY "${_QT_VERSION}/msvc2017${_PLATFORM}"))
        set(_MSVC_YEAR 2017)
    elseif(MSVC_TOOLSET_VERSION GREATER_EQUAL 140)
        set(_MSVC_YEAR 2015)
    endif()

    set(_PATH "${_QT_VERSION}/msvc${_MSVC_YEAR}${_PLATFORM}/lib/cmake/Qt6/")
    if(IS_DIRECTORY ${_PATH})
        list(APPEND CMAKE_PREFIX_PATH ${_PATH})
    endif()
endforeach()
