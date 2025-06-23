# require a sane minimum version of cmake
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

if(APPLE)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
endif()

# build debug by default with single-config generators
# must be done *before* first call to project()
set(CMAKE_BUILD_TYPE_INIT Debug)

# disable in-source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# setup base project
project(SimulationFramework)

# provide a list of config choices for single-config generators when using cmake-gui/ccmake
if(DEFINED CMAKE_BUILD_TYPE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release MinSizeRel RelWithDebInfo)
endif()

# set up sane default C++ standard settings for all projects
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set up sane warning defaults
if(MSVC)
    foreach(LANG C CXX)
        string(REGEX REPLACE "[-/]W[0-9]" "" CMAKE_${LANG}_FLAGS "${CMAKE_${LANG}_FLAGS}")
        string(APPEND CMAKE_${LANG}_FLAGS " /W4 /WX /wd4100 /wd4189 /wd4127 /wd4714 /wd4800 /D_CRT_SECURE_NO_WARNINGS")
    endforeach()
    string(APPEND CMAKE_CXX_FLAGS " /D_SCL_SECURE_NO_WARNINGS")
    foreach(TYPE EXE STATIC SHARED MODULE)
        string(APPEND CMAKE_${TYPE}_LINKER_FLAGS " /WX")
    endforeach()
else()
    foreach(LANG C CXX)
        string(APPEND CMAKE_${LANG}_FLAGS " -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter -Wno-unused-variable")
    endforeach()
    string(APPEND CMAKE_CXX_FLAGS " -Wno-attributes") # old GCC doesn't support [nodiscard]
endif()

# tell cmake where to find included .cmake files
list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")

# setup glad to generate code for OpenGL 3.3
add_subdirectory(glad/cmake)
glad_add_library(
    glad STATIC
    LOADER
    REPRODUCIBLE API gl:core=3.3
    EXTENSIONS
    GL_EXT_texture_compression_s3tc
    GL_EXT_texture_sRGB
    GL_EXT_texture_filter_anisotropic
)

# add OpenGLObjects library
add_subdirectory(OpenGLObjects)

# set up folders for better manageability in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_target_properties(
    glad OpenGLObjects
    PROPERTIES
    FOLDER Libraries
)

# setup CMAKE_PREFIX_PATH to use included copy of Eigen by default
list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/eigen/share/eigen3/cmake")
find_package(Eigen3 REQUIRED)

# find required Qt components
set(_QT_COMPONENTS Core Gui Widgets OpenGLWidgets)
include(SetupQtPrefixPath)
find_package(Qt6 6.2 REQUIRED COMPONENTS ${_QT_COMPONENTS} CONFIG)
if(Qt6_VERSION VERSION_GREATER_EQUAL 6.3)
    qt_standard_project_setup()
else()
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    find_package(Qt6 6.2 REQUIRED COMPONENTS Tools CONFIG)
endif()

# add an executable target and make it the default debug/startup project on VS
qt_add_executable(${PROJECT_NAME})
set_directory_properties(
    PROPERTIES
    VS_STARTUP_PROJECT ${PROJECT_NAME}
)

if(WIN32)
    # not a console application
    set_property(TARGET ${PROJECT_NAME} PROPERTY WIN32_EXECUTABLE ON)

    # add a pretty executable icon
    target_sources(
        ${PROJECT_NAME}
        PRIVATE
        icon.rc
    )
endif()

# set target sources
target_sources(
    ${PROJECT_NAME}
    PRIVATE
    main.cpp
    OpenGLWidget.cpp OpenGLWidget.hpp
    OpenGLRenderer.hpp
    GLMainWindow.cpp GLMainWindow.hpp GLMainWindow.ui
    ExampleRenderer.cpp ExampleRenderer.hpp
    Prog1_1Simulation.cpp Prog1_1Simulation.hpp
    Prog1_2Simulation.cpp Prog1_2Simulation.hpp
    Prog1_3Simulation.cpp Prog1_3Simulation.hpp
	Prog2_1Simulation.cpp Prog2_1Simulation.hpp
    helpers.hpp helpers.cpp
    constants.hpp
    shaders/icosphere.vert shaders/icosphere.frag
	shaders/lambertShader.vert shaders/lambertShader.frag
    shaders/skybox.vert shaders/skybox.frag
    icon.qrc
    textures.qrc
)

# set target requirements
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    Eigen3::Eigen
    glad
    OpenGLObjects
)

# set install folder
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

# link required Qt components
foreach(_COMP ${_QT_COMPONENTS})
    target_link_libraries(
        ${PROJECT_NAME}
        PRIVATE
        Qt6::${_COMP}
    )
endforeach()

# # copy/install required dlls (Windows only)
if(WIN32)
    target_include_directories(${PROJECT_NAME} PRIVATE ${Qt6Gui_PRIVATE_INCLUDE_DIRS})
    add_custom_command(
        TARGET ${PROJECT_NAME}
        COMMAND
            Qt6::windeployqt
            $<TARGET_FILE_NAME:${PROJECT_NAME}>
            --no-translations
            --no-opengl-sw
            --no-system-d3d-compiler
            --no-compiler-runtime
        WORKING_DIRECTORY $<TARGET_FILE_DIR:${PROJECT_NAME}>
        POST_BUILD
    )
endif()

# setup automatic processing of Qt files
set_target_properties(${PROJECT_NAME} PROPERTIES AUTORCC ON)
set_property(GLOBAL PROPERTY AUTOGEN_SOURCE_GROUP "Generated Files")

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}
)

source_group("Form Files" REGULAR_EXPRESSION "\\.ui$")
source_group("Resource Files" REGULAR_EXPRESSION "\\.q?rc$")

# setup shaders as qrc
include(ShadersToQRC)
shaders_to_qrc()

add_subdirectory(tests)
