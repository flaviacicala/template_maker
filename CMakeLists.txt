# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.10)

# Project's name
project(PulseProject)

# Set the C++ Standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Find the ROOT package which provides the ROOT libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree Spectrum)

include(${ROOT_USE_FILE})

# Add the directory where your header files are stored
include_directories(${PROJECT_SOURCE_DIR}/include)

# Create variables with all the sources
set(SOURCES_PULSE_EXTRACTOR_INDUCTION_1
    ${PROJECT_SOURCE_DIR}/src/pulse_extractor_induction_1.cpp
    ${PROJECT_SOURCE_DIR}/src/wave_processing.cpp
    ${PROJECT_SOURCE_DIR}/src/waveform_utils.cpp
    ${PROJECT_SOURCE_DIR}/src/WaveformManager_2.cpp
    ${PROJECT_SOURCE_DIR}/src/DataLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/heatmap_multi.cpp
)
set(SOURCES_PULSE_EXTRACTOR_INDUCTION_2
    ${PROJECT_SOURCE_DIR}/src/pulse_extractor_induction_2.cpp
    ${PROJECT_SOURCE_DIR}/src/wave_processing.cpp
    ${PROJECT_SOURCE_DIR}/src/waveform_utils.cpp
    ${PROJECT_SOURCE_DIR}/src/WaveformManager_2.cpp
    ${PROJECT_SOURCE_DIR}/src/DataLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/heatmap_multi.cpp
)
set(SOURCES_PULSE_EXTRACTOR_COLLECTION
    ${PROJECT_SOURCE_DIR}/src/pulse_extractor_collection.cpp
    ${PROJECT_SOURCE_DIR}/src/wave_processing.cpp
    ${PROJECT_SOURCE_DIR}/src/waveform_utils.cpp
    ${PROJECT_SOURCE_DIR}/src/WaveformManager_2.cpp
    ${PROJECT_SOURCE_DIR}/src/DataLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/heatmap_multi.cpp
)
set(SOURCES_TEMPLATE_MAKER
    ${PROJECT_SOURCE_DIR}/src/template_maker.cpp
    ${PROJECT_SOURCE_DIR}/src/wave_processing.cpp
    ${PROJECT_SOURCE_DIR}/src/waveform_utils.cpp
    ${PROJECT_SOURCE_DIR}/src/WaveformManager_2.cpp  
    ${PROJECT_SOURCE_DIR}/src/DataLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/heatmap_multi.cpp  
)
set(SOURCES_TEMPLATE_MAKER_INDUCTION
    ${PROJECT_SOURCE_DIR}/src/template_maker_induction.cpp
    ${PROJECT_SOURCE_DIR}/src/wave_processing.cpp
    ${PROJECT_SOURCE_DIR}/src/waveform_utils.cpp
    ${PROJECT_SOURCE_DIR}/src/WaveformManager_2.cpp  
    ${PROJECT_SOURCE_DIR}/src/DataLoader.cpp
    ${PROJECT_SOURCE_DIR}/src/heatmap_multi.cpp  
)

# Define the executables with the sources
add_executable(pulse_extractor_induction_1 ${SOURCES_PULSE_EXTRACTOR_INDUCTION_1})
add_executable(pulse_extractor_induction_2 ${SOURCES_PULSE_EXTRACTOR_INDUCTION_2})
add_executable(pulse_extractor_collection ${SOURCES_PULSE_EXTRACTOR_COLLECTION})
add_executable(template_maker ${SOURCES_TEMPLATE_MAKER})
add_executable(template_maker_induction ${SOURCES_TEMPLATE_MAKER_INDUCTION})

# Link the executables with the ROOT libraries
target_link_libraries(pulse_extractor_induction_1 PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(pulse_extractor_induction_2 PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(pulse_extractor_collection PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(template_maker PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(template_maker_induction PUBLIC ${ROOT_LIBRARIES})

