cmake_minimum_required(VERSION 3.12)
project(AdvancedAudioProcessor VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PortAudio REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(JUCE REQUIRED)
find_package(SQLite3 REQUIRED)

# Add executable
add_executable(AdvancedAudioProcessor 
    src/main.cpp
    src/AudioProcessor.cpp
)

# Include directories
target_include_directories(AdvancedAudioProcessor PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PORTAUDIO_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
    ${JUCE_INCLUDE_DIRS}
    ${SQLITE3_INCLUDE_DIRS}
)

# Link libraries
target_link_libraries(AdvancedAudioProcessor PRIVATE 
    ${PORTAUDIO_LIBRARIES}
    ${FFTW3_LIBRARIES}
    ${JUCE_LIBRARIES}
    ${SQLITE3_LIBRARIES}
)
