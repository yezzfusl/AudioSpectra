cmake_minimum_required(VERSION 3.12)
project(AdvancedAudioProcessor VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PortAudio REQUIRED)
find_package(FFTW3 REQUIRED)

# JUCE setup
add_subdirectory(JUCE)
juce_add_module(juce_dsp)

# ImGui setup
add_subdirectory(imgui)

# Add executable
add_executable(AdvancedAudioProcessor 
    src/main.cpp
    src/AudioProcessor.cpp
    src/UserInterface.cpp
)

# Include directories
target_include_directories(AdvancedAudioProcessor PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PORTAUDIO_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/backends
)

# Link libraries
target_link_libraries(AdvancedAudioProcessor PRIVATE 
    ${PORTAUDIO_LIBRARIES}
    ${FFTW3_LIBRARIES}
    juce::juce_dsp
    imgui
    glfw
    OpenGL::GL
)

# JUCE-specific settings
target_compile_definitions(AdvancedAudioProcessor PRIVATE
    JUCE_GLOBAL_MODULE_SETTINGS_INCLUDED=1
    JUCE_STANDALONE_APPLICATION=1
)
