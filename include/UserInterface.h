#pragma once

#include "AudioProcessor.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>

class UserInterface {
public:
    UserInterface(AudioProcessor& processor);
    ~UserInterface();

    bool initialize();
    void render();
    bool shouldClose();

private:
    AudioProcessor& audioProcessor;
    GLFWwindow* window;

    // EQ parameters
    float lowFreq = 100.0f;
    float midFreq = 1000.0f;
    float highFreq = 5000.0f;
    float lowGain = 0.0f;
    float midGain = 0.0f;
    float highGain = 0.0f;

    // Reverb parameters
    float reverbRoomSize = 0.5f;
    float reverbDamping = 0.5f;
    float reverbWetLevel = 0.33f;
    float reverbDryLevel = 0.4f;

    // Compressor parameters
    float compressorThreshold = -10.0f;
    float compressorRatio = 2.0f;
    float compressorAttack = 5.0f;
    float compressorRelease = 100.0f;

    void renderEQControls();
    void renderReverbControls();
    void renderCompressorControls();
};
