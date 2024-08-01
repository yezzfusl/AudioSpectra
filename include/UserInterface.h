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

    float lowFreq = 100.0f;
    float midFreq = 1000.0f;
    float highFreq = 5000.0f;
    float lowGain = 0.0f;
    float midGain = 0.0f;
    float highGain = 0.0f;

    void renderEQControls();
};
