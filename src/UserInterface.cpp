#include "../include/UserInterface.h"
#include <iostream>

UserInterface::UserInterface(AudioProcessor& processor) : audioProcessor(processor), window(nullptr) {}

UserInterface::~UserInterface() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    if (window) {
        glfwDestroyWindow(window);
    }
    glfwTerminate();
}

bool UserInterface::initialize() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return false;
    }

    window = glfwCreateWindow(800, 600, "Audio Processor UI", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return false;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    return true;
}

void UserInterface::render() {
    glfwPollEvents();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("Audio Processor Controls");
    
    if (ImGui::CollapsingHeader("Equalizer")) {
        renderEQControls();
    }
    
    if (ImGui::CollapsingHeader("Reverb")) {
        renderReverbControls();
    }
    
    if (ImGui::CollapsingHeader("Compressor")) {
        renderCompressorControls();
    }

    ImGui::End();

    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(window);
}

bool UserInterface::shouldClose() {
    return glfwWindowShouldClose(window);
}

void UserInterface::renderEQControls() {
    if (ImGui::SliderFloat("Low Frequency", &lowFreq, 20.0f, 1000.0f, "%.1f Hz")) {
        audioProcessor.setLowFrequency(lowFreq);
    }
    if (ImGui::SliderFloat("Mid Frequency", &midFreq, 200.0f, 5000.0f, "%.1f Hz")) {
        audioProcessor.setMidFrequency(midFreq);
    }
    if (ImGui::SliderFloat("High Frequency", &highFreq, 1000.0f, 20000.0f, "%.1f Hz")) {
        audioProcessor.setHighFrequency(highFreq);
    }
    if (ImGui::SliderFloat("Low Gain", &lowGain, -12.0f, 12.0f, "%.1f dB")) {
        audioProcessor.setLowGain(lowGain);
    }
    if (ImGui::SliderFloat("Mid Gain", &midGain, -12.0f, 12.0f, "%.1f dB")) {
        audioProcessor.setMidGain(midGain);
    }
    if (ImGui::SliderFloat("High Gain", &highGain, -12.0f, 12.0f, "%.1f dB")) {
        audioProcessor.setHighGain(highGain);
    }
}

void UserInterface::renderReverbControls() {
    if (ImGui::SliderFloat("Room Size", &reverbRoomSize, 0.0f, 1.0f)) {
        audioProcessor.setReverbRoomSize(reverbRoomSize);
    }
    if (ImGui::SliderFloat("Damping", &reverbDamping, 0.0f, 1.0f)) {
        audioProcessor.setReverbDamping(reverbDamping);
    }
    if (ImGui::SliderFloat("Wet Level", &reverbWetLevel, 0.0f, 1.0f)) {
        audioProcessor.setReverbWetLevel(reverbWetLevel);
    }
    if (ImGui::SliderFloat("Dry Level", &reverbDryLevel, 0.0f, 1.0f)) {
        audioProcessor.setReverbDryLevel(reverbDryLevel);
    }
}

void UserInterface::renderCompressorControls() {
    if (ImGui::SliderFloat("Threshold", &compressorThreshold, -60.0f, 0.0f, "%.1f dB")) {
        audioProcessor.setCompressorThreshold(compressorThreshold);
    }
    if (ImGui::SliderFloat("Ratio", &compressorRatio, 1.0f, 20.0f, "%.1f:1")) {
        audioProcessor.setCompressorRatio(compressorRatio);
    }
    if (ImGui::SliderFloat("Attack", &compressorAttack, 0.1f, 100.0f, "%.1f ms")) {
        audioProcessor.setCompressorAttack(compressorAttack);
    }
    if (ImGui::SliderFloat("Release", &compressorRelease, 10.0f, 1000.0f, "%.1f ms")) {
        audioProcessor.setCompressorRelease(compressorRelease);
    }
}
