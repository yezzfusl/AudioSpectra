#include "../include/AudioProcessor.h"
#include "../include/UserInterface.h"
#include <iostream>

int main() {
    std::cout << "Advanced Audio Processor" << std::endl;
    
    AudioProcessor processor;
    
    if (!processor.initialize()) {
        std::cerr << "Failed to initialize audio processor" << std::endl;
        return 1;
    }
    
    if (!processor.start()) {
        std::cerr << "Failed to start audio processor" << std::endl;
        return 1;
    }
    
    UserInterface ui(processor);
    if (!ui.initialize()) {
        std::cerr << "Failed to initialize user interface" << std::endl;
        processor.stop();
        return 1;
    }

    std::cout << "Audio processor and UI started. Adjust EQ parameters in the window." << std::endl;

    while (!ui.shouldClose()) {
        ui.render();
    }
    
    processor.stop();
    
    return 0;
}
