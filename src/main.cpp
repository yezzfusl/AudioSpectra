#include "../include/AudioProcessor.h"
#include <iostream>
#include <thread>
#include <chrono>

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
    
    std::cout << "Audio processor started. Press Enter to stop..." << std::endl;
    std::cin.get();
    
    processor.stop();
    
    return 0;
}
