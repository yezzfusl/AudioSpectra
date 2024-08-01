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
    
    std::cout << "Audio processor started with basic EQ." << std::endl;
    std::cout << "Adjusting EQ parameters..." << std::endl;

    // Adjust EQ parameters
    processor.setLowGain(3.0f);    // Boost low frequencies by 3 dB
    processor.setMidGain(-2.0f);   // Cut mid frequencies by 2 dB
    processor.setHighGain(4.0f);   // Boost high frequencies by 4 dB

    std::cout << "EQ adjusted. Audio will play for 10 seconds." << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(10));
    
    processor.stop();
    
    return 0;
}
