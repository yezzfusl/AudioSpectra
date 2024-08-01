#pragma once

#include <portaudio.h>
#include <vector>

class AudioProcessor {
public:
    AudioProcessor();
    ~AudioProcessor();

    bool initialize(int inputDevice = -1, int outputDevice = -1, 
                    double sampleRate = 44100.0, unsigned long framesPerBuffer = 256);
    bool start();
    bool stop();

private:
    static int audioCallback(const void *inputBuffer, void *outputBuffer,
                             unsigned long framesPerBuffer,
                             const PaStreamCallbackTimeInfo* timeInfo,
                             PaStreamCallbackFlags statusFlags,
                             void *userData);

    PaStream *stream;
    PaStreamParameters inputParameters;
    PaStreamParameters outputParameters;
    double sampleRate;
    unsigned long framesPerBuffer;

    std::vector<float> inputBuffer;
    std::vector<float> outputBuffer;
};
