#pragma once

#include <portaudio.h>
#include <fftw3.h>
#include <vector>
#include <complex>

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

    void performFFT();
    void performIFFT();

    PaStream *stream;
    PaStreamParameters inputParameters;
    PaStreamParameters outputParameters;
    double sampleRate;
    unsigned long framesPerBuffer;

    std::vector<float> inputBuffer;
    std::vector<float> outputBuffer;

    // FFT-related members
    fftwf_plan forwardPlan;
    fftwf_plan inversePlan;
    std::vector<float> fftInput;
    std::vector<std::complex<float>> fftOutput;
    std::vector<std::complex<float>> ifftInput;
    std::vector<float> ifftOutput;
};
