#pragma once

#include <portaudio.h>
#include <fftw3.h>
#include <vector>
#include <complex>
#include <juce_dsp/juce_dsp.h>

class AudioProcessor {
public:
    AudioProcessor();
    ~AudioProcessor();

    bool initialize(int inputDevice = -1, int outputDevice = -1, 
                    double sampleRate = 44100.0, unsigned long framesPerBuffer = 256);
    bool start();
    bool stop();

    // EQ methods
    void setLowFrequency(float freq);
    void setMidFrequency(float freq);
    void setHighFrequency(float freq);
    void setLowGain(float gain);
    void setMidGain(float gain);
    void setHighGain(float gain);

    // Reverb methods
    void setReverbRoomSize(float size);
    void setReverbDamping(float damping);
    void setReverbWetLevel(float wetLevel);
    void setReverbDryLevel(float dryLevel);

    // Compressor methods
    void setCompressorThreshold(float threshold);
    void setCompressorRatio(float ratio);
    void setCompressorAttack(float attack);
    void setCompressorRelease(float release);

private:
    static int audioCallback(const void *inputBuffer, void *outputBuffer,
                             unsigned long framesPerBuffer,
                             const PaStreamCallbackTimeInfo* timeInfo,
                             PaStreamCallbackFlags statusFlags,
                             void *userData);

    void performFFT();
    void performIFFT();
    void processAudio(const float* inputBuffer, float* outputBuffer, unsigned long framesPerBuffer);

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

    // JUCE DSP-related members
    juce::dsp::ProcessorChain<juce::dsp::IIR::Filter<float>, 
                              juce::dsp::IIR::Filter<float>, 
                              juce::dsp::IIR::Filter<float>,
                              juce::dsp::Reverb,
                              juce::dsp::Compressor<float>> processorChain;
    juce::dsp::ProcessSpec processSpec;
};
