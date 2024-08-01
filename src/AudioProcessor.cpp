#include "../include/AudioProcessor.h"
#include <iostream>
#include <cstring>
#include <algorithm>

AudioProcessor::AudioProcessor() : stream(nullptr), sampleRate(44100.0), framesPerBuffer(256),
                                   forwardPlan(nullptr), inversePlan(nullptr) {}

AudioProcessor::~AudioProcessor() {
    stop();
    Pa_Terminate();
    if (forwardPlan) fftwf_destroy_plan(forwardPlan);
    if (inversePlan) fftwf_destroy_plan(inversePlan);
    fftwf_cleanup();
}

bool AudioProcessor::initialize(int inputDevice, int outputDevice, double sampleRate, unsigned long framesPerBuffer) {
    PaError err;

    err = Pa_Initialize();
    if (err != paNoError) {
        std::cerr << "PortAudio error: " << Pa_GetErrorText(err) << std::endl;
        return false;
    }

    this->sampleRate = sampleRate;
    this->framesPerBuffer = framesPerBuffer;

    inputParameters.device = inputDevice == -1 ? Pa_GetDefaultInputDevice() : inputDevice;
    inputParameters.channelCount = 2;
    inputParameters.sampleFormat = paFloat32;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = nullptr;

    outputParameters.device = outputDevice == -1 ? Pa_GetDefaultOutputDevice() : outputDevice;
    outputParameters.channelCount = 2;
    outputParameters.sampleFormat = paFloat32;
    outputParameters.suggestedLatency = Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = nullptr;

    err = Pa_OpenStream(&stream,
                        &inputParameters,
                        &outputParameters,
                        sampleRate,
                        framesPerBuffer,
                        paClipOff,
                        audioCallback,
                        this);

    if (err != paNoError) {
        std::cerr << "PortAudio error: " << Pa_GetErrorText(err) << std::endl;
        return false;
    }

    inputBuffer.resize(framesPerBuffer * 2);
    outputBuffer.resize(framesPerBuffer * 2);

    // Initialize FFT-related buffers and plans
    fftInput.resize(framesPerBuffer);
    fftOutput.resize(framesPerBuffer / 2 + 1);
    ifftInput.resize(framesPerBuffer / 2 + 1);
    ifftOutput.resize(framesPerBuffer);

    forwardPlan = fftwf_plan_dft_r2c_1d(framesPerBuffer, fftInput.data(), 
                                        reinterpret_cast<fftwf_complex*>(fftOutput.data()), FFTW_MEASURE);
    inversePlan = fftwf_plan_dft_c2r_1d(framesPerBuffer, 
                                        reinterpret_cast<fftwf_complex*>(ifftInput.data()), 
                                        ifftOutput.data(), FFTW_MEASURE);

    return true;
}

bool AudioProcessor::start() {
    PaError err = Pa_StartStream(stream);
    if (err != paNoError) {
        std::cerr << "PortAudio error: " << Pa_GetErrorText(err) << std::endl;
        return false;
    }
    return true;
}

bool AudioProcessor::stop() {
    if (stream) {
        PaError err = Pa_StopStream(stream);
        if (err != paNoError) {
            std::cerr << "PortAudio error: " << Pa_GetErrorText(err) << std::endl;
            return false;
        }
        err = Pa_CloseStream(stream);
        if (err != paNoError) {
            std::cerr << "PortAudio error: " << Pa_GetErrorText(err) << std::endl;
            return false;
        }
        stream = nullptr;
    }
    return true;
}

void AudioProcessor::performFFT() {
    fftwf_execute(forwardPlan);
}

void AudioProcessor::performIFFT() {
    fftwf_execute(inversePlan);
    // Normalize the output
    float normFactor = 1.0f / framesPerBuffer;
    for (auto& sample : ifftOutput) {
        sample *= normFactor;
    }
}

int AudioProcessor::audioCallback(const void *inputBuffer, void *outputBuffer,
                                  unsigned long framesPerBuffer,
                                  const PaStreamCallbackTimeInfo* timeInfo,
                                  PaStreamCallbackFlags statusFlags,
                                  void *userData) {
    AudioProcessor* processor = static_cast<AudioProcessor*>(userData);
    float* in = (float*)inputBuffer;
    float* out = (float*)outputBuffer;

    // Process left channel (assuming stereo input)
    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        processor->fftInput[i] = in[i * 2];
    }

    processor->performFFT();

    // Simple spectral processing: attenuate high frequencies
    for (unsigned long i = framesPerBuffer / 4; i < processor->fftOutput.size(); i++) {
        processor->fftOutput[i] *= 0.5f;
    }

    processor->ifftInput = processor->fftOutput;
    processor->performIFFT();

    // Copy processed left channel to both output channels
    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        out[i * 2] = out[i * 2 + 1] = processor->ifftOutput[i];
    }

    return paContinue;
}
