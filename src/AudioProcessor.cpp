#include "../include/AudioProcessor.h"
#include <iostream>
#include <cstring>

AudioProcessor::AudioProcessor() : stream(nullptr), sampleRate(44100.0), framesPerBuffer(256) {}

AudioProcessor::~AudioProcessor() {
    stop();
    Pa_Terminate();
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

int AudioProcessor::audioCallback(const void *inputBuffer, void *outputBuffer,
                                  unsigned long framesPerBuffer,
                                  const PaStreamCallbackTimeInfo* timeInfo,
                                  PaStreamCallbackFlags statusFlags,
                                  void *userData) {
    AudioProcessor* processor = static_cast<AudioProcessor*>(userData);
    float* in = (float*)inputBuffer;
    float* out = (float*)outputBuffer;

    // Simple pass-through for now
    for (unsigned long i = 0; i < framesPerBuffer * 2; i++) {
        out[i] = in[i];
    }

    return paContinue;
}
