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

    // Initialize JUCE DSP
    processSpec.sampleRate = sampleRate;
    processSpec.maximumBlockSize = framesPerBuffer;
    processSpec.numChannels = 1;

    processorChain.prepare(processSpec);

    // Set up initial EQ parameters
    setLowFrequency(100.0f);
    setMidFrequency(1000.0f);
    setHighFrequency(5000.0f);
    setLowGain(0.0f);
    setMidGain(0.0f);
    setHighGain(0.0f);

    // Set up initial reverb parameters
    setReverbRoomSize(0.5f);
    setReverbDamping(0.5f);
    setReverbWetLevel(0.33f);
    setReverbDryLevel(0.4f);

    // Set up initial compressor parameters
    setCompressorThreshold(-10.0f);
    setCompressorRatio(2.0f);
    setCompressorAttack(5.0f);
    setCompressorRelease(100.0f);

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

void AudioProcessor::processAudio(const float* inputBuffer, float* outputBuffer, unsigned long framesPerBuffer) {
    juce::dsp::AudioBlock<float> block(&outputBuffer[0], 1, framesPerBuffer);
    juce::dsp::ProcessContextReplacing<float> context(block);
    processorChain.process(context);
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

    // Apply processing chain
    processor->processAudio(processor->fftInput.data(), processor->ifftOutput.data(), framesPerBuffer);

    processor->performIFFT();

    // Copy processed left channel to both output channels
    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        out[i * 2] = out[i * 2 + 1] = processor->ifftOutput[i];
    }

    return paContinue;
}

void AudioProcessor::setLowFrequency(float freq) {
    auto& filter = processorChain.get<0>();
    *filter.state = *juce::dsp::IIR::Coefficients<float>::makeLowShelf(processSpec.sampleRate, freq, 0.707f, 1.0f);
}

void AudioProcessor::setMidFrequency(float freq) {
    auto& filter = processorChain.get<1>();
    *filter.state = *juce::dsp::IIR::Coefficients<float>::makePeakFilter(processSpec.sampleRate, freq, 0.707f, 1.0f);
}

void AudioProcessor::setHighFrequency(float freq) {
    auto& filter = processorChain.get<2>();
    *filter.state = *juce::dsp::IIR::Coefficients<float>::makeHighShelf(processSpec.sampleRate, freq, 0.707f, 1.0f);
}

void AudioProcessor::setLowGain(float gain) {
    auto& filter = processorChain.get<0>();
    auto* lowShelf = dynamic_cast<juce::dsp::IIR::Coefficients<float>*>(filter.state);
    *lowShelf = *juce::dsp::IIR::Coefficients<float>::makeLowShelf(processSpec.sampleRate, lowShelf->getRawCoefficients()[4], 0.707f, std::pow(10.0f, gain / 20.0f));
}

void AudioProcessor::setMidGain(float gain) {
    auto& filter = processorChain.get<1>();
    auto* peak = dynamic_cast<juce::dsp::IIR::Coefficients<float>*>(filter.state);
    *peak = *juce::dsp::IIR::Coefficients<float>::makePeakFilter(processSpec.sampleRate, peak->getRawCoefficients()[4], 0.707f, std::pow(10.0f, gain / 20.0f));
}

void AudioProcessor::setHighGain(float gain) {
    auto& filter = processorChain.get<2>();
    auto* highShelf = dynamic_cast<juce::dsp::IIR::Coefficients<float>*>(filter.state);
    *highShelf = *juce::dsp::IIR::Coefficients<float>::makeHighShelf(processSpec.sampleRate, highShelf->getRawCoefficients()[4], 0.707f, std::pow(10.0f, gain / 20.0f));
}

void AudioProcessor::setReverbRoomSize(float size) {
    auto& reverb = processorChain.get<3>();
    auto params = reverb.getParameters();
    params.roomSize = size;
    reverb.setParameters(params);
}

void AudioProcessor::setReverbDamping(float damping) {
    auto& reverb = processorChain.get<3>();
    auto params = reverb.getParameters();
    params.damping = damping;
    reverb.setParameters(params);
}

void AudioProcessor::setReverbWetLevel(float wetLevel) {
    auto& reverb = processorChain.get<3>();
    auto params = reverb.getParameters();
    params.wetLevel = wetLevel;
    reverb.setParameters(params);
}

void AudioProcessor::setReverbDryLevel(float dryLevel) {
    auto& reverb = processorChain.get<3>();
    auto params = reverb.getParameters();
    params.dryLevel = dryLevel;
    reverb.setParameters(params);
}

void AudioProcessor::setCompressorThreshold(float threshold) {
    auto& compressor = processorChain.get<4>();
    compressor.setThreshold(threshold);
}

void AudioProcessor::setCompressorRatio(float ratio) {
    auto& compressor = processorChain.get<4>();
    compressor.setRatio(ratio);
}

void AudioProcessor::setCompressorAttack(float attack) {
    auto& compressor = processorChain.get<4>();
    compressor.setAttack(attack);
}

void AudioProcessor::setCompressorRelease(float release) {
    auto& compressor = processorChain.get<4>();
    compressor.setRelease(release);
