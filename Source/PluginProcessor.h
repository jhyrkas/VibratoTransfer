/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <numbers>
#include <JuceHeader.h>
#include "Biquad.hpp"

//==============================================================================
/**
*/
class VibratoTransferAudioProcessor  : public juce::AudioProcessor
                            #if JucePlugin_Enable_ARA
                             , public juce::AudioProcessorARAExtension
                            #endif
{
public:
    //==============================================================================
    VibratoTransferAudioProcessor();
    ~VibratoTransferAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

private:
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VibratoTransferAudioProcessor)
    // basic constants
    float fs; // set in prepareToPlay
    float T; // 1/fs, set in prepareToPlay
    int blockSize; // set in prepareToPlay
    float twopi = 2 * std::numbers::pi_v<float>;
    
    // delay business
    int del_length = 4096; // TODO: think about if this should be longer/shorter
    int del_length_mask = 4095;
    float del_buffer[4096]; // TODO: use DEFINE for del_length?
    int write_pointer = 512; // initial delay offset
    int read_pointer = 0;
    
    // f0 business
    int Nfft = 2048; // for pitch detection
    int Nfft_mask = 2047;
    // 2*NFFT because we store the FFT result in place
    float f0_buffer[2*2048]; // TODO: use DEFINE for Nfft?
    int f0_pointer = 0; // used to index into the f0 signal buffer
    int n_buffered = 0; // buffered samples after onset
    bool analyze_f0 = false; // set to true after onset TODO: figure this out
    // TODO: this fft seems to not be initialized correctly or something of that nature?
    //juce::dsp::FFT fft;
    
    // hilbert business
    Biquad hilbert_left[4] = {
        Biquad(0.161758f, 0.f, -1.f, 0.f, -0.161758f),
        Biquad(0.733029f, 0.f, -1.f, 0.f, -0.733029f),
        Biquad(0.945350f, 0.f, -1.f, 0.f, -0.945350f),
        Biquad(0.990598f, 0.f, -1.f, 0.f, -0.990598f)
    };
    Biquad hilbert_right[4] = {
        Biquad(0.479401f, 0.f, -1.f, 0.f, -0.479401f),
        Biquad(0.876218f, 0.f, -1.f, 0.f, -0.876218f),
        Biquad(0.976599f, 0.f, -1.f, 0.f, -0.976599f),
        Biquad(0.997500f, 0.f, -1.f, 0.f, -0.997500f)
    };
    float last_right_out = 0.f;
    
    // biquad business
    // TODO: figure out butter coefs and also number of biquads
    Biquad butter_chain[2] = {
        Biquad(0.f, 0.f, 0.f, 0.f, 0.f), // set these using setParams()
        Biquad(0.f, 0.f, 0.f, 0.f, 0.f) // set these using setParams()
    };
    float filter_f0 = 0.f; // set using SNAC
    bool bp_initialized = false; // set after butter chain is set
    
    // delay function business
    float last_phase = 0;
    float slowest_vibrato = 1; // Hz TODO: parameterize?
    int averaging_frames; // set in prepareToPlay using fs
    float previous_f0_sum = 0;
    float previous_f0_count = 0; // using these two to calculate previous f0 mean
    // cumsum of RFS (this is effectively the last offset of the delay function)
    float Dt = 0;
    
    // onset business TODO: figure most of this out
    float onset_level = 0.1; // -20 dB
    int onset_time_blocks; // set in prepareToPlay using fs
    bool in_vibrato = false; // set once we are analyzing vibrato
    int vibrato_blocks_processed = 0; // increment as we process, set to 0 after offset
    
    // functions I added
    float fractional_delay_read(float index);
};
