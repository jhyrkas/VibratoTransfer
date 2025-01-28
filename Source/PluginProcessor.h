/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <math.h>
#include <numbers>
#include <vector>
#include <JuceHeader.h>
#include "Biquad.hpp"
#include "Butterworth.h"
#include "d_fft_mayer.h"
#include "VibVisualizer.h"

#define V_NFFT 4096
#define V_H_NFFT 2048
#define DEL_LAG 512

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
    
    VibVisualizer& getDelayVisualizer();
    VibVisualizer& getAmpVisualizer();

    // these should probably be set using public methods instead of being public values
    float dt_scaler = 1.f;
    float amp_scaler = 1.f;
private:
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VibratoTransferAudioProcessor)
    // basic constants
    float fs; // set in prepareToPlay
    float T; // 1/fs, set in prepareToPlay
    int blockSize; // set in prepareToPlay
    float twopi = 2 * std::numbers::pi_v<float>;
    
    // delay business
    // TODO: use DEFINE for del_length?
    int del_length = V_NFFT; // TODO: think about if this should be longer/shorter
    int del_length_mask = V_NFFT - 1;
    float del_buffer[V_NFFT]; // this buffers the input signal
    int write_pointer = DEL_LAG; // initial delay offset
    float read_pointer = 0;
    
    // f0 analysis business
    // 2*NFFT because we store the FFT result in place
    float sc_buffer[V_H_NFFT]; // this buffers the sidechain / analysis signal
    float ac_buffer[V_NFFT]; // because we need the f0 signal and the ac signal to compute the norm
    int Nfft_mask = V_H_NFFT - 1;
    int sc_pointer = 0; // used to index into the f0 signal buffer
    int n_buffered = 0; // buffered samples after onset
    int snac_end_index = 0; // set in prepareToPlay
    
    // figuring out when to start and stop the delay TODO: figure this out
    bool f0_stabilized = false;
    bool filter_initialized = false;
    
    //float last_f0s[4] = {0.f, 0.f, 0.f, 0.f}; // starting with 4, totally arbitrary
    //int last_f0s_mask = 3;
    
    float last_f0s[8] = {0.f, 0.f, 0.f, 0.f,0.f, 0.f, 0.f, 0.f}; // starting with 8, totally arbitrary
    int last_f0s_mask = 7;
    
    int last_f0s_pointer = 0;
    
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
    
    // f0 isolation bandpass filter
    std::vector<Biquad> f0_bandpass;
    float filter_f0 = 0.f; // set using SNAC
    bool bp_initialized = false; // set after butter chain is set
    bool process_delay = false;
    
    // delay function business
    float last_phase = 0;
    float slowest_vibrato = 1; // Hz TODO: parameterize?
    int averaging_frames; // set in prepareToPlay using fs
    float previous_f0_sum = 0;
    int previous_f0_count = 0; // using these two to calculate previous f0 mean
    // cumsum of RFS (this is effectively the last offset of the delay function)
    //float Dt = 0;
    
    // envelope business
    std::vector<Biquad> envelopeBP;
    float last_env = 1.f;
    void initialize_env_bp(double sampleRate);
    
    // onset business TODO: figure most of this out
    float onset_level = 0.1; // -20 dB
    int onset_time_blocks; // set in prepareToPlay using fs
    int blocks_processed = 0;
        
    // functions I added
    float fractional_delay_read(float index);
    float find_f0_SNAC();
    void initialize_bp(float bp_low, float bp_high);
    bool bufferTooQuiet(auto* data, int size);
    bool f0Stable();

    // class to help us design butterworth filters
    Butterworth bp_designer;
    
    // visualizing delay
    VibVisualizer del_vis;
    VibVisualizer amp_vis;
    float dt_buffer[V_NFFT];
    float amp_buffer[V_NFFT];
    int vis_buffer_ptr = 0;
    int vis_ptr_mask = V_NFFT - 1; // redudant but it keeps the code cleaner
};
