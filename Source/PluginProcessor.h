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
#define MAX_BUF 2048
#define DEL_LAG 512
#define AMP_CNST 0.708 // -3 dB
#define NUM_F0 8

//==============================================================================
/**
*/
class VibratoTransferAudioProcessor  : public juce::AudioProcessor,
                                       public juce::ValueTree::Listener
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
    
    juce::AudioProcessorValueTreeState& getVTS() { return parameters; }

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
    float del_buffer_l[V_NFFT]; // this buffers the left input channel
    float del_buffer_r[V_NFFT]; // this buffers the right input channel
    int write_pointer = DEL_LAG; // initial delay offset
    float read_pointer = 0;
    std::vector<Biquad> dt_bandpass;
    void initialize_dt_bp();
    
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
    
    float last_f0s[NUM_F0]; // starting with 8, totally arbitrary
    int last_f0s_mask = NUM_F0-1;
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
    float hilbert_left_buf[MAX_BUF];
    float hilbert_right_buf[MAX_BUF];
    
    // f0 isolation bandpass filter
    std::vector<Biquad> f0_bandpass;
    float filter_f0 = 0.f; // set using SNAC
    bool bp_initialized = false; // set after butter chain is set
    bool process_delay = false;
    float bp_buf[MAX_BUF];
    
    // delay function business
    float last_phase = 0;
    float slowest_vibrato = 1; // Hz TODO: parameterize?
    int averaging_frames; // set in prepareToPlay using fs
    float previous_f0_sum = 0;
    int previous_f0_count = 0; // using these two to calculate previous f0 mean
    float dt_buf[MAX_BUF];
    
    // cumsum of RFS (this is effectively the last offset of the delay function)
    float Dt = 0;
    
    // envelope business
    std::vector<Biquad> envelopeBP;
    float last_env = AMP_CNST;
    void initialize_env_bp();
    float env_buf[MAX_BUF];
    float env_cf = 0.f;
    float env_cf_inc = 0.f; // set in prepareToPlay
    
    // onset business TODO: figure most of this out
    float onset_level = 0.1; // -20 dB
    int onset_time_blocks; // set in prepareToPlay using fs
    int blocks_processed = 0;
        
    // functions I added
    void fractional_delay_read(float index, float& left, float& right);
    float find_f0_SNAC();
    void initialize_bp(float bp_low, float bp_high);
    bool bufferTooQuiet(auto* data, int size);
    bool f0Stable();
    void computeAC();

    // class to help us design butterworth filters
    Butterworth b_designer;
    
    // visualizing delay
    VibVisualizer del_vis;
    VibVisualizer amp_vis;
    juce::AudioBuffer<float> dt_buffer;
    juce::AudioBuffer<float> amp_buffer;
    
    // parameters
    void valueTreePropertyChanged (juce::ValueTree &treeWhosePropertyHasChanged, const juce::Identifier &property) override;
    juce::AudioProcessorValueTreeState parameters;
    std::atomic<float>* fmScaler;
    std::atomic<float>* amScaler;
    std::atomic<float>* makeUpGain;
    // these are the actual float values
    float dt_scaler = 1.f;
    float amp_scaler = 1.f;
    float make_up_gain = 1.f;
    std::atomic<bool> updateParams { false };
    
    void resetProcessing();
    
    // private arrays needed for mayer FFT
    float coswrk[20]=
        {
         .00000000000000000000000000000000000000000000000000,
         .70710678118654752440084436210484903928483593768847,
         .92387953251128675612818318939678828682241662586364,
         .98078528040323044912618223613423903697393373089333,
         .99518472667219688624483695310947992157547486872985,
         .99879545620517239271477160475910069444320361470461,
         .99969881869620422011576564966617219685006108125772,
         .99992470183914454092164649119638322435060646880221,
         .99998117528260114265699043772856771617391725094433,
         .99999529380957617151158012570011989955298763362218,
         .99999882345170190992902571017152601904826792288976,
         .99999970586288221916022821773876567711626389934930,
         .99999992646571785114473148070738785694820115568892,
         .99999998161642929380834691540290971450507605124278,
         .99999999540410731289097193313960614895889430318945,
         .99999999885102682756267330779455410840053741619428
        };
    float sinwrk[20]=
        {
         1.0000000000000000000000000000000000000000000000000,
         .70710678118654752440084436210484903928483593768846,
         .38268343236508977172845998403039886676134456248561,
         .19509032201612826784828486847702224092769161775195,
         .09801714032956060199419556388864184586113667316749,
         .04906767432741801425495497694268265831474536302574,
         .02454122852291228803173452945928292506546611923944,
         .01227153828571992607940826195100321214037231959176,
         .00613588464915447535964023459037258091705788631738,
         .00306795676296597627014536549091984251894461021344,
         .00153398018628476561230369715026407907995486457522,
         .00076699031874270452693856835794857664314091945205,
         .00038349518757139558907246168118138126339502603495,
         .00019174759731070330743990956198900093346887403385,
         .00009587379909597734587051721097647635118706561284,
         .00004793689960306688454900399049465887274686668768
        };
};
