/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
VibratoTransferAudioProcessor::VibratoTransferAudioProcessor()
:
#ifndef JucePlugin_PreferredChannelConfigurations
     AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       ),
#endif
peakingFilter(0.f, 0.f, 0.f, 0.f, 0.f)
{
    // zero out delay buffer
    memset(del_buffer, 0, del_length * sizeof(float));
    memset(f0_buffer, 0, V_H_NFFT * sizeof(float));
    memset(ac_buffer, 0, V_NFFT*sizeof(float));
}

VibratoTransferAudioProcessor::~VibratoTransferAudioProcessor()
{
}

//==============================================================================
const juce::String VibratoTransferAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool VibratoTransferAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool VibratoTransferAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool VibratoTransferAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double VibratoTransferAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int VibratoTransferAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int VibratoTransferAudioProcessor::getCurrentProgram()
{
    return 0;
}

void VibratoTransferAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String VibratoTransferAudioProcessor::getProgramName (int index)
{
    return {};
}

void VibratoTransferAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void VibratoTransferAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    fs = (float) sampleRate;
    T = 1.f / fs;
    blockSize = samplesPerBlock;
    averaging_frames = int((fs / slowest_vibrato) / samplesPerBlock);
    // 100 ms - parameterize?
    onset_time_blocks = int(0.1 * fs/samplesPerBlock);
    snac_end_index = int(fs/50); // min f0 is 50 Hz
    // TODO: memset delay and f0 bufs to 0?
}

void VibratoTransferAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool VibratoTransferAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

// TODO: right now we're just assuming everything is mono but the logic should
// TODO: be redone a bit if we handle stereo or beyond.
// TODO: also, one of these channels should be the sidechain?
void VibratoTransferAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.
    
    // TODO: to get things started off, let's just start by buffering and make
    // TODO: sure that we don't crash
    int totalChans = totalNumInputChannels > 0 ? 1 : 0;
    for (int channel = 0; channel < totalChans; ++channel)
    {
        auto* channelData = buffer.getWritePointer (channel);
        
        for (int i = 0; i < blockSize; ++i) {
            // STEP 1: buffer f0 signal (eventually from the sidechain)
            f0_buffer[f0_pointer] = channelData[i];
            ac_buffer[f0_pointer] = channelData[i];
            f0_pointer = (f0_pointer + 1) & Nfft_mask;
            ++n_buffered;
            
            // STEP 1.1: if we have buffered enough, do f0 analysis
            // TODO: is there any way to flag this outside of the for-loop?
            if (n_buffered == V_H_NFFT) {
                n_buffered = 0; // reset
                
                // autocorrelation and SNAC
                mayer_realfft(V_NFFT,ac_buffer);
                ac_buffer[0] *= ac_buffer[0]; // DC
                ac_buffer[V_H_NFFT] *= ac_buffer[V_H_NFFT]; // Nyquist
                for (int k = 1; k < V_H_NFFT; ++k) {
                    ac_buffer[k] = (ac_buffer[k] * ac_buffer[k] + ac_buffer[V_NFFT-k] * ac_buffer[V_NFFT-k]); // real
                    ac_buffer[V_NFFT-k] = 0.0; // imaginary
                }
                mayer_realifft(V_NFFT, ac_buffer);
                
                // peak finding
                float f0 = find_f0_SNAC();
                last_f0s[last_f0s_pointer] = f0;
                last_f0s_pointer = (last_f0s_pointer+1) & last_f0s_mask;
                
                // STEP 1.2 TODO: if f0s have stabilized and the filter isn't initialized, initialize it
                // if (something)
                //        initialize_bp(0.95* f0, 1.05 * f0)
                
                // now reset the buffer
                memset(ac_buffer, 0, V_NFFT * sizeof(float));
                // TODO: reset f0 buffer? if the logic all holds then we shouldn't ever need to reset it
            }
            
            // STEP 2: buffer delay signal (buffer the input)
            del_buffer[write_pointer] = channelData[i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            
            // STEP 3: read from the delay line
            // read (right now do nothing...there should be a 512 sample delay)
            channelData[i] = fractional_delay_read(read_pointer);
            
            // STEP 4: calculate the next delay based on the output of the bandpass filter and f0 analysis
            // if (something)
            //      bandpass
            //      hilbert
            //      angle -> unwrap
            //      diff for if
            //      delay offset = 1 - (if/w0)
            
            // PLACEHOLDER doing a regular vibrato just to make sure that read pointer is working
            read_pointer = (read_pointer + 1) - (0.1f * sinf(sin_phase));
            read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            sin_phase += (twopi * T * 1);
            sin_phase = sin_phase > twopi ? sin_phase - twopi : sin_phase;
        }
    }
}

//==============================================================================
bool VibratoTransferAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* VibratoTransferAudioProcessor::createEditor()
{
    return new VibratoTransferAudioProcessorEditor (*this);
}

//==============================================================================
void VibratoTransferAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void VibratoTransferAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new VibratoTransferAudioProcessor();
}

// functions I added
float VibratoTransferAudioProcessor::fractional_delay_read(float index) {
    // TODO: should we ever check that index is in range, or should that be done
    // TODO: outside of the function?
    int low = int(index);
    int high = low == del_length - 1 ? 0 : low+1;
    float high_frac = index - low;
    return (1-high_frac) * del_buffer[low] + high_frac * del_buffer[high];
}

float VibratoTransferAudioProcessor::find_f0_SNAC() {
    
    // TODO: any way to avoid divisions here?
    double norm = 2 * ac_buffer[0];
    //double norm = 2 * ac_buffer[0]; // website recommends this be a double
    double decr = (1.0 / (V_H_NFFT-1))*0.2; // triangular window from 1 to 0.8
    ac_buffer[0] = (2 * ac_buffer[0]) / norm;
    for (int i = 1; i < V_H_NFFT; ++i) {
        norm = norm - (f0_buffer[i-1]*f0_buffer[i-1] + f0_buffer[V_H_NFFT-i]*f0_buffer[V_H_NFFT-i]);
        ac_buffer[i] = (1-i*decr)*(2 * ac_buffer[i]) / norm;
    }
    
    int peak_index = 1;
    float peak_max = 0.f;
    for (int i = 2; i < snac_end_index-1; ++i) {
        if (ac_buffer[i-1] < ac_buffer[i] && ac_buffer[i+1] < ac_buffer[i] && ac_buffer[i] > peak_max) {
            peak_index = i;
            peak_max = ac_buffer[i];
        }
    }
    
    return fs / peak_index;
}

// doing this now for the peaking filter, change later for butter chain
void VibratoTransferAudioProcessor::initialize_bp(float bp_low, float bp_high) {
    // more or less copied from the audio bashing
    float f0 = (bp_low + bp_high) * 0.5;
    float bw = (bp_high - bp_low);
    float w0 = (twopi * f0) / fs;
    float bwr = (twopi * bw) / fs;
    float gain = 1.f / sqrt(2); // don't really need to compute every time
    float beta = (sqrt(1.f - (gain*gain))/gain) * tanf(bwr*0.5f);
    float norm_gain = 1.f/(1+beta);
    float ng_cos_w0_neg2 = -2.f*norm_gain*cosf(w0);
    peakingFilter.setParams(1.f-norm_gain, 0.f, norm_gain-1.f, ng_cos_w0_neg2, 2.f*norm_gain-1.f);
    bp_initialized = true;
}
