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
fft(11) // order is the exponent (i.e. 2^11 = 2048)
{
    // zero out delay buffer
    memset(del_buffer, 0, del_length * sizeof(float));
    memset(f0_buffer, 0, Nfft * sizeof(float));
    memset(ac_buffer, 0, 2 * Nfft*sizeof(float));
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
            // buffer f0 signal
            
            f0_buffer[f0_pointer] = channelData[i];
            ac_buffer[f0_pointer] = channelData[i];
            f0_pointer = (f0_pointer + 1) & Nfft_mask;
            // i would like to check n_buffered at the block level, not the
            // loop level...is it guaranteed that blockSize is a factor of Nfft?
            // probably not, so we should think about that
            ++n_buffered;
            
            if (n_buffered == Nfft) {
                n_buffered = 0;
                fft.performFrequencyOnlyForwardTransform(ac_buffer);
                // autocorrelation
                for (int j = 0; j < Nfft; ++j) {
                    ac_buffer[2*j] = ac_buffer[j]*ac_buffer[j] + (ac_buffer[j+1] * ac_buffer[j+1]);
                    ac_buffer[2*j+1] = 0.f;
                }
                fft.performRealOnlyInverseTransform(ac_buffer);
                // peak finding
                float f0 = find_f0_SNAC();
                // now reset the buffer
                memset(ac_buffer, 0, 2*Nfft * sizeof(float));
                // TODO: reset f0 buffer? if the logic all holds then we shouldn't ever need to reset it
            }
            
            // buffer delay signal
            del_buffer[write_pointer] = channelData[i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            
            // read (right now do nothing...there should be a 512 sample delay)
            channelData[i] = fractional_delay_read(read_pointer);
            
            // trying to do a regular vibrato, make sure that read pointer is always in range
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
    double norm = 2 * f0_buffer[0]; // website recommends this be a double
    double decr = (1.0 / (Nfft-1))*0.2; // triangular window from 1 to 0.8
    ac_buffer[0] = (2 * ac_buffer[0]) / norm;
    for (int i = 1; i < Nfft; ++i) {
        norm = norm - (f0_buffer[i-1]*f0_buffer[i-1] + f0_buffer[Nfft-i]*f0_buffer[Nfft-i]);
        ac_buffer[i] = (i*decr)*(2 * ac_buffer[i]) / norm;
    }
    
    int peak_index = 1;
    for (int i = 2; i < snac_end_index-1; ++i) {
        if (ac_buffer[i-1] < ac_buffer[i] && ac_buffer[i+1] < ac_buffer[i] && ac_buffer[i] > ac_buffer[peak_index]) {
            peak_index = i;
        }
    }
    
    return fs / peak_index;
}
