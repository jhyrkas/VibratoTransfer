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
bp_designer(),
del_vis(),
amp_vis(),
dt_buffer(1, MAX_BUF),
amp_buffer(1, MAX_BUF)
//envelopeBP() // ButterBP
//envelopeBP(0, 0, 0, 0, 0) // Biquad
{
    // zero out delay buffer
    memset(del_buffer, 0, del_length * sizeof(float));
    memset(sc_buffer, 0, V_H_NFFT * sizeof(float));
    memset(ac_buffer, 0, V_NFFT*sizeof(float));
    memset(phase_buf, 0, MAX_BUF*sizeof(float));
    memset(env_buf, 0, MAX_BUF*sizeof(float));
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
    
    // TODO: double check, should we clear these here? problably in case fs changes
    for (Biquad quad : f0_bandpass) {
        f0_bandpass.clear();
    }
    hilbert_left[0].clear(); hilbert_left[1].clear(); hilbert_left[2].clear(); hilbert_left[3].clear();
    hilbert_right[0].clear(); hilbert_right[1].clear(); hilbert_right[2].clear(); hilbert_right[3].clear();
    process_delay = false;
    bp_initialized = false;
    //envelopeBP.setParams(1.f, 10.f, fs); // 1 Hz thru 10 Hz, this does a buffer clear
    initialize_env_bp(); // this does a buffer clear
    
    del_vis.clear();
    del_vis.setSamplesPerBlock(samplesPerBlock);
    amp_vis.clear();
    amp_vis.setSamplesPerBlock(samplesPerBlock);
    dt_buffer.setSize(1, blockSize);
    amp_buffer.setSize(1, blockSize);
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
    
    int totalChans = totalNumInputChannels > 0 ? 1 : 0; // TODO: get main and side-chain channels
    
    // analysis loop - side chain (right channel for testing)
    for (int channel = 1; channel < 2; ++channel)
    {
        auto* channelData = buffer.getWritePointer (channel);
        for (int i = 0; i < blockSize; ++i) {
            // STEP 1: buffer f0 signal (eventually from the sidechain)
            sc_buffer[sc_pointer] = channelData[i];
            ac_buffer[sc_pointer] = channelData[i];
            sc_pointer = (sc_pointer + 1) & Nfft_mask;
            ++n_buffered;
            
            // STEP 2: if we have buffered enough, do f0 analysis
            // TODO: is there any way to flag this outside of the for-loop?
            if (n_buffered == V_H_NFFT) {
                n_buffered = 0; // reset
                float f0 = 0.f;
                
                // if analysis signal is too quiet, don't even bother
                if (bufferTooQuiet(ac_buffer, V_H_NFFT)) {
                    if (bp_initialized) {
                        bp_initialized = false;
                    }
                    process_delay = false;
                    last_f0s[last_f0s_pointer] = f0;
                    last_f0s_pointer = (last_f0s_pointer+1) & last_f0s_mask;
                    previous_f0_sum = 0;
                    previous_f0_count = 0;
                    for (Biquad quad : f0_bandpass) {
                        quad.clear();
                    }
                    for (Biquad quad : envelopeBP) {
                        quad.clear();
                    }
                    hilbert_left[0].clear(); hilbert_left[1].clear(); hilbert_left[2].clear(); hilbert_left[3].clear();
                    hilbert_right[0].clear(); hilbert_right[1].clear(); hilbert_right[2].clear(); hilbert_right[3].clear();
                    blocks_processed = 0;
                // do autocorrelation and SNAC
                } else {
                    mayer_realfft(V_NFFT,ac_buffer); // FFT
                    ac_buffer[0] *= ac_buffer[0]; // process DC
                    ac_buffer[V_H_NFFT] *= ac_buffer[V_H_NFFT]; // process Nyquist
                    // process positive frequencies
                    for (int k = 1; k < V_H_NFFT; ++k) {
                        ac_buffer[k] = (ac_buffer[k] * ac_buffer[k] + ac_buffer[V_NFFT-k] * ac_buffer[V_NFFT-k]); // real
                        ac_buffer[V_NFFT-k] = 0.0; // imaginary
                    }
                    mayer_realifft(V_NFFT, ac_buffer); // IFFT
                    
                    // peak finding
                    f0 = find_f0_SNAC();
                    last_f0s[last_f0s_pointer] = f0;
                    last_f0s_pointer = (last_f0s_pointer + 1) & last_f0s_mask;
                    // if f0 is stable, process bp if necessary and set delay processing to true
                    if (f0Stable()) {
                        // TODO: it seems for some signals this is too narrow!
                        if (!bp_initialized) {initialize_bp(0.90*f0, 1.10*f0);}
                        if (previous_f0_count < averaging_frames) {
                            previous_f0_sum += f0;
                            previous_f0_count += 1;
                        }
                        process_delay = true;
                    }
                }
                // now reset the buffer
                memset(ac_buffer, 0, V_NFFT * sizeof(float));
                // TODO: reset f0 buffer? if the logic all holds then we shouldn't ever need to reset it
            }
        }
    }
    
    // processing loop - main input (left channel for testing)
    auto* inputData = buffer.getWritePointer (0);
    auto* scData = buffer.getWritePointer(1);
    
    process_delay = process_delay && !bufferTooQuiet(scData, blockSize);
    // putting the if statement outside of the processing loop
    if (process_delay) {
        float w0 = T * twopi * (previous_f0_sum / previous_f0_count);
        // STEP 1: buffer the input
        for (int i = 0; i < blockSize; ++i) {
            del_buffer[write_pointer] = inputData[i];
            write_pointer = (write_pointer + 1) & del_length_mask;
        }
        
        // STEP 2: f0 bandpass filter and perform the hilbert transform
        // filters initialized with Butterworth cascade in reverse order
        f0_bandpass[1].processBlock(scData, bp_buf, blockSize);
        f0_bandpass[0].processBlockInPlace(bp_buf, blockSize);
        hilbert_left[0].processBlock(bp_buf, hilbert_left_buf, blockSize);
        hilbert_left[1].processBlockInPlace(hilbert_left_buf, blockSize);
        hilbert_left[2].processBlockInPlace(hilbert_left_buf, blockSize);
        hilbert_left[3].processBlockInPlace(hilbert_left_buf, blockSize);
        hilbert_right[0].processBlock(bp_buf, hilbert_right_buf, blockSize);
        hilbert_right[1].processBlockInPlace(hilbert_right_buf, blockSize);
        hilbert_right[2].processBlockInPlace(hilbert_right_buf, blockSize);
        hilbert_right[3].processBlockInPlace(hilbert_right_buf, blockSize);
        
        // STEP 3: calculate the phase, envelope and filter envelope
        phase_buf[0] = atan2f(last_right_out, hilbert_left_buf[0]);
        env_buf[0] = sqrtf(hilbert_left_buf[0]*hilbert_left_buf[0]+last_right_out*last_right_out);
        for (int i = 1; i < blockSize; ++i) {
            phase_buf[i] = atan2f(hilbert_right_buf[i-1], hilbert_left_buf[i]);
            env_buf[i] = sqrtf(hilbert_left_buf[i]*hilbert_left_buf[i] +
                               hilbert_right_buf[i-1]*hilbert_right_buf[i-1]);
        }
        last_right_out = hilbert_right_buf[blockSize-1];
        envelopeBP[1].processBlockInPlace(env_buf, blockSize);
        envelopeBP[0].processBlockInPlace(env_buf, blockSize);
        
        // STEP 4: perform transfer
        bool performTransfer = blocks_processed >= onset_time_blocks;
        for (int i = 0; i < blockSize; ++i) {
            // TODO: add a clip around last_env
            float env = performTransfer ?
                        (i == 0 ? last_env : AMP_CNST + amp_scaler*env_buf[i-1]) :
                        AMP_CNST;
            inputData[i] = make_up_gain * env * fractional_delay_read(read_pointer);
            //inputData[i] = make_up_gain * fractional_delay_read(read_pointer);
            
            float inst_freq = i == 0 ? phase_buf[0] - last_phase : phase_buf[i] - phase_buf[i-1];
            // phase wrapping - too aggressive?
            while (inst_freq < 0) {inst_freq += twopi;}
            while (inst_freq > twopi) {inst_freq -= twopi;}
            float dt = 1 - (inst_freq/w0);
            read_pointer = performTransfer ? (read_pointer + 1) - dt_scaler*dt : (read_pointer + 1);
            read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            
            // visualize dt
            dt_buffer.setSample(0, i, performTransfer ? dt_scaler*dt : 0.f);
            amp_buffer.setSample(0, i, performTransfer ? AMP_CNST + amp_scaler*env : AMP_CNST);
        }
        // post loop cleanup
        last_phase = phase_buf[blockSize-1];
        last_env = env_buf[blockSize-1];
        blocks_processed += 1;
        last_env = performTransfer ? AMP_CNST + amp_scaler*env_buf[blockSize-1] : AMP_CNST;
    // not processing delay because sidechain is quiet or unstable
    } else {
        // see if we need to subtly reset the read pointer
        float cur_lag = (write_pointer - read_pointer);
        while (cur_lag < 0) {cur_lag += del_length;}
        while (cur_lag > del_length) {cur_lag -= del_length;}
        float dt_diff = DEL_LAG - cur_lag;
        // try to catch up, but not too fast
        float dt = fmin(fabsf(dt_diff/blockSize), 0.002);
        dt = dt_diff > 0 ? dt : -dt;
        // TODO: right here, keep track of last envelope, start ramping back to 1
        float e_diff = AMP_CNST - last_env;
        float env_incr = fmin(fabsf(e_diff/blockSize),0.001);
        env_incr = e_diff < 0 ? -env_incr : env_incr;
        
        for (int i = 0; i < blockSize; ++i) {
            // STEP 1: buffer the input
            del_buffer[write_pointer] = inputData[i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            // STEP 2: output from the buffer while catching up dt and envelope if necessary
            inputData[i] = make_up_gain * last_env * fractional_delay_read(read_pointer);
            // read pointer can still be fractional
            read_pointer = read_pointer + 1 - dt;
            read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            
            // visualize dt
            dt_buffer.setSample(0, i, dt);
            amp_buffer.setSample(0, i, last_env);

            last_env += env_incr;
        }
    }
    del_vis.pushBuffer(dt_buffer);
    amp_vis.pushBuffer(amp_buffer);
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
        norm = norm - (sc_buffer[i-1]*sc_buffer[i-1] + sc_buffer[V_H_NFFT-i]*sc_buffer[V_H_NFFT-i]);
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

bool VibratoTransferAudioProcessor::bufferTooQuiet(auto* data, int size) {
    float l_thresh = 0.001f; // -60 dB TODO: parameterize?
    int s_thresh = int(0.1*V_H_NFFT); // 10% of signal above thresh //TODO: better cheap env follower
    int samps_above_thresh = 0;
    for (int i = 0; i < size; ++i) {
        auto samp = data[i] < 0 ? -data[i] : data[i];
        samps_above_thresh = samp > l_thresh ? samps_above_thresh + 1 : samps_above_thresh;
    }
    return samps_above_thresh < s_thresh;
}

bool VibratoTransferAudioProcessor::f0Stable() {
    // have we done at least four non-zero analyses
    float thresh = 0.05 * last_f0s[last_f0s_mask]; // this is also maybe too narrow
    if (thresh == 0.f) { return false; }
    // really bad
    bool stable = fabsf(last_f0s[0] - last_f0s[last_f0s_mask]) < thresh;
    for (int i = 1; i < last_f0s_mask; ++i) {
        stable = stable && fabsf(last_f0s[i] - last_f0s[i-1]) < thresh;
    }
    return stable;
}

void VibratoTransferAudioProcessor::initialize_bp(float bp_low, float bp_high) {
    double gain = 0.0;
    bp_initialized = bp_designer.bandPass(fs, bp_low, bp_high, 2, f0_bandpass, gain);
    if (bp_initialized) {
        f0_bandpass[1].gainCorrectNumerator(gain);
    }
}

// using fixed audio rate envelopes
void VibratoTransferAudioProcessor::initialize_env_bp() {
    double gain = 0.0;
    bool initialized = bp_designer.bandPass(fs, 1.0, 10.0, 2, envelopeBP, gain);
    if (initialized) {
        envelopeBP[1].gainCorrectNumerator(gain);
    }
}

VibVisualizer& VibratoTransferAudioProcessor::getDelayVisualizer() {
    return del_vis;
}

VibVisualizer& VibratoTransferAudioProcessor::getAmpVisualizer() {
    return amp_vis;
}

/*
// OLD CODE: peaking biquad
void VibratoTransferAudioProcessor::initialize_env_bp() {
    float f0 = 5.f;
    float bw = 4.f;
    float w0 = (twopi * f0) / fs;
    float bwr = (twopi * bw) / fs;
    float gain = 1.f / sqrt(2); // don't really need to compute every time
    float beta = (sqrt(1.f - (gain*gain))/gain) * tanf(bwr*0.5f);
    float norm_gain = 1.f/(1+beta);
    float ng_cos_w0_neg2 = -2.f*norm_gain*cosf(w0);
    envelopeBP.setParams(1.f-norm_gain, 0.f, norm_gain-1.f, ng_cos_w0_neg2, 2.f*norm_gain-1.f);
}
*/
