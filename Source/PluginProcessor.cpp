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
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                       .withInput("Sidechain", juce::AudioChannelSet::stereo(), true)
                       ),
#endif
b_designer(),
del_vis(),
amp_vis(),
dt_buffer(1, MAX_BUF),
amp_buffer(1, MAX_BUF),
parameters(*this, nullptr, juce::Identifier("VibratoTransfer"),
           {
            std::make_unique<juce::AudioParameterFloat>
            (juce::ParameterID { "fmScaler",  1 }, "FM Scaler", 0.0f, 2.0f, 1.f),
            std::make_unique<juce::AudioParameterFloat>
            (juce::ParameterID { "amScaler",  1 }, "AM Scaler", 0.0f, 10.0f, 1.f),
            std::make_unique<juce::AudioParameterFloat>
            (juce::ParameterID { "makeUpGain",  1 }, "Make-up Gain", -6.f, 6.f, 0.f)
            }
)
{
    // zero out delay buffer
    memset(del_buffer_l, 0, del_length * sizeof(float));
    memset(del_buffer_r, 0, del_length * sizeof(float));
    memset(sc_buffer, 0, V_H_NFFT * sizeof(float));
    memset(ac_buffer, 0, V_NFFT*sizeof(float));
    memset(dt_buf, 0, MAX_BUF*sizeof(float));
    memset(env_buf, 0, MAX_BUF*sizeof(float));
    
    fmScaler = parameters.getRawParameterValue("fmScaler");
    amScaler = parameters.getRawParameterValue("amScaler");
    makeUpGain = parameters.getRawParameterValue("makeUpGain");
    parameters.state.addListener(this);
    
    // old code, different parameter setup
    // TODO: look into this normalize range business
    // https://docs.juce.com/master/tutorial_audio_parameter.html
    //addParameter (fmScaler = new juce::AudioParameterFloat ({ "fmScaler", 1 }, "FM Scaler", 0.0f, 2.0f, 1.f));
    //addParameter (amScaler = new juce::AudioParameterFloat ({ "amScaler", 1 }, "AM Scaler", 0.0f, 10.0f, 1.f));
    //addParameter (makeUpGain = new juce::AudioParameterFloat ({ "makeUpGain", 1 }, "Make-up Gain", -6.f, 6.f, 0.f));
    
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
    return false;
}

bool VibratoTransferAudioProcessor::producesMidi() const
{
    return false;
}

bool VibratoTransferAudioProcessor::isMidiEffect() const
{
    return false;
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
    initialize_dt_hp();
    
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
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;

    // sidechain can be mono or stereo, so don't check?
    return true;
}
#endif

void VibratoTransferAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();
    
    // parameter updates
    if (updateParams) {
        dt_scaler = fmScaler->load();
        amp_scaler = amScaler->load();
        make_up_gain = powf(10.f, makeUpGain->load()/20.f);
        updateParams = false;
    }
    
    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());
    
    juce::AudioSampleBuffer mainInputOutput = getBusBuffer(buffer, true, 0);
    juce::AudioSampleBuffer sideChainInput  = getBusBuffer(buffer, true, 1);
    
    bool performTransfer = blocks_processed >= onset_time_blocks;
    int sc_channels = sideChainInput.getNumChannels();
    float sc_norm = 1.f / sc_channels;
    auto* sc_buffers = sideChainInput.getArrayOfWritePointers();
    // analysis loop - side chain (right channel for testing)

    // STEP 1: buffer sidechain signal
    for (int i = 0; i < blockSize; ++i) {
        // sum to mono, normalize
        for (int channel = 0; channel < sc_channels; ++channel) {
            float samp = sc_norm * sc_buffers[channel][i];
            sc_buffer[sc_pointer] += samp;
            ac_buffer[sc_pointer] += samp;
        }
        bp_buf[i] = sc_buffer[sc_pointer]; // store in bandpass buffer
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
                    resetProcessing();
                }
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
                } else {
                    if (bp_initialized) {
                        resetProcessing();
                    }
                }
            }
            // now reset the buffers
            memset(ac_buffer, 0, V_NFFT * sizeof(float));
            memset(sc_buffer, 0, V_H_NFFT * sizeof(float));
        }
    }
    // processing loop - main input (left channel for testing)
    auto* inputData = mainInputOutput.getArrayOfWritePointers();
    
    process_delay = process_delay && !bufferTooQuiet(bp_buf, blockSize);
#if 1 // block processing
    // putting the if statement outside of the processing loop
    if (process_delay) {
        // STEP 1: f0 bandpass filter and perform the hilbert transform
        // filters initialized with Butterworth cascade in reverse order
        f0_bandpass[1].processBlockInPlace(bp_buf, blockSize); // contains sidechain at first
        f0_bandpass[0].processBlockInPlace(bp_buf, blockSize);
        hilbert_left[0].processBlock(bp_buf, hilbert_left_buf, blockSize);
        hilbert_left[1].processBlockInPlace(hilbert_left_buf, blockSize);
        hilbert_left[2].processBlockInPlace(hilbert_left_buf, blockSize);
        hilbert_left[3].processBlockInPlace(hilbert_left_buf, blockSize);
        hilbert_right[0].processBlock(bp_buf, hilbert_right_buf, blockSize);
        hilbert_right[1].processBlockInPlace(hilbert_right_buf, blockSize);
        hilbert_right[2].processBlockInPlace(hilbert_right_buf, blockSize);
        hilbert_right[3].processBlockInPlace(hilbert_right_buf, blockSize);
        
        // STEP 2: calculate the phase, dt, envelope and filter envelope
        float w0 = T * twopi * (previous_f0_sum / previous_f0_count);
        float last_phase_tmp = last_phase;
        last_phase = atan2f(last_right_out, hilbert_left_buf[0]);
        float inst_freq = last_phase - last_phase_tmp;
        // phase wrapping - too aggressive?
        while (inst_freq < 0) {inst_freq += twopi;}
        while (inst_freq > twopi) {inst_freq -= twopi;}
        dt_buf[0] = performTransfer ? 1 - (inst_freq/w0) : 0;
        env_buf[0] = sqrtf(hilbert_left_buf[0]*hilbert_left_buf[0]+last_right_out*last_right_out);
        for (int i = 1; i < blockSize; ++i) {
            last_phase_tmp = last_phase;
            last_phase = atan2f(hilbert_right_buf[i-1], hilbert_left_buf[i]);
            inst_freq = last_phase - last_phase_tmp;
            // phase wrapping - too aggressive?
            while (inst_freq < 0) {inst_freq += twopi;}
            while (inst_freq > twopi) {inst_freq -= twopi;}
            dt_buf[i] = performTransfer ? 1 - (inst_freq/w0) : 0;
            
            env_buf[i] = sqrtf(hilbert_left_buf[i]*hilbert_left_buf[i] +
                               hilbert_right_buf[i-1]*hilbert_right_buf[i-1]);
        }
        last_right_out = hilbert_right_buf[blockSize-1];
        envelopeBP[1].processBlockInPlace(env_buf, blockSize);
        envelopeBP[0].processBlockInPlace(env_buf, blockSize);
        dt_highpass[0].processBlockInPlace(dt_buf, blockSize);
        
        // STEP 3: perform transfer
        for (int i = 0; i < blockSize; ++i) {
            // buffer input
            del_buffer_l[write_pointer] = inputData[0][i];
            del_buffer_r[write_pointer] = inputData[1][i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            
            // TODO: add a clip around last_env
            float env = performTransfer ?
                        (i == 0 ? AMP_CNST + amp_scaler*last_env : AMP_CNST + amp_scaler*env_buf[i-1]) :
                        AMP_CNST;
            // output from buffer
            float left,right;
            fractional_delay_read(read_pointer, left, right);
            inputData[0][i] = make_up_gain * env * left;
            inputData[1][i] = make_up_gain * env * right;

            read_pointer = performTransfer ? (read_pointer + 1) - dt_scaler*dt_buf[i] : (read_pointer + 1);
            //read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            // slow but double checking
            read_pointer = fmod(read_pointer, del_length);
            
            // read pointer too close to write pointer
            if (read_pointer - write_pointer <= 1.f && read_pointer - write_pointer >= 0.f) {
                if (dt_buf[i] > 0) {
                    bool tmp = true; // no-op: delay exceeded delay length
                } else {
                    bool tmp = false; // no-op: delay caught up to write head
                }
            }
            
            // visualize dt
            dt_buffer.setSample(0, i, performTransfer ? dt_scaler*dt_buf[i] : 0.f);
            amp_buffer.setSample(0, i, performTransfer ? env : AMP_CNST);
        }
        // post loop cleanup
        last_env = env_buf[blockSize-1];
        blocks_processed += 1;
        last_env = env_buf[blockSize-1];
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
        float e_diff = AMP_CNST - last_env;
        float env_incr = fmin(fabsf(e_diff/blockSize),0.001);
        env_incr = e_diff < 0 ? -env_incr : env_incr;
        
        for (int i = 0; i < blockSize; ++i) {
            // buffer the input
            del_buffer_l[write_pointer] = inputData[0][i];
            del_buffer_r[write_pointer] = inputData[1][i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            // output from the buffer while catching up dt and envelope if necessary
            float left,right;
            fractional_delay_read(read_pointer, left, right);
            inputData[0][i] = make_up_gain * last_env * left;
            inputData[1][i] = make_up_gain * last_env * right;
            // read pointer can still be fractional
            read_pointer = read_pointer + 1 - dt;
            read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            
            // visualize dt
            dt_buffer.setSample(0, i, dt);
            amp_buffer.setSample(0, i, last_env);

            last_env += env_incr;
        }
    }
#else // sample-by-sample processing
    // putting the if statement outside of the processing loop
    if (process_delay) {
        float w0 = T * twopi * (previous_f0_sum / previous_f0_count);
        for (int i = 0; i < blockSize; ++i) {
            // STEP 1: buffer the input
            del_buffer_l[write_pointer] = inputData[0][i];
            del_buffer_r[write_pointer] = inputData[1][i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            
            // STEP 2: read from the delay line
            float left,right;
            fractional_delay_read(read_pointer, left, right);
            inputData[0][i] = make_up_gain * last_env * left;
            inputData[1][i] = make_up_gain * last_env * right;
            
            // STEP 3: calculate the next delay based on the output of the bandpass filter and f0 analysis
            // filters initialized with Butterworth cascade in reverse order
            float bp_sample = f0_bandpass[0].processSample(f0_bandpass[1].processSample(bp_buf[i]));
            float hb_left = hilbert_left[3].processSample(hilbert_left[2].processSample(hilbert_left[1].processSample(hilbert_left[0].processSample(bp_sample))));
            float hb_right = hilbert_right[3].processSample(hilbert_right[2].processSample(hilbert_right[1].processSample(hilbert_right[0].processSample(bp_sample))));
            float curr_phase = atan2f(last_right_out, hb_left);
            float inst_freq = curr_phase - last_phase;
            // phase wrapping - too aggressive?
            while (inst_freq < 0) {inst_freq += twopi;}
            while (inst_freq > twopi) {inst_freq -= twopi;}
            last_phase = curr_phase;
            float dt = performTransfer ? dt_highpass[0].processSample(1 - (inst_freq/w0)) : 0.f;
            
            // filters initialized with Butterworth cascade in reverse order
            float envBPout = envelopeBP[0].processSample(
                                envelopeBP[1].processSample(
                                    sqrt(hb_left*hb_left+last_right_out*last_right_out)
                                                            )
                                                         );
            // TODO: add a clip around last_env
            last_env = performTransfer ? AMP_CNST + amp_scaler*envBPout : AMP_CNST;
            read_pointer = (read_pointer + 1) - dt_scaler*dt;
            read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            last_right_out = hb_right;
            
            // visualize dt
            dt_buffer.setSample(0, i, performTransfer ? dt_scaler*dt : 0.f);
            amp_buffer.setSample(0, i, performTransfer ? AMP_CNST + amp_scaler*envBPout : AMP_CNST);
        }
        blocks_processed += 1;
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
        
        float e_diff = AMP_CNST - last_env;
        float env_incr = fmin(fabsf(e_diff/blockSize),0.001);
        env_incr = e_diff < 0 ? -env_incr : env_incr;
        
        for (int i = 0; i < blockSize; ++i) {
            // STEP 1: buffer the input
            del_buffer_l[write_pointer] = inputData[0][i];
            del_buffer_r[write_pointer] = inputData[1][i];
            write_pointer = (write_pointer + 1) & del_length_mask;
            // STEP 2: output from the buffer while catching up dt and envelope if necessary
            float left,right;
            fractional_delay_read(read_pointer, left, right);
            inputData[0][i] = make_up_gain * last_env * left;
            inputData[1][i] = make_up_gain * last_env * right;
            // read pointer can still be fractional
            read_pointer = read_pointer + 1 - dt;
            read_pointer = ((int)read_pointer & del_length_mask) + (read_pointer - (int)read_pointer);
            last_env += env_incr;
            
            // visualize dt
            dt_buffer.setSample(0, i, dt);
            amp_buffer.setSample(0, i, last_env);
        }
    }
#endif
    del_vis.pushBuffer(dt_buffer);
    amp_vis.pushBuffer(amp_buffer);
    
    for (int i = 1; i < blockSize; ++i) {
        if (fabsf(dt_buffer.getSample(0, i) - dt_buffer.getSample(0, i-1)) > 0.01) {
            bool tmp = true; // no-op
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

// output values get put into left and right
void VibratoTransferAudioProcessor::fractional_delay_read(float index, float& left, float& right) {
    // TODO: should we ever check that index is in range, or should that be done
    // TODO: outside of the function?
    int low = int(index);
    int high = low == del_length - 1 ? 0 : low+1;
    float high_frac = index - low;
    left = (1-high_frac) * del_buffer_l[low] + high_frac * del_buffer_l[high];
    right = (1-high_frac) * del_buffer_r[low] + high_frac * del_buffer_r[high];
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
    
    // parabolic interpolation, based on FFT interp
    // https://ccrma.stanford.edu/~jos/sasp/Quadratic_Interpolation_Spectral_Peaks.html
    float alpha = ac_buffer[peak_index-1];
    float beta = ac_buffer[peak_index];
    float gamma = ac_buffer[peak_index+1];
    float peak_index_interp = peak_index + 0.5 * (alpha - gamma)/(alpha - 2*beta + gamma);
    
    //return fs / peak_index;
    return fs / peak_index_interp;
}

bool VibratoTransferAudioProcessor::bufferTooQuiet(auto* data, int size) {
    float l_thresh = 0.001f; // -60 dB TODO: parameterize?
    int s_thresh = int(0.1*V_H_NFFT); // 10% of signal above thresh //TODO: better cheap env follower
    int samps_above_thresh = 0;
    for (int i = 0; i < size; ++i) {
        auto samp = data[i] < 0 ? -data[i] : data[i];
        samps_above_thresh = samp > l_thresh ? samps_above_thresh + 1 : samps_above_thresh;
    }
    
    // if we are currently processing and fall below the threshold
    if (process_delay && samps_above_thresh < s_thresh) {
        bool tmp = samps_above_thresh < s_thresh; // no-op for breakpoint
    }
    
    return samps_above_thresh < s_thresh;
}

bool VibratoTransferAudioProcessor::f0Stable() {
    // have we done at least four non-zero analyses
    float thresh = 0.1 * last_f0s[last_f0s_mask]; // this is also maybe too narrow
    if (thresh == 0.f) { return false; }
    // really bad
    int stableCount = fabsf(last_f0s[0] - last_f0s[last_f0s_mask]) < thresh;
    for (int i = 1; i < last_f0s_mask; ++i) {
        stableCount += fabsf(last_f0s[i] - last_f0s[i-1]) < thresh;
    }
    
    // if we are currently processing and f0 becomes unstable
    if (process_delay && stableCount < last_f0s_mask - 1) {
        bool tmp = false; // no-op for breakpoint
    }
    
    // six of eight within threshold
    return stableCount >= last_f0s_mask - 1;
}

void VibratoTransferAudioProcessor::initialize_bp(float bp_low, float bp_high) {
    double gain = 0.0;
    bp_initialized = b_designer.bandPass(fs, bp_low, bp_high, 2, f0_bandpass, gain);
    if (bp_initialized) {
        f0_bandpass[1].gainCorrectNumerator(gain);
    }
}

void VibratoTransferAudioProcessor::initialize_env_bp() {
    double gain = 0.0;
    bool initialized = b_designer.bandPass(fs, 1.0, 10.0, 2, envelopeBP, gain);
    if (initialized) {
        envelopeBP[1].gainCorrectNumerator(gain);
    }
}

void VibratoTransferAudioProcessor::initialize_dt_hp() {
    double gain = 0.0;
    bool initialized = b_designer.hiPass(fs, 1.0, 0.f, 2, dt_highpass, gain);
    if (initialized) {
        dt_highpass[0].gainCorrectNumerator(gain);
    }
}

VibVisualizer& VibratoTransferAudioProcessor::getDelayVisualizer() {
    return del_vis;
}

VibVisualizer& VibratoTransferAudioProcessor::getAmpVisualizer() {
    return amp_vis;
}

void VibratoTransferAudioProcessor::valueTreePropertyChanged(juce::ValueTree &treeWhosePropertyHasChanged, const juce::Identifier &property) {
    updateParams = true;
}

void VibratoTransferAudioProcessor::resetProcessing() {
    bp_initialized = false;
    process_delay = false;
    last_f0s_pointer = 0;
    memset(last_f0s, 0, NUM_F0 * sizeof(float));
    previous_f0_sum = 0;
    previous_f0_count = 0;
    for (Biquad quad : f0_bandpass) {
        quad.clear();
    }
    for (Biquad quad : envelopeBP) {
        quad.clear();
    }
    for (Biquad quad : dt_highpass) {
        quad.clear();
    }
    hilbert_left[0].clear(); hilbert_left[1].clear(); hilbert_left[2].clear(); hilbert_left[3].clear();
    hilbert_right[0].clear(); hilbert_right[1].clear(); hilbert_right[2].clear(); hilbert_right[3].clear();
    blocks_processed = 0;
}
