/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class VibratoTransferAudioProcessorEditor  : public juce::AudioProcessorEditor, private juce::Slider::Listener
{
public:
    VibratoTransferAudioProcessorEditor (VibratoTransferAudioProcessor&);
    ~VibratoTransferAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    VibratoTransferAudioProcessor& audioProcessor;
    juce::Slider ampSlider;
    juce::Slider dtSlider;
    void sliderValueChanged (juce::Slider* slider) override; // for Slider::Listener

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VibratoTransferAudioProcessorEditor)
};
