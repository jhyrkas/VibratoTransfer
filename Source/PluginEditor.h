/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include "math.h" // for log
#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class VibratoTransferAudioProcessorEditor  : public juce::AudioProcessorEditor
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
    juce::Slider mugSlider; // make up gain
    
    juce::Label amLabel;
    juce::Label fmLabel;
    juce::Label mugLabel;
    juce::Label amVizLabel;
    juce::Label fmVizLabel;
    
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fmScalerAttachment;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> amScalerAttachment;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> makeUpGainAttachment;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VibratoTransferAudioProcessorEditor)
};
