/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
VibratoTransferAudioProcessorEditor::VibratoTransferAudioProcessorEditor (VibratoTransferAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (600, 300);
    
    // code modified from https://docs.juce.com/master/tutorial_code_basic_plugin.html
    
    // these define the parameters of our slider object
    ampSlider.setSliderStyle (juce::Slider::LinearBarVertical);
    ampSlider.setTextBoxStyle (juce::Slider::TextBoxAbove, false, 90, ampSlider.getTextBoxHeight());
    ampSlider.setValue(1.0);
    amLabel.setText("AM", juce::dontSendNotification);
    amLabel.attachToComponent(&ampSlider, false);
 
    dtSlider.setSliderStyle (juce::Slider::LinearBarVertical);
    dtSlider.setTextBoxStyle (juce::Slider::TextBoxAbove, false, 90, dtSlider.getTextBoxHeight());
    dtSlider.setValue(1.0);
    fmLabel.setText("FM", juce::dontSendNotification);
    fmLabel.attachToComponent(&dtSlider, false);
    
    mugSlider.setSliderStyle (juce::Slider::LinearBarVertical);
    mugSlider.setTextBoxStyle (juce::Slider::TextBoxAbove, false, 90, mugSlider.getTextBoxHeight());
    mugSlider.setTextValueSuffix (" dB");
    mugSlider.setValue(0.0);
    mugLabel.setText("Gain", juce::dontSendNotification);
    mugLabel.attachToComponent(&mugSlider, false);
    
    fmScalerAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>(audioProcessor.getVTS(), "fmScaler", dtSlider);
    amScalerAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>(audioProcessor.getVTS(), "amScaler", ampSlider);
    makeUpGainAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>(audioProcessor.getVTS(), "makeUpGain", mugSlider);
 
    addAndMakeVisible (&ampSlider);
    addAndMakeVisible (&dtSlider);
    addAndMakeVisible (&mugSlider);
    
    addAndMakeVisible(audioProcessor.getDelayVisualizer());
    addAndMakeVisible(audioProcessor.getAmpVisualizer());
    
    fmVizLabel.setText("Pitch Shift", juce::dontSendNotification);
    fmVizLabel.attachToComponent(&(audioProcessor.getDelayVisualizer()), false);
    amVizLabel.setText("AM Envelope", juce::dontSendNotification);
    amVizLabel.attachToComponent(&(audioProcessor.getAmpVisualizer()), false);
}

VibratoTransferAudioProcessorEditor::~VibratoTransferAudioProcessorEditor()
{
}

//==============================================================================
void VibratoTransferAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setColour (juce::Colours::white);
    g.setFont (15.0f);
    //g.drawFittedText ("VibratoTransfer", getLocalBounds(), juce::Justification::centred, 1);
}

void VibratoTransferAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
    //ampSlider.setBounds (30, 30, 20, getHeight() - 60);
    //dtSlider.setBounds (60, 30, 20, getHeight() - 60);
    //mugSlider.setBounds (90, 30, 20, getHeight() - 60);
    ampSlider.setBounds (5, 30, 50, getHeight() - 60);
    dtSlider.setBounds (60, 30, 50, getHeight() - 60);
    mugSlider.setBounds (115, 30, 50, getHeight() - 60);
    audioProcessor.getDelayVisualizer().setBounds(170, 30, getWidth() - 190, getHeight()/2 - 40);
    audioProcessor.getAmpVisualizer().setBounds(170, getHeight()/2 + 10, getWidth() - 190, getHeight()/2 - 40);
}

// keeping this code around to reference later, delete after it has been put in the right spot
// audioProcessor.make_up_gain = powf(10.f, mugSlider.getValue()/20.f);
