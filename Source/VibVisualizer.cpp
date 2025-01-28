/*
  ==============================================================================

    DelayVisualizer.cpp
    Created: 28 Jan 2025 10:46:58am
    Author:  Jeremy Hyrkas

  ==============================================================================
*/

#include "VibVisualizer.h"

VibVisualizer::VibVisualizer() : juce::AudioVisualiserComponent(1)
{
    // TODO: pick right sizes for these, can we set samplesperblock based on preparetoplay in processor?
    setBufferSize(512);
    setSamplesPerBlock(64);
    setColours(juce::Colours::black, juce::Colours::indianred);
}

VibVisualizer::~VibVisualizer()
{
    
}
