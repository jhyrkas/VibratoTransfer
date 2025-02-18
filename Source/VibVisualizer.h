/*
  ==============================================================================

    VibVisualizer.h
    Created: 28 Jan 2025 10:46:58am
    Author:  Jeremy Hyrkas

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>

//#include <juce_AudioVisualiserComponent.h>

class VibVisualizer : public juce::AudioVisualiserComponent {
public:
    VibVisualizer();
    virtual ~VibVisualizer() override;
private:
};
