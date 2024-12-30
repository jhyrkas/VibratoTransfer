//
//  Biquad.hpp
//  VibratoTransfer - VST3
//
//  Created by author on 11/18/24.
//

#ifndef Biquad_hpp
#define Biquad_hpp

class Biquad
{
public:
    Biquad(float _b0, float _b1, float _b2, float _a1, float _a2);
    virtual ~Biquad();
    // TODO: discuss - for simplicity it is easier to define this sample-by-sample,
    // TODO: but is it more efficient to have each biquad process a block at a time?
    // TODO: basically, this is multiple for-loops versus many function calls
    virtual float processSample(float sample);
    virtual void setParams(float _b0, float _b1, float _b2, float _a1, float _a2);
    virtual void clear();
protected:
    float b0, b1, b2, a1, a2;
    // TODO: reimplement in canonical form for fewer floats kept around
    float b_zi[2];
    float a_zi[2];
};

#endif /* Biquad_hpp */
