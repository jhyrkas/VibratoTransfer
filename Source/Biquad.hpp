//
//  Biquad.hpp
//  VibratoTransfer - VST3
//
//  Created by Jeremy Hyrkas on 11/18/24.
//

#ifndef Biquad_hpp
#define Biquad_hpp

class Biquad
{
public:
    Biquad(float b_coefs[3], float a_coefs[2]);
    virtual ~Biquad();
    // TODO: discuss - for simplicity it is easier to define this sample-by-sample,
    // TODO: but is it more efficient to have each biquad process a block at a time?
    // TODO: basically, this is multiple for-loops versus many function calls
    virtual float processSample(float sample);
private:
    float b[3];
    float a[2];
    // TODO: reimplement in canonical form for fewer floats kept around
    float b_zi[2];
    float a_zi[2];
};

#endif /* Biquad_hpp */
