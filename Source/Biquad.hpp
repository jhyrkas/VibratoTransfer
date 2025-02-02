//
//  Biquad.h
//  VibratoTransfer - VST3
//
//  Created by author on 11/18/24.
//

#ifndef Biquad_hpp
#define Biquad_hpp

class Biquad
{
public:
    Biquad(); // empty constructor for templates
    Biquad(float _b0, float _b1, float _b2, float _a1, float _a2);
    virtual ~Biquad();
    // TODO: discuss - for simplicity it is easier to define this sample-by-sample,
    // TODO: but is it more efficient to have each biquad process a block at a time?
    // TODO: basically, this is multiple for-loops versus many function calls
    virtual float processSample(float sample);
    virtual void processBlock (float* in_samples, float* out_samples, int num_samples);
    virtual void processBlockInPlace (float* in_samples, int num_samples);
    virtual void setParams(float _b0, float _b1, float _b2, float _a1, float _a2);
    virtual void clear();
    
    // additions from Ruoho Ruotsi's Butterworth Filter Design codebase
    // https://github.com/ruohoruotsi/Butterworth-Filter-Design
    // Ruotsi's code implements a Biquad object that is largely redundant with functionality
    // here, but has some methods that we need
    
    // Coefficients for a DF2T biquad section.
    void DF2TBiquad(double B0, double B1, double B2,
                        double A0, double A1, double A2);
    
    // gain correct
    void gainCorrectNumerator(double gain);
protected:
    float b0, b1, b2, a1, a2;
    // TODO: reimplement in canonical form for fewer floats kept around
    float b_zi[2];
    float a_zi[2];
};

#endif /* Biquad_hpp */
