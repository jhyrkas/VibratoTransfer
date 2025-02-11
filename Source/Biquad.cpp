//
//  Biquad.cpp
//  VibratoTransfer - VST3
//
//  Created by author on 11/18/24.
//

#include "Biquad.hpp"

Biquad::Biquad() :
b0(1), b1(0), b2(0), a1(0), a2(0)
{
    b_zi[0] = 0; b_zi[1] = 0;
    a_zi[0] = 0; a_zi[1] = 0;
}

Biquad::Biquad(float _b0, float _b1, float _b2, float _a1, float _a2) :
b0(_b0), b1(_b1), b2(_b2), a1(_a1), a2(_a2)
{
    b_zi[0] = 0; b_zi[1] = 0;
    a_zi[0] = 0; a_zi[1] = 0;
}

Biquad::~Biquad() {
    // no-op, as far as I know everything is on the stack
}

float Biquad::processSample(float sample) {
    float outp = b0 * sample + b1 * b_zi[0] + b2 * b_zi[1]
        - a1 * a_zi[0] - a2 * a_zi[1];
    b_zi[1] = b_zi[0];
    b_zi[0] = sample;
    a_zi[1] = a_zi[0];
    a_zi[0] = outp;
    return outp;
}

void Biquad::processBlock(float* in_samples, float* out_samples, int num_samples) {
    if (num_samples < 2) return; // this is being ridiculously safe
    // first sample
    out_samples[0] = b0 * in_samples[0] + b1 * b_zi[0] + b2 * b_zi[1] - a1 * a_zi[0] - a2 * a_zi[1];
    // second sample
    out_samples[1] = b0 * in_samples[1] + b1 * in_samples[0] + b2 * b_zi[0] - a1 * out_samples[0] - a2 * a_zi[0];
    // rest of samples
    for (int i = 2; i < num_samples; ++i) {
        out_samples[i] = b0 * in_samples[i] + b1 * in_samples[i-1] + b2 * in_samples[i-2]
                        - a1 * out_samples[i-1] - a2 * out_samples[i-2];
    }
    // set b_zi and a_zi
    b_zi[0] = in_samples[num_samples-1];
    b_zi[1] = in_samples[num_samples-2];
    a_zi[0] = out_samples[num_samples-1];
    a_zi[1] = out_samples[num_samples-2];
}

// making this a separate method in case the compiler optimizer does anything wild up above
// under the assumption that in_samples and out_samples must be different
void Biquad::processBlockInPlace (float* in_samples, int num_samples) {
    if (num_samples < 2) return; // this is being ridiculously safe
    // first sample
    float sample = in_samples[0];
    in_samples[0] = b0 * sample + b1 * b_zi[0] + b2 * b_zi[1] - a1 * a_zi[0] - a2 * a_zi[1];
    b_zi[1] = b_zi[0];
    b_zi[0] = sample;
    // second sample
    sample = in_samples[1];
    in_samples[1] = b0 * sample + b1 * b_zi[0] + b2 * b_zi[1] - a1 * in_samples[0] - a2 * a_zi[0];
    b_zi[1] = b_zi[0];
    b_zi[0] = sample;
    for (int i = 2; i < num_samples; ++i) {
        float sample = in_samples[i];
        in_samples[i] = b0 * sample + b1 * b_zi[0] + b2 * b_zi[1]
                        - a1 * in_samples[i-1] - a2 * in_samples[i-2];
        b_zi[1] = b_zi[0];
        b_zi[0] = sample;
    }
    // set a_zi
    a_zi[0] = in_samples[num_samples-1];
    a_zi[1] = in_samples[num_samples-2];
}

// TODO: maybe in the future we allow the zi to stick around?
void Biquad::setParams(float _b0, float _b1, float _b2, float _a1, float _a2) {
    b0 = _b0;
    b1 = _b1;
    b2 = _b2;
    a1 = _a1;
    a2 = _a2;
    b_zi[0] = 0; b_zi[1] = 0;
    a_zi[0] = 0; a_zi[1] = 0;
}

void Biquad::clear() {
    b_zi[0] = 0; b_zi[1] = 0;
    a_zi[0] = 0; a_zi[1] = 0;
}

// converting from df2t form to direct form 1
void Biquad::DF2TBiquad(double B0, double B1, double B2,
                        double A0, double A1, double A2) {
    b0 = B0  / A0;
    b1 = B1  / A0;
    b2 = B2  / A0;
    a1 = (A1) / A0;
    a2 = (A2) / A0;
    b_zi[0] = 0; b_zi[1] = 0;
    a_zi[0] = 0; a_zi[1] = 0;
}

void Biquad::gainCorrectNumerator(double gain) {
    b0 *= gain;
    b1 *= gain;
    b2 *= gain;
}
