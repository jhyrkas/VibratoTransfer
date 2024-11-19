//
//  Biquad.cpp
//  VibratoTransfer - VST3
//
//  Created by Jeremy Hyrkas on 11/18/24.
//

#include "Biquad.hpp"
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
    float outp = b0 * sample + b1 * b_zi[0] + b2 * b_zi[1] -
    a1 * a_zi[0] - a2 * a_zi[1];
    b_zi[1] = b_zi[0];
    b_zi[0] = sample;
    a_zi[1] = a_zi[0];
    a_zi[0] = outp;
    return outp;
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
