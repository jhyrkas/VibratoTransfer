//
//  Biquad.cpp
//  VibratoTransfer - VST3
//
//  Created by Jeremy Hyrkas on 11/18/24.
//

#include "Biquad.hpp"
Biquad::Biquad(float b_coefs[3], float a_coefs[2]) {
    b[0] = b_coefs[0];
    b[1] = b_coefs[1];
    b[2] = b_coefs[2];
    a[0] = a_coefs[0];
    a[1] = a_coefs[1];
    b_zi[0] = 0; b_zi[1] = 0;
    a_zi[0] = 0; a_zi[1] = 0;
}

Biquad::~Biquad() {
    // no-op, as far as I know everything is on the stack
}

float Biquad::processSample(float sample) {
    float outp = b[0] * sample + b[1] * b_zi[0] + b[2] * b_zi[1] -
    a[0] * a_zi[0] - a[1] * a_zi[1];
    b_zi[1] = b_zi[0];
    b_zi[0] = sample;
    a_zi[1] = a_zi[0];
    a_zi[0] = outp;
    return outp;
}
