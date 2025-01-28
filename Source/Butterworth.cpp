/*
 
 This file is part of Butterworth Filter Design, a pair C++ classes and an
 accompanying suite of unit tests for designing high order Butterworth IIR &
 EQ filters using the bilinear transform.
 The generated filter coefficients are split out into cascaded biquad sections,
 for easy use in your garden variety biquad or second-order section (SOS).
 
 Reference: http://en.wikipedia.org/wiki/Butterworth_filter
 http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 
 
 Copyright (C) 2013,  iroro orife
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 */


#import <iostream>
#include <math.h>
#include "Butterworth.h"

#define LOG_OUTPUT 0    // Enable output logging

using namespace std;

#pragma mark - Public

//******************************************************************************
// Lowpass analogue prototype. Places Butterworth poles evenly around
// the S-plane unit circle.
//
// Reference: MATLAB buttap(filterOrder)
//******************************************************************************

vector <complex_double>
Butterworth::prototypeAnalogLowPass(int filterOrder){
    
    vector <complex_double> poles;
    
    for(uint32_t k = 0; k < (filterOrder + 1) / 2; k++){
        double theta = (double)(2 * k + 1) * M_PI / (2 * filterOrder);
        double real = -sin(theta);
        double imag = cos(theta);
        poles.push_back(complex_double(real,  imag));
        poles.push_back(complex_double(real, -imag)); // conjugate
    }
    
    return poles;
}


//******************************************************************************
// Generate Butterworth coefficients, the main public method
//
// Reference: MATLAB butter(n, Wn, ...)
//            http://en.wikipedia.org/wiki/Butterworth_filter
//******************************************************************************

bool Butterworth::coefficients(FILTER_TYPE filter, const double fs, const double freq1_cutoff,
                               const double freq2_cutoff, const int filterOrder,
                               vector <Biquad> & coeffs,
                               double & overallGain){
    
    //******************************************************************************
    // Init internal state based on filter design requirements
    
    zeros.resize(2 * filterOrder);
    poles.resize(2 * filterOrder);
    
    f1 = freq1_cutoff;
    f2 = freq2_cutoff;
    
    Wc = 0;  // Omega cutoff = passband edge freq
    bw = 0;
    
    
    //******************************************************************************
    // Prewarp
    
    f1 = 2 * tan(M_PI * f1 / fs);
    f2 = 2 * tan(M_PI * f2 / fs);
    
    //******************************************************************************
    // Design basic S-plane poles-only analogue LP prototype
    
    // Get zeros & poles of prototype analogue low pass.
    vector <complex_double> tempPoles = prototypeAnalogLowPass(filterOrder);
    
    // Copy tmppole into poles
    int index = 0;
    for(vector <complex_double>::iterator itr = tempPoles.begin(); itr != tempPoles.end(); itr++){
        poles[index] = *itr;
        index++;
    }
    
    numPoles = (int)tempPoles.size();
    numZeros = 0;       // butterworth LP prototype has no zeros.
    gain = 1.0;         // always 1 for the butterworth prototype lowpass.
    
    
    //******************************************************************************
    // Convert prototype to target filter type (LP/HP/BP/BS) - S-plane
    
    // Re-orient BP/BS corner frequencies if necessary
    if(f1 > f2){
        double temp = f2;
        f2 = f1;
        f1 = temp;
    }
    
    // Cutoff Wc = f2
    switch(filter){
            
        case kLoPass:
            convert2lopass();
            break;
            
        case kHiPass:
            convert2hipass();
            break;
            
        case kBandPass:
            convert2bandpass();
            break;
            
        case kBandStop:
            convert2bandstop();
            break;
            
        default: {
            return false;
        }
    }
    
    
    //******************************************************************************
    // SANITY CHECK: Ensure poles are in the left half of the S-plane
    for(uint32_t i = 0; i < numPoles; i++){
        if(poles[i].real() > 0){
            return false;
        }
    }
    
    
    //******************************************************************************
    // Map zeros & poles from S-plane to Z-plane
    
    nba = 0;
    ba = new double[2 * std::max(numPoles, numZeros) + 5];
    preBLTgain = gain;
    
    if(!s2Z()){
        return false;
    }
    
    
    //******************************************************************************
    // Split up Z-plane poles and zeros into SOS
    
    if(!zp2SOS()){
        return false;
    }
    
    // correct the overall gain
    if(filter == kLoPass || filter == kBandPass){ // pre-blt is okay for S-plane
        ba[0] = preBLTgain * (preBLTgain / gain); // 2nd term is how much BLT boosts,
    }
    else if(filter == kHiPass || kBandStop){ // HF gain != DC gain
        ba[0] = 1 / ba[0];
    }
    
    
    //******************************************************************************
    // Init biquad chain with coefficients from SOS
    
    overallGain = ba[0];
    int numFilters = filterOrder / 2;
    if(filter == kBandPass || filter == kBandStop){
        numFilters = filterOrder; // we have double the # of biquad sections
        
        // IOHAVOC filterOrder is never used again? figure this out FIXME
        // filterOrder *= 2;
    }
    
    coeffs.resize(numFilters);
    for(uint32_t i = 0; i < numFilters; i++){
        (coeffs)[i].DF2TBiquad(1.0,                     // b0
                               ba[4 * i + 1],           // b1
                               ba[4 * i + 2],           // b2
                               1.0,                     // a0
                               ba[4 * i + 3],           // a1
                               ba[4 * i + 4]);          // a2
    }
    
    return true;
}

#pragma mark - Filter design utility methods

//******************************************************************************
//
// Z = (2 + S) / (2 - S) is the S-plane to Z-plane bilinear transform
//
// Reference: http://en.wikipedia.org/wiki/Bilinear_transform
//
//******************************************************************************

double Butterworth::blt(complex_double & sz){
    
    complex_double two(2.0, 0);
    complex_double s = sz;      // sz is an input(S-plane) & output(Z-plane) arg
    sz = (two + s) / (two - s);
    
    // return the gain
    return abs((two - s));
}


//******************************************************************************
//
// Convert poles & zeros from S-plane to Z-plane via Bilinear Tranform (BLT)
//
//******************************************************************************

bool Butterworth::s2Z(){
    
    // blt zeros
    for(uint32_t i = 0; i < numZeros; i++){
        gain /= blt(zeros[i]);
    }
    
    // blt poles
    for(uint32_t i = 0; i < numPoles; i++){
        gain *= blt(poles[i]);
    }
    
    return true;
}


//******************************************************************************
//
// Convert filter poles and zeros to second-order sections
//
// Reference: http://www.mathworks.com/help/signal/ref/zp2sos.html
//
//******************************************************************************

bool Butterworth::zp2SOS(){
    
    int filterOrder = std::max(numZeros, numPoles);
    complex_double * zerosTempVec = new complex_double[filterOrder];
    complex_double * polesTempVec = new complex_double[filterOrder];
    
    // Copy
    for(uint32_t i = 0; i < numZeros; i++){
        zerosTempVec[i] = zeros[i];
    }
    
    // Add zeros at -1, so if S-plane degenerate case where
    // numZeros = 0 will map to -1 in Z-plane.
    for(uint32_t i = numZeros; i < filterOrder; i++){
        zerosTempVec[i] = complex_double(-1, 0);
    }
    
    // Copy
    for(uint32_t i = 0; i < numPoles; i++){
        polesTempVec[i] = poles[i];
    }
    
    ba[0] = gain; // store gain
    
    int numSOS = 0;
    for(uint32_t i = 0; i + 1 < filterOrder; i += 2, numSOS++){
        ba[4 * numSOS + 1] = -(zerosTempVec[i] + zerosTempVec[i + 1]).real();
        ba[4 * numSOS + 2] =  (zerosTempVec[i] * zerosTempVec[i + 1]).real();
        ba[4 * numSOS + 3] = -(polesTempVec[i] + polesTempVec[i + 1]).real();
        ba[4 * numSOS + 4] =  (polesTempVec[i] * polesTempVec[i + 1]).real();
    }
    
    // Odd filter order thus one pair of poles/zeros remains
    if(filterOrder % 2 == 1){
        ba[4 * numSOS + 1] = -zerosTempVec[filterOrder - 1].real();
        ba[4 * numSOS + 2] = ba[4 * numSOS + 4] = 0;
        ba[4 * numSOS + 3] = -polesTempVec[filterOrder - 1].real();
        numSOS++;
    }
    
    // Set output param
    nba = 1 + 4 * numSOS;
    
    delete[] zerosTempVec;
    delete[] polesTempVec;
    return true;
}




#pragma mark - Analog lowpss prototype conversion methods

//******************************************************************************
// Convert analog lowpass prototype poles to lowpass
//******************************************************************************

void Butterworth::convert2lopass(){
    Wc = f2;                                   // critical frequency
    
    gain *= pow(Wc, numPoles);
    
    numZeros = 0;                              // poles only
    for(uint32_t i = 0; i < numPoles; i++){    // scale poles by the cutoff Wc
        poles[i] = Wc * poles[i];
    }
}


//******************************************************************************
// Convert lowpass poles & zeros to highpass
// with Wc = f2, use:  hp_S = Wc / lp_S;
//******************************************************************************

void Butterworth::convert2hipass(){
    Wc = f2;   // Critical frequency
    
    // Calculate gain
    complex_double prodz(1.0, 0.0);
    complex_double prodp(1.0, 0.0);
    
    for(uint32_t i = 0; i < numZeros; i++){
        prodz *= -zeros[i];
    }
    
    for(uint32_t i = 0; i < numPoles; i++){
        prodp *= -poles[i];
    }
    
    gain *= prodz.real() / prodp.real();
    
    // Convert LP poles to HP
    for(uint32_t i = 0; i < numPoles; i++){
        if(abs(poles[i])){
            poles[i] = complex_double(Wc) / poles[i]; //  hp_S = Wc / lp_S;
            
        }
    }
    
    // Init with zeros, no non-zero values to convert
    numZeros = numPoles;
    for(uint32_t i = 0; i < numZeros; i++){
        zeros[i] = complex_double(0.0);
    }
}


//******************************************************************************
// Convert lowpass poles to bandpass
// use:  bp_S = 0.5 * lp_S * BW +
//                   0.5 * sqrt ( BW^2 * lp_S^2 - 4*Wc^2 )
// where   BW = W2 - W1
//            Wc^2 = W2 * W1
//******************************************************************************

void Butterworth::convert2bandpass(){
    bw = f2 - f1;
    Wc = sqrt(f1 * f2);
    
    // Calculate bandpass gain
    gain *= pow(bw, numPoles - numZeros);
    
    // Convert LP poles to BP: these two sets of for-loops result in an ordered
    // list of poles and their complex conjugates
    vector <complex_double> tempPoles;
    for(uint32_t i = 0; i < numPoles; i++){    // First set of poles + conjugates
        if(abs(poles[i])){
            complex_double firstterm = 0.5 * poles[i] * bw;
            complex_double secondterm = 0.5 * sqrt((bw * bw) * (poles[i] * poles[i]) - 4 * Wc * Wc);
            tempPoles.push_back(firstterm + secondterm);
        }
    }
    
    for(uint32_t i = 0; i < numPoles; i++){    // Second set of poles + conjugates
        if(abs(poles[i])){
            complex_double firstterm = 0.5 * poles[i] * bw;
            complex_double secondterm = 0.5 * sqrt((bw * bw) * (poles[i] * poles[i]) - 4 * Wc * Wc);
            tempPoles.push_back(firstterm - secondterm);      // complex conjugate
        }
    }
    
    // Init zeros, no non-zero values to convert
    numZeros = numPoles;
    for(uint32_t i = 0; i < numZeros; i++){
        zeros[i] = complex_double(0.0);
    }
    
    // Copy converted poles to output array
    int index = 0;
    for(vector <complex_double>::iterator itr = tempPoles.begin(); itr != tempPoles.end(); itr++){
        poles[index] = *itr;
        index++;
    }
    numPoles = (int)tempPoles.size();
}


//******************************************************************************
// Convert lowpass poles to bandstop
// use:  bs_S = 0.5 * BW / lp_S +
//                   0.5 * sqrt ( BW^2 / lp_S^2 - 4*Wc^2 )
// where   BW = W2 - W1
//         Wc^2 = W2 * W1
//******************************************************************************

void Butterworth::convert2bandstop(){
    bw = f2 - f1;
    Wc = sqrt(f1 * f2);
    
    // Compute gain
    complex_double prodz(1.0, 0.0);
    complex_double prodp(1.0, 0.0);
    for(uint32_t i = 0; i < numZeros; i++){
        prodz *= -zeros[i];
    }
    
    for(uint32_t i = 0; i < numPoles; i++){
        prodp *= -poles[i];
    }
    
    gain *= prodz.real() / prodp.real();
    
    // Convert LP zeros to band stop
    numZeros = numPoles;
    vector <complex_double> ztmp;
    for(uint32_t i = 0; i < numZeros; i++){
        ztmp.push_back(complex_double(0.0,  Wc));
        ztmp.push_back(complex_double(0.0, -Wc));         // complex conjugate
    }
    
    vector <complex_double> tempPoles;
    for(uint32_t i = 0; i < numPoles; i++){    // First set of poles + conjugates
        if(abs(poles[i])){
            complex_double term1 = 0.5 * bw / poles[i];
            complex_double term2 = 0.5 * sqrt((bw * bw) / (poles[i] * poles[i]) - (4 * Wc * Wc));
            tempPoles.push_back(term1 + term2);
        }
    }
    
    for(uint32_t i = 0; i < numPoles; i++){    // Second set of poles + conjugates
        if(abs(poles[i])){
            complex_double term1 = 0.5 * bw / poles[i];
            complex_double term2 = 0.5 * sqrt((bw * bw) / (poles[i] * poles[i]) - (4 * Wc * Wc));
            tempPoles.push_back(term1 - term2);       // complex conjugate
        }
    }
    
    // Copy converted zeros to output array
    int index = 0;
    for(vector <complex_double>::iterator itr = ztmp.begin(); itr != ztmp.end(); itr++){
        zeros[index] = *itr;
        index++;
    }
    numZeros = (int)ztmp.size();
    
    // Copy converted poles to output array
    index = 0;
    for(vector <complex_double>::iterator itr = tempPoles.begin(); itr != tempPoles.end(); itr++){
        poles[index] = *itr;
        index++;
    }
    numPoles = (int)tempPoles.size();
}
