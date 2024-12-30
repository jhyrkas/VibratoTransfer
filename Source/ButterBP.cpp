/*
  ==============================================================================

    ButterBP.cpp
    Created: 9 Dec 2024 9:12:40am
    Author:  author

  ==============================================================================
*/

#include "ButterBP.h"

ButterBP::ButterBP() {
    b_zi[0] = 0; b_zi[1] = 0; b_zi[2] = 0; b_zi[3] = 0;
    a_zi[0] = 0; a_zi[1] = 0; a_zi[2] = 0; a_zi[3] = 0;
}

ButterBP::~ButterBP() {
    // no-op
}

void ButterBP::setParams(float bpL, float bpH, float fs) {
    b_zi[0] = 0; b_zi[1] = 0; b_zi[2] = 0; b_zi[3] = 0;
    a_zi[0] = 0; a_zi[1] = 0; a_zi[2] = 0; a_zi[3] = 0;
    
    float wLow = (2*bpL)/fs;
    float wHigh = (2*bpH)/fs;
    std::vector<double> a = ComputeDenCoeffs(wLow, wHigh);
    std::vector<double> b = ComputeNumCoeffs(wLow, wHigh, a);
    b0 = b[0]; b1 = b[1]; b2 = b[2]; b3 = b[3]; b4 = b[4];
    a1 = a[1]; a2 = a[2]; a3 = a[3]; a4 = a[4]; // a0 = 1
}

void ButterBP::clear() {
    b_zi[0] = 0; b_zi[1] = 0; b_zi[2] = 0; b_zi[3] = 0;
    a_zi[0] = 0; a_zi[1] = 0; a_zi[2] = 0; a_zi[3] = 0;
}

// this is a very inefficient implementation, and also risks numerical errors due to no SOS
float ButterBP::processSample(float sample) {
    float outp = b0 * sample + b1 * b_zi[0] + b2 * b_zi[1] + b3 * b_zi[2] + b4 * b_zi[3] -
                 a1 * a_zi[0] - a2 * a_zi[1] - a3 * a_zi[2] - a4 * a_zi[3];
    
    // horrible swapping
    b_zi[3] = b_zi[2]; b_zi[2] = b_zi[1]; b_zi[1] = b_zi[0]; b_zi[0] = sample;
    a_zi[3] = a_zi[2]; a_zi[2] = a_zi[1]; a_zi[1] = a_zi[0]; a_zi[0] = outp;
    
    return outp;
}

// -----------------
// functions below adapted from Edson Niu's C++ implementation of the butter filter in MATLAB
// code from: https://github.com/nxsEdson/Butterworth-Filter/blob/master/butterworth.cpp
// TODO: if this code sticks around, we can switch from vector to arrays of size 5
std::vector<double> ButterBP::ComputeDenCoeffs(double Lcutoff, double Ucutoff)
{
    int k;            // loop variables
    double theta;     // PI * (Ucutoff - Lcutoff) / 2.0
    double cp;        // cosine of phi
    double st;        // sine of theta
    double ct;        // cosine of theta
    double s2t;       // sine of 2*theta
    double c2t;       // cosine 0f 2*theta
    std::vector<double> RCoeffs(2 * FILTERORDER);     // z^-2 coefficients
    std::vector<double> TCoeffs(2 * FILTERORDER);     // z^-1 coefficients
    std::vector<double> DenomCoeffs;     // dk coefficients
    double PoleAngle;      // pole angle
    double SinPoleAngle;     // sine of pole angle
    double CosPoleAngle;     // cosine of pole angle
    double a;         // workspace variables

    cp = cos(M_PI * (Ucutoff + Lcutoff) / 2.0);
    theta = M_PI * (Ucutoff - Lcutoff) / 2.0;
    st = sin(theta);
    ct = cos(theta);
    s2t = 2.0*st*ct;        // sine of 2*theta
    c2t = 2.0*ct*ct - 1.0;  // cosine of 2*theta

    for (k = 0; k < FILTERORDER; ++k)
    {
        PoleAngle = M_PI * (double)(2 * k + 1) / (double)(2 * FILTERORDER);
        SinPoleAngle = sin(PoleAngle);
        CosPoleAngle = cos(PoleAngle);
        a = 1.0 + s2t*SinPoleAngle;
        RCoeffs[2 * k] = c2t / a;
        RCoeffs[2 * k + 1] = s2t*CosPoleAngle / a;
        TCoeffs[2 * k] = -2.0*cp*(ct + st*SinPoleAngle) / a;
        TCoeffs[2 * k + 1] = -2.0*cp*st*CosPoleAngle / a;
    }

    DenomCoeffs = TrinomialMultiply(TCoeffs, RCoeffs);

    DenomCoeffs[1] = DenomCoeffs[0];
    DenomCoeffs[0] = 1.0;
    for (k = 3; k <= 2 * FILTERORDER; ++k)
        DenomCoeffs[k] = DenomCoeffs[2 * k - 2];

    // this originally seemed to always make a one coefficient larger than b
    // changing from i > FILTERORDER*2+1 to i > FILTERORDER+2, but I don't feel great about it
    for (int i = DenomCoeffs.size() - 1; i > FILTERORDER * 2; i--)
        DenomCoeffs.pop_back();

    return DenomCoeffs;
}

std::vector<double> ButterBP::TrinomialMultiply(std::vector<double> b, std::vector<double> c)
{
    int i, j;
    std::vector<double> RetVal(4 * FILTERORDER);

    RetVal[2] = c[0];
    RetVal[3] = c[1];
    RetVal[0] = b[0];
    RetVal[1] = b[1];

    for (i = 1; i < FILTERORDER; ++i)
    {
        RetVal[2 * (2 * i + 1)] += c[2 * i] * RetVal[2 * (2 * i - 1)] - c[2 * i + 1] * RetVal[2 * (2 * i - 1) + 1];
        RetVal[2 * (2 * i + 1) + 1] += c[2 * i] * RetVal[2 * (2 * i - 1) + 1] + c[2 * i + 1] * RetVal[2 * (2 * i - 1)];

        for (j = 2 * i; j > 1; --j)
        {
            RetVal[2 * j] += b[2 * i] * RetVal[2 * (j - 1)] - b[2 * i + 1] * RetVal[2 * (j - 1) + 1] +
                c[2 * i] * RetVal[2 * (j - 2)] - c[2 * i + 1] * RetVal[2 * (j - 2) + 1];
            RetVal[2 * j + 1] += b[2 * i] * RetVal[2 * (j - 1) + 1] + b[2 * i + 1] * RetVal[2 * (j - 1)] +
                c[2 * i] * RetVal[2 * (j - 2) + 1] + c[2 * i + 1] * RetVal[2 * (j - 2)];
        }

        RetVal[2] += b[2 * i] * RetVal[0] - b[2 * i + 1] * RetVal[1] + c[2 * i];
        RetVal[3] += b[2 * i] * RetVal[1] + b[2 * i + 1] * RetVal[0] + c[2 * i + 1];
        RetVal[0] += b[2 * i];
        RetVal[1] += b[2 * i + 1];
    }

    return RetVal;
}

std::vector<double> ButterBP::ComputeNumCoeffs(double Lcutoff, double Ucutoff, std::vector<double> DenC)
{
    std::vector<double> TCoeffs;
    std::vector<double> NumCoeffs(2 * FILTERORDER + 1);
    std::vector<std::complex<double>> NormalizedKernel(2 * FILTERORDER + 1);

    std::vector<double> Numbers;
    for (double n = 0; n < FILTERORDER * 2 + 1; n++)
        Numbers.push_back(n);
    int i;

    TCoeffs = ComputeHP();

    for (i = 0; i < FILTERORDER; ++i)
    {
        NumCoeffs[2 * i] = TCoeffs[i];
        NumCoeffs[2 * i + 1] = 0.0;
    }
    NumCoeffs[2 * FILTERORDER] = TCoeffs[FILTERORDER];

    double cp[2];
    double Bw, Wn;
    cp[0] = 2 * 2.0*tan(M_PI * Lcutoff / 2.0);
    cp[1] = 2 * 2.0*tan(M_PI * Ucutoff / 2.0);

    Bw = cp[1] - cp[0];
    //center frequency
    Wn = sqrt(cp[0] * cp[1]);
    Wn = 2 * atan2(Wn, 4);
    const std::complex<double> result = std::complex<double>(-1, 0);

    for (int k = 0; k< FILTERORDER * 2 + 1; k++)
    {
        NormalizedKernel[k] = std::exp(-sqrt(result)*Wn*Numbers[k]);
    }
    double b = 0;
    double den = 0;
    for (int d = 0; d < FILTERORDER * 2 + 1; d++)
    {
        b += real(NormalizedKernel[d] * NumCoeffs[d]);
        den += real(NormalizedKernel[d] * DenC[d]);
    }
    for (int c = 0; c < FILTERORDER * 2 + 1; c++)
    {
        NumCoeffs[c] = (NumCoeffs[c] * den) / b;
    }

    for (int i = NumCoeffs.size() - 1; i > FILTERORDER * 2 + 1; i--)
        NumCoeffs.pop_back();

    return NumCoeffs;
}

std::vector<double> ButterBP::ComputeLP()
{
    std::vector<double> NumCoeffs(FILTERORDER + 1);
    int m;
    int i;

    NumCoeffs[0] = 1;
    NumCoeffs[1] = FILTERORDER;
    m = FILTERORDER / 2;
    for (i = 2; i <= m; ++i)
    {
        NumCoeffs[i] = (double)(FILTERORDER - i + 1)*NumCoeffs[i - 1] / i;
        NumCoeffs[FILTERORDER - i] = NumCoeffs[i];
    }
    NumCoeffs[FILTERORDER - 1] = FILTERORDER;
    NumCoeffs[FILTERORDER] = 1;

    return NumCoeffs;
}

std::vector<double> ButterBP::ComputeHP()
{
    std::vector<double> NumCoeffs;
    int i;

    NumCoeffs = ComputeLP();

    for (i = 0; i <= FILTERORDER; ++i)
        if (i % 2) NumCoeffs[i] = -NumCoeffs[i];

    return NumCoeffs;
}
