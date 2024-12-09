/*
  ==============================================================================

    ButterBP.h
    Created: 9 Dec 2024 9:12:40am
    Author:  Jeremy Hyrkas

  ==============================================================================
*/

#ifndef ButterBP_hpp
#define ButterBP_hpp

#include<complex>
#include<vector>

#define FILTERORDER 2

class ButterBP
{
public:
    ButterBP(float bpL, float bpH, float fs);
    virtual ~ButterBP();
    virtual float processSample(float sample);
    virtual void setParams(float bpL, float bpH, float fs);
protected:
    // making these doubles because we aren't using a SOS...hopefully this helps
    double b0, b1, b2, b3, b4, a1, a2, a3, a4;
    // TODO: reimplement in canonical form for fewer floats kept around
    double b_zi[4];
    double a_zi[4];
    
    std::vector<double> ComputeDenCoeffs(double Lcutoff, double Ucutoff);

    std::vector<double> TrinomialMultiply(std::vector<double> b, std::vector<double> c);

    std::vector<double> ComputeNumCoeffs(double Lcutoff, double Ucutoff, std::vector<double> DenC);

    std::vector<double> ComputeLP();

    std::vector<double> ComputeHP();
};
#endif /* ButterBP_hpp */
