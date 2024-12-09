//
//  atans.hpp
//  VibratoTransfer
//
//  Created by Jeremy Hyrkas on 11/27/24.
//

/*inline float atan_approx1(float x)
{
  float a1  =  0.99997726f;
  float a3  = -0.33262347f;
  float a5  =  0.19354346f;
  float a7  = -0.11643287f;
  float a9  =  0.05265332f;
  float a11 = -0.01172120f;


  float x_sq = x*x;
  return
    x * (a1 + x_sq * (a3 + x_sq * (a5 + x_sq * (a7 + x_sq * (a9 + x_sq * a11)))));
}*/
inline float atan_approx2(float x)
{


   float s = x*x;
   return
  ((-0.0464964749f * s + 0.15931422f) * s - 0.327622764f) * s * x + x;
}

inline float fatan2f(float y, float x)
{
    int32_t swap = fabsf(x) < fabsf(y);
    float atan_input = (swap ? x : y) / (swap ? y : x);


    // Approximate atan 5000
    float res = atan_approx2(atan_input);
    // If swapped, adjust atan output 34000
    res = swap ? (atan_input >= 0.0f ? 1.570796326794897f : -1.570796326794897f) - res : res;
    // Adjust quadrants
    if (x < 0.0f)
      res += (y >= 0.0f ? 3.141592653589793f : -3.141592653589793f);
    return(res);
}
