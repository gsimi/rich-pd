
#define EXP_A (1048576/0.693147180559945286)

#define EXP_C 60801

inline double fastexp(double y)
{
  union
  {
    double d;
    //#ifdef LITTLE_ENDIAN
    struct { int j, i; } n;
    //#elseif
    //    struct { int i, j; } n;
    //#endif
  }
  _eco;
  
  _eco.n.i = (int)(EXP_A*(y)) + (1072693248 - EXP_C);
  _eco.n.j = 0;

  return _eco.d;
}


Double_t fastBesselI1(Double_t x)
{
   // Compute the modified Bessel function I_1(x) for any real x.
   //
   //  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   //     Applied Mathematics Series vol. 55 (1964), Washington.
   //
   //--- NvE 12-mar-2000 UU-SAP Utrecht

   // Parameters of the polynomial approximation
   const Double_t p1=0.5,          p2=0.87890594,   p3=0.51498869,
                  p4=0.15084934,   p5=2.658733e-2,  p6=3.01532e-3,  p7=3.2411e-4;

   const Double_t q1= 0.39894228,  q2=-3.988024e-2, q3=-3.62018e-3,
                  q4= 1.63801e-3,  q5=-1.031555e-2, q6= 2.282967e-2,
                  q7=-2.895312e-2, q8= 1.787654e-2, q9=-4.20059e-3;

   const Double_t k1 = 3.75;
   Double_t ax = TMath::Abs(x);

   Double_t y=0, result=0;

   if (ax < k1) {
      Double_t xx = x/k1;
      y = xx*xx;
      result = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = k1/ax;
      result = (fastexp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
      if (x < 0) result = -result;
   }
   return result;
}
