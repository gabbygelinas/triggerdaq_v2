// wfsuppress.cxx - waveform suppression

#include "wfsuppress.h"
#include <stdio.h>

WfSuppress::WfSuppress() // ctor
{

}
WfSuppress::~WfSuppress() // dtor
{

}

void WfSuppress::Init(int16_t a, int16_t t)
{
  aa = a;
  ab = a;
  ac = a;
  ad = a;
  a0 = a;
  a1 = a;
  a2 = a;
  a3 = a;
  a4 = a;
  a5 = a;
  a6 = a;
  a7 = a;
  a8 = a;
  a9 = a;
  a10 = a;
  a11 = a;
  a12 = a;
  a13 = a;
  a14 = a;
  a15 = a;
  abase = a;

  threshold = t;

  amp = 0;
  ampMin = 0;
  ampMax = 0;

  amin = a;
  amax = a;

  clipped = false;
  above_threshold = false;
  below_threshold = false;
  keep = false;
}

bool WfSuppress::Add(int16_t a)
{
  // this clock

  // NB: check for 16-bit overflow!!!

   //abase = (uint16_t)
   //  ((
   //    (uint32_t)a0 + (uint32_t)a1 + (uint32_t)a2 + (uint32_t)a3
   //    + (uint32_t)a4 + (uint32_t)a5 + (uint32_t)a6 + (uint32_t)a7
   //    + (uint32_t)a8 + (uint32_t)a9 + (uint32_t)a10 + (uint32_t)a11
   //    + (uint32_t)a12 + (uint32_t)a13 + (uint32_t)a14 + (uint32_t)a15
   //    ) >> (uint32_t)4);

   // 8 samples
   //abase = (uint16_t)(((uint32_t)a0 + (uint32_t)a1 + (uint32_t)a2 + (uint32_t)a3 + (uint32_t)a4 + (uint32_t)a5 + (uint32_t)a6 + (uint32_t)a7) >> (uint32_t)3);

   // 4 samples
   abase = (uint16_t)(((uint32_t)a0 + (uint32_t)a1 + (uint32_t)a2 + (uint32_t)a3) >> (uint32_t)2);

  //printf("%d %d %d %d %d %d %d %d, abase %d\n", a0, a1, a2, a3, a4, a5, a6, a7, abase);

  //a &= 0xFFF; // 12-bit of ADC data

#define MASK 0x000F

  //abase &= ~MASK;

  amp = a - abase;

  bool xclipped = false; // (a == -2048) || (a == 2047);
  bool above = (amp >= threshold);
  bool below = (amp <= -threshold);

  if (amp < ampMin)
    ampMin = amp;
  if (amp > ampMax)
    ampMax = amp;

  if (a < amin)
    amin = a;
  if (a > amax)
    amax = a;

  clipped |= xclipped;
  above_threshold |= above;
  below_threshold |= below;
  keep |= (above||below||xclipped);

  // next clock

  a0 = a1;
  a1 = a2;
  a2 = a3;
  a3 = a; // a4;
  a4 = a5;
  a5 = a6;
  a6 = a7;
  a7 = a8; // aa; // a;
  a8 = a1;
  a9 = a2;
  a10 = a3;
  a11 = a4;
  a12 = a5;
  a13 = a6;
  a14 = a7;
  a15 = aa; // a;

  aa = ab;
  ab = ac;
  ac = ad;
  ad = a; // & ~MASK;

  return above||below||xclipped;
}

void WfSuppress::Print() const
{
   printf("thr %d, adc %d..%d, range %d, base %4d, amp %4d..%4d, clip %d keep %d %d %d", threshold, amin, amax, amax-amin, abase, ampMin, ampMax, clipped, below_threshold, above_threshold, keep);
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
