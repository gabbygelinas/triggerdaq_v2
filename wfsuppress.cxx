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
  a0 = a;
  a1 = a;
  a2 = a;
  a3 = a;
  a4 = a;
  a5 = a;
  a6 = a;
  a7 = a;
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
  abase = (uint16_t)(((uint32_t)a0 + (uint32_t)a1 + (uint32_t)a2 + (uint32_t)a3 + (uint32_t)a4 + (uint32_t)a5 + (uint32_t)a6 + (uint32_t)a7) >> (uint32_t)3);

  //printf("%d %d %d %d %d %d %d %d, abase %d\n", a0, a1, a2, a3, a4, a5, a6, a7, abase);

  //a &= 0xFFF; // 12-bit of ADC data

  amp = a - abase;

  bool xclipped = (a == -2048) || (a == 2047);
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
  a3 = a4;
  a4 = a5;
  a5 = a6;
  a6 = a7;
  a7 = a;

  return above||below||xclipped;
}

void WfSuppress::Print() const
{
  printf("thr %d, adc %d..%d, base %4d, amp %4d..%4d, clip %d keep %d %d %d", threshold, amin, amax, abase, ampMin, ampMax, clipped, below_threshold, above_threshold, keep);
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
