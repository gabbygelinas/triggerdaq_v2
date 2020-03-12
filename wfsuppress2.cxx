// wfsuppress.cxx - waveform suppression

#include "wfsuppress2.h"
#include <stdio.h>

WfSuppress2::WfSuppress2() // ctor
{

}
WfSuppress2::~WfSuppress2() // dtor
{

}

void WfSuppress2::Init(int16_t a, int16_t t)
{
   fCounter = 0;
   fSum = 0;

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

bool WfSuppress2::Add(int16_t a)
{
   // this clock

   fCounter++;
   fSum += a;

   if (fCounter == 1)
      abase = fSum;
   else if (fCounter == 2)
      abase = fSum>>1; // div by 2
   else if (fCounter == 4)
      abase = fSum>>2; // div by 4
   else if (fCounter == 8)
      abase = fSum>>3; // div by 8
   else if (fCounter == 16)
      abase = fSum>>4; // div by 16
   else if (fCounter == 32)
      abase = fSum>>5; // div by 32
   else if (fCounter == 64)
      abase = fSum>>6; // div by 64
   else if (fCounter == 128)
      abase = fSum>>7; // div by 128

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
   
   return above||below||xclipped;
}

void WfSuppress2::Print() const
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
