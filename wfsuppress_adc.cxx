// wfsuppress_adc.cxx - waveform suppression as implemented in ADC firmware

#include "wfsuppress_adc.h"
#include <stdio.h>

WfSuppressAdc::WfSuppressAdc() // ctor
{

}
WfSuppressAdc::~WfSuppressAdc() // dtor
{

}

void WfSuppressAdc::Config(int threshold, int keep_more, int nsamples)
{
   fThreshold = threshold;
   fKeepMore  = keep_more;
   fNumSamples = nsamples;
}

void WfSuppressAdc::Reset()
{
   fCounter = 0;

   fBaselineCounter = 0;
   fBaselineSum = 0;
   fBaseline = 0;
   fBaselineReady = false;

   fAdcValue = 0;

   fTrigPos = false;
   fTrigNeg = false;
   fTrig = false;

   fKeepBit = false;
   fKeepLast = 0;
   fAdcLast  = 0;
   fAdcMaxPos = 0;
   fAdcMaxNeg = 0;
}

bool WfSuppressAdc::Add(int adc_stream)
{
   fCounter++;

   if (fBaselineCounter < 64) {
      fBaselineSum += adc_stream;
      fBaselineCounter ++;
      if (fBaselineCounter == 64) {
         fBaseline = fBaselineSum/64;
         fBaselineReady = true;
      }
      fAdcValue = 0;
   } else {
      fAdcValue = adc_stream - fBaseline;
   }
   
   //bool xclipped = false; // (a == -2048) || (a == 2047);
   fTrigPos = (fAdcValue >= fThreshold);
   fTrigNeg = (fAdcValue <= -fThreshold);
   fTrig = fTrigPos | fTrigNeg;

   // data suppression does not see the first 64 or so samples
   // and it does not see the last few samples. KO 2020-NOV-19
   if (fBaselineReady && fCounter > 66 && fCounter <= fNumSamples - 6) {
      if (fAdcValue > fAdcMaxPos)
         fAdcMaxPos = fAdcValue;
      if (fAdcValue < fAdcMaxNeg)
         fAdcMaxNeg = fAdcValue;

      if (fTrig) {
         fKeepBit |= true;
         if (fCounter < fNumSamples - 7) {
            fKeepLast = (fCounter+1)/2 + 1;
            fAdcLast  = fAdcValue;
         }
      }
   }
   
   return fTrig;
}

std::string WfSuppressAdc::PrintToString() const
{
   char buf[1024];
   //printf("thr %d, adc %d..%d, range %d, base %4d, amp %4d..%4d, clip %d keep %d %d %d", fThreshold, amin, amax, amax-amin, abase, ampMin, ampMax, clipped, below_threshold, above_threshold, keep);
   sprintf(buf, "%3d: b %d %d %d, a %4d, t %d %d %d", fCounter, fBaselineCounter, fBaselineSum, fBaseline, fAdcValue, fTrigPos, fTrigNeg, fTrig);
   return buf;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
