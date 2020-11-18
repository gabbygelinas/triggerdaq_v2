// wfsuppress_adc.cxx - waveform suppression as implemented in ADC firmware

#include "wfsuppress_adc.h"
#include <stdio.h>

WfSuppressAdc::WfSuppressAdc() // ctor
{

}
WfSuppressAdc::~WfSuppressAdc() // dtor
{

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
}

bool WfSuppressAdc::Add(int adc_stream)
{
   fCounter++;

   if (fCounter < 5) {
      fAdcValue = 0;
   } else {
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
   }

   
   //bool xclipped = false; // (a == -2048) || (a == 2047);
   fTrigPos = (fAdcValue >= fThreshold);
   fTrigNeg = (fAdcValue <= -fThreshold);
   fTrig = fTrigPos | fTrigNeg;
   
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
