// wfsuppress_pwb.cxx - waveform suppression as implemented in PWB firmware

#include "wfsuppress_pwb.h"
#include <stdio.h>

WfSuppressPwb::WfSuppressPwb() // ctor
{

}
WfSuppressPwb::~WfSuppressPwb() // dtor
{

}

void WfSuppressPwb::Reset()
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

bool WfSuppressPwb::Add(int adc_stream)
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

std::string WfSuppressPwb::PrintToString() const
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
