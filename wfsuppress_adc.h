// wfsuppress_adc.h - waveform suppression as implemented in the ADC firmware

#include <string>

class WfSuppressAdc
{
 public:
   int fThreshold = 0;
   int fCounter = 0;
   int fBaselineCounter = 0;
   int fBaselineSum = 0;
   bool fBaselineReady = false;
   int fBaseline = 0;
   int fAdcValue = 0;
   bool fTrigPos = false;
   bool fTrigNeg = false;
   bool fTrig = false;
   
 public:
   WfSuppressAdc(); // ctor
   ~WfSuppressAdc(); // dtor
   
 public:
   void Reset();
   bool Add(int adc_stream);
   std::string PrintToString() const;
};

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
