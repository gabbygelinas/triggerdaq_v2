// wfsuppress_adc.h - waveform suppression as implemented in the ADC firmware

#include <string>

class WfSuppressAdc
{
 public: // configuration
   int fThreshold = 0;
   int fKeepMore  = 0;
   int fNumSamples = 0;

 public: // internal state
   int fCounter   = 0;
   int fBaselineCounter = 0;
   int fBaselineSum = 0;

 public: // output
   bool fBaselineReady = false;
   int fBaseline = 0;
   int fAdcValue = 0;
   bool fTrigPos = false;
   bool fTrigNeg = false;
   bool fTrig = false;
   bool fKeepBit = false;
   int  fKeepLast = 0;
   int  fAdcLast  = 0; // ADC value corresponding to fKeepLast

   int  fAdcMaxPos = 0; // biggest ADC value, positive polarity
   int  fAdcMaxNeg = 0; // biggest ADC value, negative polarity
   
 public:
   WfSuppressAdc(); // ctor
   ~WfSuppressAdc(); // dtor
   
 public:
   void Config(int threshold, int keep_more, int nsamples);
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
