// unpack_cb.h

#ifndef INCLUDE_UNPACK_CB_H
#define INCLUDE_UNPACK_CB_H

#include <cstddef> // size_t
#include <stdint.h> // uint32_t
#include <vector> // std::vector

#define CB_HIT_FLAG_TE (1<<0)

struct CbHit
{
   double   time = 0;
   uint32_t epoch = 0;
   uint32_t timestamp = 0;
   uint32_t channel = 0;
   uint32_t flags = 0;
};

typedef std::vector<CbHit> CbHits;

struct CbScaler
{
   uint32_t raw = 0;
   double   epoch = 0;
   double   rate = 0;
   double   sum = 0;
};

typedef std::vector<CbScaler> CbLatchedScalers;

typedef std::vector<CbLatchedScalers> CbScalers;

class CbUnpack
{
 public:
   bool fVerbose = false;
   double   fCbTsFreq = 10*1e6; // 10 MHz
   int fKludge = 0;
   //  fKludge value 0 is chronobox fw 0x62608957
   //  fKludge value 1 is cbtrg fw 0x618b790b

 public:
   CbUnpack(uint32_t num_inputs); // ctor
   ~CbUnpack(); // dtor

 public:
   void Unpack(const uint32_t* data, size_t nwords, CbHits* hits, CbScalers* scalers);

 public: // public data
   bool     fFailed = false;
   uint32_t fNumScalers = 0;
   uint32_t fNumInputs = 0;
   bool     fWaitForEpoch0 = true;
   
 public: // internal state: scalers
   bool fInScalersPacket = false;
   std::vector<uint32_t> fScalers;
   std::vector<uint32_t> fScalersPrev;
   std::vector<uint32_t> fScalersIncr;
   std::vector<double> fScalersEpoch;

 public: // internal state: timestamps
   std::vector<uint32_t> fChanLastTimestamp;
   std::vector<uint32_t> fChanEpoch;
   uint32_t fCurrentEpoch = 0;
   uint32_t fCurrentTsRangeMin = 0;
   uint32_t fCurrentTsRangeMax = 0;
   bool fWaitingForData = true;

 public: // internal functions
   void SaveScalers(CbScalers* scalers);
};

void PrintCbHits(const CbHits& hits);
void PrintCbScalers(const CbScalers& scalers);

#endif
/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
