// unpack_cb.h

#ifndef INCLUDE_UNPACK_CB_H
#define INCLUDE_UNPACK_CB_H

#include <cstddef> // size_t
#include <stdint.h> // uint32_t
#include <vector> // std::vector

#define CB_HIT_FLAG_TE (1<<0)

struct CbHit
{
   double   time;
   uint32_t epoch;
   uint32_t timestamp;
   uint32_t channel;
   uint32_t flags;
};

typedef std::vector<CbHit> CbHits;

struct CbScaler
{
   uint32_t raw;
   double   epoch;
   double   rate;
   double   sum;
};

typedef std::vector<CbScaler> CbLatchedScalers;

typedef std::vector<CbLatchedScalers> CbScalers;

class CbUnpack
{
 public:
   bool fVerbose = false;
   double   fCbTsFreq = 10*1e6; // 10 MHz
   //bool     fCbEpochFromReset = true;

 public:
   CbUnpack(uint32_t num_inputs); // ctor
   ~CbUnpack(); // dtor

 public:
   void Unpack(const uint32_t* data, size_t nwords, CbHits* hits, CbScalers* scalers, int kluge_offset = 0);

 public: // public data
   bool     fFailed = false;
   uint32_t fNumScalers = 0;
   uint32_t fNumInputs = 0;
   std::vector<uint32_t> fChanLastTimestamp;
   std::vector<uint32_t> fChanEpoch;
   //uint32_t fEpoch = 0;
   //uint32_t fEpochCounter = 0;
   
 public: // internal state
   bool fInScalersPacket = false;
   std::vector<uint32_t> fScalers;
   std::vector<uint32_t> fScalersPrev;
   std::vector<uint32_t> fScalersIncr;
   std::vector<double> fScalersEpoch;

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
