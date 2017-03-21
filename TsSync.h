//
// TsSync.h
// Timestamp synchronization
// K.Olchanski
//

#ifndef TsSyncH
#define TsSyncH

#include <stdio.h>
#include <stdint.h>
#include <vector>

struct TsSyncEntry
{
   uint32_t ts;
   int epoch;
   double time;

   TsSyncEntry(uint32_t xts, int xepoch, double xtime); // ctor
};

class TsSyncModule
{
public:
   double   fFreqHz;
   int      fEpoch;
   uint32_t fFirstTs;
   uint32_t fPrevTs;
   uint32_t fLastTs;
   double   fOffsetSec;
   double   fPrevTimeSec;
   double   fLastTimeSec;
   double   fEps;
   int      fSyncedWith;
   bool     fOverflow;
   unsigned fBufMax;

   std::vector<TsSyncEntry> fBuf;

public:
   TsSyncModule(); // ctor
   void Print() const;
   double GetTime(uint32_t ts, int epoch) const;
   void Add(uint32_t ts);
   void Retime();
   double GetDt(unsigned j);
   unsigned FindDt(double dt);
};

class TsSync
{
public:
   std::vector<TsSyncModule> fModules;
   bool fSyncOk;
   bool fOverflow;
   bool fTrace;

public:
   TsSync(); // ctor
   ~TsSync(); // dtor
   void Configure(unsigned i, double freq_hz, int buf_max);
   void CheckSync(unsigned ii, unsigned i);
   void Check(unsigned inew);
   void Add(unsigned i, uint32_t ts);
   void Dump() const;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
