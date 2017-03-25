//
// AgEVB.h
//
// ALPHA-g event builder
// K.Olchanski
//

#ifndef AgEVB_H
#define AgEVB_H

#include "AgEvent.h"
#include "TsSync.h"

#include <vector>
#include <deque>

struct AgEvbBuf
{
   Alpha16Event* a16;
   FeamEvent* feam;

   uint32_t ts;
   int epoch;
   double time;
   double timeIncr;
};

class AgEVB
{
public:
   TsSync fSync;
   int    fCounter;
   std::deque<AgEvbBuf*> fBuf[2];
   std::deque<AgEvent*> fEvents;
   unsigned fMaxSkew;
   unsigned fMaxDead;
   double fLastA16Time;
   double fLastFeamTime;
   double fMaxDt;
   double fMinDt;
   double fEpsSec;

   AgEVB(double a16_ts_freq, double feam_ts_freq, double eps_sec, int max_skew, int max_dead); // ctor
   ~AgEVB(); // dtor
   void AddAlpha16Event(Alpha16Event *e);
   void AddFeamEvent(FeamEvent *e);
   AgEvent* FindEvent(double t);
   void CheckEvent(AgEvent *e);
   void Build(int index, AgEvbBuf *m);
   void Build();
   void Print() const;
   AgEvent* Get();
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


