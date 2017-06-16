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
 public: // settings
   unsigned fMaxSkew;
   unsigned fMaxDead;
   double   fEpsSec;
   bool     fClockDrift;

 public: // event builder state
   TsSync fSync;
   int    fCounter;
   std::deque<AgEvbBuf*> fBuf[2];
   std::deque<AgEvent*> fEvents;
   double fLastA16Time;
   double fLastFeamTime;

 public: // diagnostics
   double fMaxDt;
   double fMinDt;

 public: // counters
   int fCount = 0;
   int fCountComplete   = 0;
   int fCountError      = 0;
   int fCountIncomplete = 0;
   int fCountIncompleteA16  = 0;
   int fCountIncompleteFeam = 0;
   int fCountA16 = 0;
   int fCountFeam = 0;
   int fCountRejectedA16 = 0;
   int fCountRejectedFeam = 0;

 public: // member functions
   AgEVB(double a16_ts_freq, double feam_ts_freq, double eps_sec, int max_skew, int max_dead, bool clock_drift); // ctor
   ~AgEVB(); // dtor
   void AddAlpha16Event(Alpha16Event *e);
   void AddFeamEvent(FeamEvent *e);
   AgEvent* FindEvent(double t);
   void CheckEvent(AgEvent *e);
   void Build(int index, AgEvbBuf *m);
   void Build();
   void Print() const;
   void UpdateCounters(const AgEvent* e);
   AgEvent* Get();
   AgEvent* GetLastEvent();
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


