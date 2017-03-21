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

class AgEVB
{
public:
   TsSync fSync;
   int fCounter;
   std::deque<Alpha16Event*> fBuf0;
   std::deque<FeamEvent*> fBuf1;
   std::deque<AgEvent*> fEvents;

   AgEVB(double a16_ts_freq, double feam_ts_freq); // ctor
   void AddAlpha16Event(Alpha16Event *e);
   void AddFeamEvent(FeamEvent *e);
   AgEvent* FindEvent(double t);
   void CheckEvent(AgEvent *e);
   void BuildAlpha16(Alpha16Event *e);
   void BuildFeam(FeamEvent *e);
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


