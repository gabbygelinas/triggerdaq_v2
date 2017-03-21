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
   std::deque<FeamModuleData*> fBuf;
   std::deque<AgEvent*> fEvents;

   AgEVB(double a16_ts_freq, double feam_ts_freq); // ctor
   AgEvent* FindEvent(double t);
   void CheckEvent(AgEvent *e);
   void AddAlpha16(Alpha16Event *e);
   void AddFeam(FeamEvent *e);
   void Build();
   void AddPacket(int ifeam, const FeamPacket* p, const char* ptr, int size);
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


