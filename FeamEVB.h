//
// Unpacking FEAM data
// K.Olchanski
//

#ifndef FeamEVB_H
#define FeamEVB_H

#include "Feam.h"
#include "TsSync.h"

#include <vector>
#include <deque>

class FeamEVB
{
 public:
   unsigned fNumModules;
   double fEpsSec;
   TsSync fSync;
   int fCounter;
   std::vector<FeamAsm*> fAsm;
   std::vector<FeamModuleData*> fData;
   std::deque<FeamModuleData*> fBuf;
   std::deque<FeamEvent*> fEvents;
   double fMaxDt;
   double fMinDt;
   int fCountComplete;
   int fCountIncomplete;
   int fCountDuplicate;
   int fCountError;
   int fCountDropped;

 public:
   FeamEVB(int num_modules, double ts_freq, double eps_sec); // ctor
   ~FeamEVB(); // dtor

   FeamEvent* FindEvent(double t);
   void CheckFeam(FeamEvent *e);
   void AddFeam(int position, FeamModuleData *m);
   void Build();
   void BuildLastEvent();
   void AddPacket(const char* bank, int position, const FeamPacket* p, const char* ptr, int size);
   void Finalize(int position);
   void Print() const;
   FeamEvent* Get();
   FeamEvent* GetLastEvent();
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


