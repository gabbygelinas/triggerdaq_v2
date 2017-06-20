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
 public: // config
   unsigned fNumModules;
   double fEpsSec;

 public: // state
   TsSync fSync;
   int fCounter = 0;
   std::vector<FeamAsm*> fAsm;
   std::deque<FeamModuleData*> fBuf;
   std::deque<FeamEvent*> fEvents;

 public: // diagnostics and counters
   double fMaxDt = 0;
   double fMinDt = 0;
   int fCountComplete = 0;
   int fCountIncomplete = 0;
   int fCountDuplicate = 0;
   int fCountError = 0;

 public:
   FeamEVB(int num_modules, double ts_freq, double eps_sec); // ctor
   ~FeamEVB(); // dtor

   FeamEvent* FindEvent(double t);
   void CheckFeam(FeamEvent *e);
   void AddFeam(int position, FeamModuleData *m);
   void Build(bool force_build = false);
   void BuildLastEvent();
   void AddPacket(const char* bank, int position, const FeamPacket* p, const char* ptr, int size);
   void Flush(int position);
   void Finalize(int position, FeamModuleData* m);
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


