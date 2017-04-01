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
   TsSync fSync;
   int fCounter;
   std::vector<FeamModuleData*> fData;
   std::deque<FeamModuleData*> fBuf;
   std::deque<FeamEvent*> fEvents;

   FeamEVB(int num_modules, double ts_freq); // ctor
   FeamEvent* FindEvent(double t);
   void CheckFeam(FeamEvent *e);
   void AddFeam(int position, FeamModuleData *m);
   void Build();
   void AddPacket(const char* bank, int position, const FeamPacket* p, const char* ptr, int size);
   void Print() const;
   FeamEvent* Get();
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


