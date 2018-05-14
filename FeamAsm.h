//
// Unpacking FEAM data
// K.Olchanski
//

#ifndef FeamAsm_H
#define FeamAsm_H

#include "Feam.h"

#include <vector>

class FeamAsm
{
 public: // config

 public: // state
   int fCounter = 0;
   std::vector<FeamModuleAsm*> fAsm;
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
   FeamAsm(); // ctor
   ~FeamAsm(); // dtor

   //FeamEvent* FindEvent(double t);
   void CheckFeam(FeamEvent *e);
   void AddFeam(int position, FeamModuleData *m);
   void Build(bool force_build = false);
   void BuildLastEvent();
   void AddPacket(int imodule, int icolumn, int iring, int format, const FeamPacket* p, const char* ptr, int size);
   void Flush(int position);
   void Finalize(int position, FeamModuleData* m);
   void Print() const;
   //FeamEvent* Get();
   //FeamEvent* GetLastEvent();
   void BuildEvent(FeamEvent* event);
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


