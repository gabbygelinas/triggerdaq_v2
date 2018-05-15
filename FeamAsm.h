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
   double fConfMaxDt = 0.000000100; // 100 ns

 public: // state
   int fCounter = 0;
   std::vector<FeamModuleAsm*> fAsm;

 public: // diagnostics and counters
   int fCountComplete = 0;
   int fCountIncomplete = 0;
   int fCountError = 0;
   double fMaxDt = 0;

 public:
   FeamAsm(); // ctor
   ~FeamAsm(); // dtor

   void AddPacket(int imodule, int icolumn, int iring, int format, const FeamPacket* p, const char* ptr, int size);
   void Print() const;
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


