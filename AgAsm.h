//
// AgAsm.h
//
// ALPHA-g event assembler
// K.Olchanski
//

#ifndef AgAsm_H
#define AgAsm_H

#include "midasio.h" // TMEvent
#include "AgEvent.h"
#include "TrgAsm.h"
#include "PwbAsm.h"
#include "FeamAsm.h"

class AgAsm
{
 public: // settings

 public: // event builder state
   int    fCounter;
   double fLastEventTime = 0;

 public: // diagnostics

 public: // counters
   int fCount = 0;
   int fCountComplete   = 0;
   int fCountError      = 0;
   int fCountIncomplete = 0;

 public: // member functions
   AgAsm(); // ctor
   ~AgAsm(); // dtor
   AgEvent* UnpackEvent(TMEvent* me);

 public: // internal data
   TrgAsm* fTrgAsm = NULL;
   std::vector<std::string> fAdcMap;
   Alpha16Asm* fAdcAsm = NULL;
   PwbModuleMap* fPwbMap = NULL;
   PwbAsm* fPwbAsm = NULL;
   FeamAsm* fFeamAsm = NULL;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


