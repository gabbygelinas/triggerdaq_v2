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
   double fConfMaxDt = 0.000000100; // max timestamp deviation in sec

 public: // event builder state
   int    fCounter = 0;
   double fLastEventTime = 0;

 public: // diagnostics
   double fTrgMaxDt = 0;
   double fAdcMaxDt = 0;
   double fPwbMaxDt = 0;

 public: // counters
   int fCountComplete   = 0;
   int fCountCompleteWithError = 0;
   int fCountIncomplete = 0;
   int fCountIncompleteWithError = 0;

 public: // member functions
   AgAsm(); // ctor
   ~AgAsm(); // dtor
   AgEvent* UnpackEvent(TMEvent* me);
   void Print() const;

 public: // internal data
   TrgAsm* fTrgAsm = NULL;
   std::vector<std::string> fAdcMap;
   std::vector<std::string> fFeamBanks;
   Alpha16Asm* fAdcAsm = NULL;
   PwbModuleMap* fPwbModuleMap = NULL;
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


