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
#include "Feam.h"
//#include "FeamAsm.h"
#include "Tdc.h"
#include "ncfm.h"
#include "Sim.h"

class AgAsm
{
 public: // settings
   double fConfMaxDt = 0.000000100; // max timestamp deviation in sec
   int fConfAdc32Rev = -1; // fmc-adc32 board revision number

 public: // event builder state
   int    fCounter = 0;
   double fLastEventTime = 0;

 public: // diagnostics
   double fTrgMaxDt = 0;
   double fAdcMaxDt = 0;
   double fPwbMaxDt = 0;
   double fTdcMaxDt = 0;

 public: // counters
   int fCountComplete   = 0;
   int fCountCompleteWithError = 0;
   int fCountIncomplete = 0;
   int fCountIncompleteWithError = 0;

   int fCountMissingTrg = 0;
   int fCountMissingAdc = 0;
   int fCountMissingPwb = 0;
   int fCountMissingTdc = 0;

   int fCountErrorTrg = 0;
   int fCountErrorAdc = 0;
   int fCountErrorPwb = 0;
   int fCountErrorTdc = 0;

   int fCountLoneTdc = 0;

 public: // count missing triggers
   int fCountMissingTrgTrig = 0;
   int fCountMissingAdcTrig = 0;
   int fCountMissingPwbTrig = 0;
   int fCountMissingTdcTrig = 0;

 public: // member functions
   AgAsm(); // ctor
   ~AgAsm(); // dtor
   void Reset();
   void BeginRun(int runno);
   AgEvent* UnpackEvent(TMEvent* me);
   void Print() const;

 public: // internal data
   Ncfm* fCfm = NULL;
   TrgAsm* fTrgAsm = NULL;
   std::vector<std::string> fAdcMap;
   std::vector<std::string> fPwbBanks;
   Alpha16Asm* fAdcAsm = NULL;
   PwbModuleMap* fPwbModuleMap = NULL;
   PwbAsm* fPwbAsm = NULL;
   //FeamAsm* fFeamAsm = NULL;
   TdcAsm* fTdcAsm = NULL;
   SimAsm* fSimAsm = NULL;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


