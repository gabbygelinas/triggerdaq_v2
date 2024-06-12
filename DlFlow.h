//
// DlFlow.h
//
// manalyzer flow objects for DARK-LIGHT events
// K.Olchanski
//

#ifndef DlFlow_H
#define DlFlow_H

#include "unpack_cb.h"
#include "manalyzer.h"

class CbHitsFlow: public TAFlowEvent
{
public:
   std::string fCbBankName;
   int fCbIndex = 0; // 0=cbtrg, 1=cb01, etc
   int fNumInputs = 0;
   CbHits fHits;
   
public:
   CbHitsFlow(TAFlowEvent* flow) // ctor
      : TAFlowEvent(flow)
   {
   }
};

class CbScalersFlow: public TAFlowEvent
{
public:
   std::string fCbBankName;
   int fCbIndex = 0; // 0=cbtrg, 1=cb01, etc
   int fNumInputs = 0;
   CbLatchedScalers fScalers;
   
public:
   CbScalersFlow(TAFlowEvent* flow) // ctor
      : TAFlowEvent(flow)
   {
   }
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


