//
// unpack_cb_module.cxx
//
// K.Olchanski
//

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "manalyzer.h"
#include "midasio.h"

#include "AgFlow.h"
#include "unpack_cb.h"

#include "ncfm.h"

#ifdef _TIME_ANALYSIS_
#include "AnalysisTimer.h"
#endif

class UnpackFlags
{
public:
   bool fPrint = false;
   bool fVerboseCbtrg = false;
   bool fVerboseCb01 = false;
   bool fVerboseCb02 = false;
   bool fVerboseCb03 = false;
   bool fVerboseCb04 = false;
};

class UnpackCbModule: public TARunObject
{
public:
   UnpackFlags* fFlags = NULL;

   std::vector<std::string> fCbBanks;
   std::vector<CbUnpack*>   fCbUnpack;

   bool fTrace = false;
   
   UnpackCbModule(TARunInfo* runinfo, UnpackFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("UnpackCbModule::ctor!\n");

      fModuleName = "unpack_cb_module";
      fFlags   = flags;
   }

   ~UnpackCbModule()
   {
      if (fTrace)
         printf("UnpackCbModule::dtor!\n");

      fCbBanks.clear();
      for (size_t i=0; i<fCbUnpack.size(); i++) {
         if (fCbUnpack[i]) {
            delete fCbUnpack[i];
            fCbUnpack[i] = NULL;
         }
      }
      fCbUnpack.clear();
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackCbModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      fCbBanks.push_back("CBFT");
      fCbBanks.push_back("CBF1");
      fCbBanks.push_back("CBF2");
      fCbBanks.push_back("CBF3");
      fCbBanks.push_back("CBF4");

      fCbUnpack.push_back(new CbUnpack(23));
      fCbUnpack.push_back(new CbUnpack(59));
      fCbUnpack.push_back(new CbUnpack(59));
      fCbUnpack.push_back(new CbUnpack(59));
      fCbUnpack.push_back(new CbUnpack(59));

      fCbUnpack[0]->fKludge = 1;
   }

   void PreEndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackCbModule::PreEndRun, run %d\n", runinfo->fRunNo);
   }
   
   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackCbModule::EndRun, run %d\n", runinfo->fRunNo);
   }
   
   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackCbModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackCbModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void UnpackCbFifo(TARunInfo* runinfo, int ibank, CbUnpack* cb, TMEvent* event, TMBank* cbbank)
   {
      int nwords = cbbank->data_size/4; // byte to uint32_t words
      uint32_t* cbdata = (uint32_t*)event->GetBankData(cbbank);
      
      CbHitsFlow* hits_flow = new CbHitsFlow(NULL);
      CbScalers scalers;

      hits_flow->fCbBankName = fCbBanks[ibank];
      hits_flow->fCbIndex = ibank;
      hits_flow->fNumInputs = cb->fNumInputs;
      
      //if (ibank==0)
      //   cb->fVerbose = true;
      //else
      //   cb->fVerbose = false;
      //
      //printf("%s: %d words\n", fCbBanks[ibank].c_str(), nwords);

      if (ibank==0 && fFlags->fVerboseCbtrg)
         cb->fVerbose = true;
      if (ibank==1 && fFlags->fVerboseCb01)
         cb->fVerbose = true;
      if (ibank==2 && fFlags->fVerboseCb02)
         cb->fVerbose = true;
      if (ibank==3 && fFlags->fVerboseCb03)
         cb->fVerbose = true;
      if (ibank==4 && fFlags->fVerboseCb04)
         cb->fVerbose = true;

      cb->Unpack(cbdata, nwords, &hits_flow->fHits, &scalers);

      if (fFlags->fPrint) {
         printf("%s: %5d words, %5zu hits, %1zu scalers\n", fCbBanks[ibank].c_str(), nwords, hits_flow->fHits.size(), scalers.size());
      }
      
      //if (hits.size() > 0) {
      //   if (0) {
      //      printf("Data from %s: ", name);
      //      PrintCbHits(hits);
      //   }
      //}
      
      //if (ibank == 0) {
      //for (size_t i=0; i<hits.size(); i++) {
      //   printf("%s: hit %zu, time %.6f sec, %d+%d, channel %2d (%d)\n", cbbank->name.c_str(), i, hits[i].time, hits[i].timestamp, hits[i].epoch, hits[i].channel, (hits[i].flags&CB_HIT_FLAG_TE));
      //}
      //}
      
      //if (scalers.size() > 0) {
      //   if (print) {
      //      printf("Data from %s: ", name);
      //      PrintCbScalers(scalers);
      //   }
      //}

      if (hits_flow->fHits.empty())
         delete hits_flow;
      else
         runinfo->AddToFlowQueue(hits_flow);

      // bank of chronobox data can contain more than one block of scalers,
      // each block goes into it's own flow event.

      for (size_t i=0; i<scalers.size(); i++) {
         CbScalersFlow* scalers_flow = new CbScalersFlow(NULL);
         scalers_flow->fCbBankName = fCbBanks[ibank];
         scalers_flow->fCbIndex = ibank;
         scalers_flow->fNumInputs = cb->fNumInputs;
         scalers_flow->fScalers = scalers[i];
         runinfo->AddToFlowQueue(scalers_flow);
      }
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id == 1) {
         static bool once = false;
         if (!once)
            printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
         once = true;
      }

      if (event->event_id != 4)
         return flow;

      printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      for (size_t ibank=0; ibank<fCbBanks.size(); ibank++) {
         TMBank* cbbank = event->FindBank(fCbBanks[ibank].c_str());
         if (cbbank) {
            UnpackCbFifo(runinfo, ibank, fCbUnpack[ibank], event, cbbank);
         }
      }

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("UnpackCbModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class UnpackCbModuleFactory: public TAFactory
{
public:
   UnpackFlags fFlags;

public:
   void Usage()
   {
      printf("UnpackCbModuleFactory flags:\n");
      printf("--print-cb -- print chronobox unpacked data\n");
      printf("--verbose-cbtrg -- print chronobox cbtrg raw data\n");
      printf("--verbose-cb01  -- print chronobox cb01 raw data\n");
      printf("--verbose-cb02  -- print chronobox cb02 raw data\n");
      printf("--verbose-cb03  -- print chronobox cb03 raw data\n");
      printf("--verbose-cb04  -- print chronobox cb04 raw data\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("UnpackCbModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--print-cb") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--verbose-cbtrg") {
            fFlags.fVerboseCbtrg = true;
         }
         if (args[i] == "--verbose-cb01") {
            fFlags.fVerboseCb01 = true;
         }
         if (args[i] == "--verbose-cb02") {
            fFlags.fVerboseCb02 = true;
         }
         if (args[i] == "--verbose-cb03") {
            fFlags.fVerboseCb03 = true;
         }
         if (args[i] == "--verbose-cb04") {
            fFlags.fVerboseCb04 = true;
         }
      }
   }

   void Finish()
   {
      printf("UnpackCbModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("UnpackCbModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new UnpackCbModule(runinfo, &fFlags);
   }
};

static TARegister tar(new UnpackCbModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
