//
// MIDAS analyzer for chronobox data
//
// K.Olchanski
//

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>

#include "manalyzer.h"
#include "midasio.h"

#include "unpack_cb.h"

class CbkoModule: public TARunObject
{
public:
   const double kTsFreq = 10.0e6; // 10 MHz

   uint32_t fT0 = 0;
   uint32_t fTE = 0;
   int fCountWC = 0;
   int fCountMissingWC = 0;

   std::vector<CbUnpack*> fCbUnpack;
   std::vector<std::string> fCbBanks;
   std::vector<std::vector<std::vector<CbHit>>> fCbHits;

   int fRunEventCounter = 0;

   bool fPrint = true;

   CbkoModule(TARunInfo* runinfo) : TARunObject(runinfo)
   {
      fModuleName="cbko_module";
      fRunEventCounter = 0;

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

      fCbHits.resize(fCbBanks.size());
      for (size_t i=0; i<fCbBanks.size(); i++) {
         fCbHits[i].resize(fCbUnpack[i]->fNumInputs);
      }
   }

   ~CbkoModule()
   {
      // empty
   }
  
   void BeginRun(TARunInfo* runinfo)
   {
      printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      fRunEventCounter = 0;
   }

   bool Check(size_t itrgchan, size_t ichan)
   {
      bool all_ok = true;
      size_t all_count = 0;
      printf("Check cbtrg chan %zu and chronobox channel %zu:\n", itrgchan, ichan);
      assert(fCbBanks.size() == 5);
      std::vector<size_t> ihit;
      ihit.resize(fCbBanks.size());
      size_t itrgbank = 0;
      printf("%s: chan %zu: %zu hits\n", fCbBanks[itrgbank].c_str(), itrgchan, fCbHits[itrgbank][itrgchan].size());
      for (size_t ibank=1; ibank<=4; ibank++) {
         printf("%s: chan %zu: %zu hits\n", fCbBanks[ibank].c_str(), ichan, fCbHits[ibank][ichan].size());
      }
      while (1) {
         std::vector<double> tv;
         if (1) {
            if (ihit[itrgbank] >= fCbHits[itrgbank][itrgchan].size())
               break;
            while (fCbHits[itrgbank][itrgchan][ihit[itrgbank]].flags&CB_HIT_FLAG_TE) {
               if (ihit[itrgbank] >= fCbHits[itrgbank][itrgchan].size())
                  break;
               ihit[itrgbank]++;
            }
            if (ihit[itrgbank] >= fCbHits[itrgbank][itrgchan].size())
               break;
            tv.push_back(fCbHits[itrgbank][itrgchan][ihit[itrgbank]].time);
            ihit[itrgbank]++;
         }
         for (size_t ibank=1; ibank<=4; ibank++) {
            if (ihit[ibank] >= fCbHits[ibank][ichan].size())
               break;
            while (fCbHits[ibank][ichan][ihit[ibank]].flags&CB_HIT_FLAG_TE) {
               if (ihit[ibank] >= fCbHits[ibank][ichan].size())
                  break;
               ihit[ibank]++;
            }
            if (ihit[ibank] >= fCbHits[ibank][ichan].size())
               break;
            tv.push_back(fCbHits[ibank][ichan][ihit[ibank]].time);
            ihit[ibank]++;
         }
         if (tv.size() != 5)
            break;

         bool ok = true;
         for (size_t i=0; i<5; i++) {
            if (fabs(tv[i]-tv[0]) > 0.000009)
               ok = false;
         }
         all_ok &= ok;
         all_count++;

         if (!all_ok) {
            for (size_t i=0; i<5; i++) {
               printf(" %.6f", tv[i]);
            }
            printf(" ok %d\n", ok);
         }
      }

      printf("Check cbtrg chan %zu and chronobox channel %zu, %zu hits, ok %d\n", itrgchan, ichan, all_count, all_ok);
      return all_ok;
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("EndRun, run %d\n", runinfo->fRunNo);
      printf("Counted %d events in run %d\n", fRunEventCounter, runinfo->fRunNo);
      double elapsed = fTE - fT0;
      double minutes = elapsed/60.0;
      double hours = minutes/60.0;
      printf("Elapsed time %d -> %d is %.0f sec or %f minutes or %f hours\n", fT0, fTE, elapsed, minutes, hours);
      if (fCountMissingWC) {
         printf("Wraparounds: %d, missing %d, cannot compute TS frequency\n", fCountWC, fCountMissingWC);
      } else {
         double tsbits = (1<<24);
         printf("Wraparounds: %d, approx rate %f Hz, ts freq %.1f Hz\n", fCountWC, fCountWC/elapsed, 0.5*fCountWC/elapsed*tsbits);
      }

      for (size_t ibank=0; ibank<fCbBanks.size(); ibank++) {
         for (size_t ichan=0; ichan<fCbHits[ibank].size(); ichan++) {
            if (fCbHits[ibank][ichan].size() > 0) {
               printf("%s: chan %zu: %zu hits\n", fCbBanks[ibank].c_str(), ichan, fCbHits[ibank][ichan].size());
            }
         }
      }

      for (size_t ibank=0; ibank<fCbBanks.size(); ibank++) {
         for (size_t ichan=0; ichan<fCbHits[ibank].size(); ichan++) {
            if (ibank==0 && ichan!=0 && ichan!=3 && ichan!=4)
               continue;
            if (ibank>0 && ichan != 33)
               continue;
            printf("%s: chan %zu: %zu hits\n", fCbBanks[ibank].c_str(), ichan, fCbHits[ibank][ichan].size());
            //if (ibank==0)
            //   continue;
            for (size_t ihit=0; ihit<fCbHits[ibank][ichan].size(); ihit++) {
               const CbHit* phit = &fCbHits[ibank][ichan][ihit];
               if (phit->flags & CB_HIT_FLAG_TE)
                  continue;
               printf("%s: hit %zu, time %.6f sec, %d+%d, channel %2d (%d)\n", fCbBanks[ibank].c_str(), ihit, phit->time, phit->timestamp, phit->epoch, phit->channel, (phit->flags&CB_HIT_FLAG_TE));
            }
         }
      }

      Check(0, 33);
      Check(1, 33);
      Check(2, 33);
      Check(3, 33);
      Check(3, 36);
   }

   void AnalyzeCbFifo(int ibank, CbUnpack* cb, TMEvent* event, TMBank* cbbank)
   {
      int nwords = cbbank->data_size/4; // byte to uint32_t words
      uint32_t* cbdata = (uint32_t*)event->GetBankData(cbbank);

      CbHits hits;
      CbScalers scalers;

      cb->Unpack(cbdata, nwords, &hits, &scalers);

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

      for (size_t i=0; i<hits.size(); i++) {
         fCbHits[ibank][hits[i].channel].push_back(hits[i]);
      }

      //if (scalers.size() > 0) {
      //   if (print) {
      //      printf("Data from %s: ", name);
      //      PrintCbScalers(scalers);
      //   }
      //}
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      return flow;

      //      printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
      //event->FindAllBanks();
      //printf("Event: %s\n", event->HeaderToString().c_str());
      //printf("Banks: %s\n", event->BankListToString().c_str());

      //if (event->serial_number == 0) {
      //   printf("AAA\n");
      //   //fUnpackCb01.fEpochFromReset = true;
      //   exit(1);
      //}

      fRunEventCounter++;

      event->FindAllBanks();

      for (size_t ibank=0; ibank<fCbBanks.size(); ibank++) {
         TMBank* cbbank = event->FindBank(fCbBanks[ibank].c_str());
         if (cbbank)
            AnalyzeCbFifo(ibank, fCbUnpack[ibank], event, cbbank);
      }

      //cbbank = event->FindBank("CBF2");
      //if (cbbank)
      //   AnalyzeCbFifo("cb02", 2, &fUnpackCb02, event, cbbank);

      //cbbank = event->FindBank("CBF3");
      //if (cbbank)
      //   AnalyzeCbFifo("cb03", 3, &fUnpackCb03, event, cbbank);

      //cbbank = event->FindBank("CBF4");
      //if (cbbank)
      //   AnalyzeCbFifo("cb04", 4, &fUnpackCb04, event, cbbank);

      // TRG chronobox data
      //cbbank = event->FindBank("CBFT");
      //if (cbbank)
      //   AnalyzeCbFifo("cbtrg", 0, &fUnpackCbTRG, event, cbbank);

      return flow;
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      return flow;
   }
};

class CbkoModuleFactory: public TAFactory
{
public:
   void Usage()
   {
      //      printf("\tCbModuleFactory Usage:\n");
      //printf("\t--print-cb-data : print chronobox raw data\n");
      //#ifdef HAVE_ROOT
      //printf("\t--dumpchronojson write out the chronobox channel list as json\n");
      //printf("\t--loadchronojson filename.json override odb channel list with json\n");
      //#endif
   }

   void Init(const std::vector<std::string> &args)
   {
      //printf("Init!\n");
      //printf("Arguments:\n");
      //for (unsigned i=0; i<args.size(); i++) {
      //   printf("arg[%d]: [%s]\n", i, args[i].c_str());
      //   if (args[i] == "--print-cb-data") {
      //      fFlags.fPrint = true;
      //   }
      //#ifdef HAVE_ROOT
      //   if (args[i] == "--dumpchronojson")
      //      fFlags.fDumpJsonChannelNames = true;
      //   if (args[i] == "--loadchronojson")
      //   {
      //      fFlags.fLoadJsonChannelNames = true;
      //      i++;
      //      fFlags.fLoadJsonFile=args[i];
      //   }
      //#endif
      //}
   }
   
   void Finish()
   {
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      return new CbkoModule(runinfo);
   }
};

static TARegister tar(new CbkoModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
