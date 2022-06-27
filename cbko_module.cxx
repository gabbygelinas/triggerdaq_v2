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

   bool KillDupes(size_t ibank, size_t ichan)
   {
      bool all_ok = true;
      size_t nhits = fCbHits[ibank][ichan].size();
      if (nhits < 1)
         return all_ok;
      //printf("%s: chan %zu: %zu hits\n", fCbBanks[ibank].c_str(), ichan, fCbHits[ibank][ichan].size());
      double time_prev = 0;
      double min_dt = 99999;
      size_t count = 0;
      for (size_t ihit=0; ihit<nhits-1; ihit++) {
         if (fCbHits[ibank][ichan][ihit].flags&CB_HIT_FLAG_TE)
            continue;
         double time = fCbHits[ibank][ichan][ihit].time;
         double dt = time-time_prev;
         //if (ibank==0 && ichan==0)
            //   printf("dt %f\n", dt);
         time_prev = time;
         if (dt < 0) {
            printf("time goes backwards!\n");
            all_ok = false;
            break;
         }
         
         if (dt < 0.000001) {
            //printf("%s: chan %zu: hit %zu: duplicate time %.6f\n", fCbBanks[ibank].c_str(), ichan, ihit, time);
            fCbHits[ibank][ichan][ihit].flags = -1;
            count++;
            continue;
         }
         
         if (dt < min_dt)
            min_dt = dt;
      }
      printf("%s: chan %zu: %zu hits, min_dt %.6f, %zu duplicates\n", fCbBanks[ibank].c_str(), ichan, fCbHits[ibank][ichan].size(), min_dt, count);
      return all_ok;
   }

   bool KillDupes()
   {
      bool all_ok = true;
      printf("Check duplicates:\n");
      for (size_t ibank=0; ibank<fCbBanks.size(); ibank++) {
         size_t nchan = fCbHits[ibank].size();
         for (size_t ichan=0; ichan<nchan; ichan++) {
            all_ok |= KillDupes(ibank, ichan);
         }
      }
      printf("Check duplicates: ok %d\n", all_ok);
      return all_ok;
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

   bool Check4(size_t ichan)
   {
      bool all_ok = true;
      size_t all_count = 0;
      printf("Check chronobox channel %zu:\n", ichan);
      assert(fCbBanks.size() == 5);
      std::vector<size_t> ihit;
      ihit.resize(fCbBanks.size());
      for (size_t ibank=1; ibank<=4; ibank++) {
         printf("%s: chan %zu: %zu hits\n", fCbBanks[ibank].c_str(), ichan, fCbHits[ibank][ichan].size());
      }
      while (1) {
         std::vector<double> tv;
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
         if (tv.size() != 4)
            break;

         bool ok = true;
         for (size_t i=0; i<4; i++) {
            if (fabs(tv[i]-tv[0]) > 0.000001)
               ok = false;
         }
         all_ok &= ok;
         all_count++;

         if (!all_ok) {
            for (size_t i=0; i<4; i++) {
               printf(" %.6f", tv[i]);
            }
            printf(" ok %d\n", ok);
         }
      }

      printf("Check chronobox channel %zu, %zu hits, ok %d\n", ichan, all_count, all_ok);
      return all_ok;
   }

   bool Check2(size_t ibank1, size_t ichan1, size_t ibank2, size_t ichan2)
   {
      //bool print = true;
      bool print = false;
      //if (ibank1 == 0) print = true;
      bool all_ok = true;
      size_t all_count = 0;
      if (print) {
         printf("Check bank %zu channel %zu against bank %zu channel %zu:\n", ibank1, ichan1, ibank2, ichan2);
         printf("%s: chan %zu: %zu hits, first time %.6f\n", fCbBanks[ibank1].c_str(), ichan1, fCbHits[ibank1][ichan1].size(), fCbHits[ibank1][ichan1][0].time);
         printf("%s: chan %zu: %zu hits, first time %.6f\n", fCbBanks[ibank2].c_str(), ichan2, fCbHits[ibank2][ichan2].size(), fCbHits[ibank2][ichan2][0].time);
      }
      size_t n1 = fCbHits[ibank1][ichan1].size();
      size_t n2 = fCbHits[ibank2][ichan2].size();

      size_t i1 = 0;
      size_t i2 = 0;

      double drift = 0;

      size_t missing1 = 0;
      size_t missing2 = 0;

      size_t dupe1 = 0;

      double t1prev = -1;
      double t2prev = -1;
      
      while (1) {
         if (i1 >= n1)
            break;
         if (i2 >= n2)
            break;

         while (fCbHits[ibank1][ichan1][i1].flags&CB_HIT_FLAG_TE) {
            i1++;
            if (i1 >= n1)
               break;
         }

         while (fCbHits[ibank2][ichan2][i2].flags&CB_HIT_FLAG_TE) {
            i2++;
            if (i2 >= n2)
               break;
         }

         if (i1 >= n1)
            break;
         if (i2 >= n2)
            break;

         double t1 = fCbHits[ibank1][ichan1][i1].time;
         double t2 = fCbHits[ibank2][ichan2][i2].time;

         if (fabs(t1 - t1prev) < 0.000001) {
            if (print)
               printf("hit %zu %zu time %.6f dupe\n", i1, i2, t1);
            dupe1++;
            all_ok = false;
            i1++;
            t1prev = t1;
            continue;
         }

         if (fabs(t1-t2-drift) < 0.000002) {
            if (print)
               printf("hit %zu %zu time %.6f %.6f match (diff %.6f, drift %.6f)\n", i1, i2, t1, t2, t1-t2, drift);
            drift = t1-t2;
            all_count++;
            i1++;
            i2++;
            t1prev = t1;
            t2prev = t2;
         } else if (t1 > t2) {
            if (print)
               printf("hit %zu %zu time %.6f %.6f mismatch\n", i1, i2, t1, t2);
            // bank1 ahead of bank2, missing a hit? check next hits in bank2
            missing1++;
            all_ok = false;
            i2++;
            all_count++;
         } else if (t1 < t2) {
            if (print)
               printf("hit %zu %zu time %.6f %.6f mismatch\n", i1, i2, t1, t2);
            // bank2 ahead of bank1, missing a hit? check next hits in bank1
            missing2++;
            all_ok = false;
            i1++;
            all_count++;
         } else {
            printf("hit %zu %zu time %.6f %.6f mismatch, stopping\n", i1, i2, t1, t2);
            all_ok = false;
            break;
         }
      }

      if (fabs(drift) > 0.000001)
         all_ok = false;

      printf("Check bank %s channel %2zu against bank %s channel %2zu: matching %zu, bank %s: missing %zu, dupe %zu, missing in bank %s: %zu, drift %.6f. ok %d\n", fCbBanks[ibank1].c_str(), ichan1, fCbBanks[ibank2].c_str(), ichan2, all_count, fCbBanks[ibank1].c_str(), missing1, dupe1, fCbBanks[ibank2].c_str(), missing2, drift, all_ok);

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

      if (0) {
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
      }

      //KillDupes();

      //KillDupes(0, 0);
      //KillDupes(0, 1);
      //KillDupes(0, 2);
      //KillDupes(0, 3);

      Check4(33);
      Check4(36);

      //Check(0, 33);
      //Check(1, 33);
      //Check(2, 33);
      //Check(3, 33);
      //Check(3, 36);

      Check2(1, 33, 1, 36);
      Check2(1, 33, 2, 33);
      Check2(0, 0, 1, 33);
      Check2(0, 1, 1, 33);
      Check2(0, 2, 1, 33);
      Check2(0, 3, 1, 33);
   }

   void AnalyzeCbFifo(int ibank, CbUnpack* cb, TMEvent* event, TMBank* cbbank)
   {
      int nwords = cbbank->data_size/4; // byte to uint32_t words
      uint32_t* cbdata = (uint32_t*)event->GetBankData(cbbank);

      CbHits hits;
      CbScalers scalers;

      if (ibank==0)
         cb->fKludge = 1;

      //if (ibank==0)
      //   cb->fVerbose = true;
      //else
      //   cb->fVerbose = false;

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
