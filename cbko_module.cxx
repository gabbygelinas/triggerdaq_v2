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

#include "AgFlow.h"

class Flags
{
public:
   bool fPrint = false;     // print something for every event
   bool fPrintHits = false; // print all hits
   bool fCheck = false;     // enable consistency checks
};

class CbkoModule: public TARunObject
{
public:
   Flags* fFlags = NULL;

   std::vector<std::string> fCbBanks;
   std::vector<std::vector<std::vector<CbHit>>> fCbHits;

   bool fCheckHitsOk = true;

   CbkoModule(TARunInfo* runinfo, Flags* flags) : TARunObject(runinfo)
   {
      fModuleName="cbko_module";
      fFlags = flags;
   }

   ~CbkoModule()
   {
      // empty
   }
  
   void BeginRun(TARunInfo* runinfo)
   {
      //printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
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
      if (fCbBanks.size() != 5) {
         fprintf(stderr, "CbkoModule::Check: Error: Unexpected number of chronobox banks: %zu, expected %d\n", fCbBanks.size(), 5);
         return false;
      }
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
      if (fCbBanks.size() != 5) {
         fprintf(stderr, "CbkoModule::Check: Error: Unexpected number of chronobox banks: %zu, expected %d\n", fCbBanks.size(), 5);
         return false;
      }
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

         if (0 && !all_ok) {
            for (size_t i=0; i<4; i++) {
               printf(" %.6f", tv[i]);
            }
            printf(" ok %d\n", ok);
         }
      }

      printf("Check chronobox channel %zu across cb01..04: %zu hits, ok %d\n", ichan, all_count, all_ok);
      return all_ok;
   }

   bool Check2(size_t ibank1, size_t ichan1, size_t ibank2, size_t ichan2)
   {
      if (ibank1>=fCbBanks.size() || ibank2>=fCbBanks.size()) {
         fprintf(stderr, "CbkoModule::Check2: No chronobox data for banks %zu and %zu\n", ibank1, ibank2);
         return false;
      }

      if (ichan1>=fCbHits[ibank1].size()) {
         fprintf(stderr, "CbkoModule::Check2: No chronobox hits in bank %zu channel %zu\n", ibank1, ichan1);
         return false;
      }

      if (ichan2>=fCbHits[ibank2].size()) {
         fprintf(stderr, "CbkoModule::Check2: No chronobox hits in bank %zu channel %zu\n", ibank2, ichan2);
         return false;
      }

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
      //double t2prev = -1;
      
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
               printf("hit %zu %zu time %.6f %.6f match (diff %.6f, timestamp drift %.6f)\n", i1, i2, t1, t2, t1-t2, drift);
            drift = t1-t2;
            all_count++;
            i1++;
            i2++;
            t1prev = t1;
            //t2prev = t2;
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

      printf("Check bank %s channel %2zu against bank %s channel %2zu: matching %zu, bank %s: missing %zu, dupe %zu, missing in bank %s: %zu, timestamp drift %.6f. ok %d\n", fCbBanks[ibank1].c_str(), ichan1, fCbBanks[ibank2].c_str(), ichan2, all_count, fCbBanks[ibank1].c_str(), missing1, dupe1, fCbBanks[ibank2].c_str(), missing2, drift, all_ok);

      return all_ok;
   }

   void EndRun(TARunInfo* runinfo)
   {
      bool ok = true;

      if (!fCheckHitsOk) {
         ok = false;
         printf("CbkoModule::EndRun: Chronobox data self-consistency check: CheckHits() FAILED\n");
      } else {
         printf("CbkoModule::EndRun: Chronobox data self-consistency check: CheckHits() ok!\n");
      }

      if (fFlags->fCheck) {
         printf("CbkoModule::EndRun: List of chronobox channels with hits:\n");
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

         printf("CbkoModule::EndRun: Chronobox cb01..04 synchronization cross check:\n");

         bool okc = true;
         okc &= Check4(33);
         okc &= Check4(36);

         if (!okc)
            printf("CbkoModule::EndRun: Chronobox cb01..04 synchronization cross check: FAILED!\n");

         ok &= okc;
         
         //Check(0, 33);
         //Check(1, 33);
         //Check(2, 33);
         //Check(3, 33);
         //Check(3, 36);

         printf("CbkoModule::EndRun: Chronobox check:\n");

         ok &= Check2(1, 33, 1, 36);
         ok &= Check2(1, 33, 2, 33);
         ok &= Check2(0, 0, 1, 33);
         ok &= Check2(0, 1, 1, 33);
         ok &= Check2(0, 2, 1, 33);
         ok &= Check2(0, 3, 1, 33);
         
         if (ok)
            printf("CbkoModule::EndRun: Chronobox checks ok!\n");
         else
            printf("CbkoModule::EndRun: Chronobox checks FAILED!\n");
      }
   }

   bool CheckHitsRange(const char* bank_name, const CbHits& hits, size_t first, size_t last)
   {
      bool ok = true;

      double min_time = hits[first].time;
      double max_time = hits[first].time;
      int min_epoch = hits[first].epoch;
      int max_epoch = hits[first].epoch;
      
      for (size_t i=first; i<last; i++) {
         double time = hits[i].time;
         int epoch = hits[i].epoch;
         if (time < min_time) min_time = time;
         if (time > max_time) max_time = time;
         if (epoch < min_epoch) min_epoch = epoch;
         if (epoch > max_epoch) max_epoch = epoch;
      }
      
      if (max_time - min_time > 1.10) {
         ok = false;
      }

      if (!ok) {
         //printf("range %zu..%zu: ", first, last);
         printf("CbkoModule::CheckHitsRange: %s, hits: %6zu, time: %.6f..%.6f (diff %.6f) sec, epoch: %d..%d (diff %d), inconsistent time or epoch!\n", bank_name, last-first, min_time, max_time, max_time-min_time, min_epoch, max_epoch, max_epoch-min_epoch);
      }
      
      if (!ok) {
         for (size_t i=first; i<last; i++) {
            //int channel = hits[i].channel;
            //double time = hits[i].time;
            //int epoch = hits[i].epoch;
            
            //printf("%s: hit %zu, time %.6f sec, %d+%d, channel %2d (%d)\n", bank_name, i, hits[i].time, hits[i].timestamp, hits[i].epoch, hits[i].channel, (hits[i].flags&CB_HIT_FLAG_TE));
         }
      }

      return ok;
   }

   bool CheckHits(const char* bank_name, const CbHits& hits)
   {
      if (hits.empty())
         return true;

      bool ok = true;

      if (hits.size() < 10000) {
         ok &= CheckHitsRange(bank_name, hits, 0, hits.size());
      } else {
         size_t first = 0;
         while (1) {
            bool done = false;
            size_t last = first + 10000;
            if (last > hits.size()) {
               last = hits.size();
               done = true;
            }
            ok &= CheckHitsRange(bank_name, hits, first, last);
            if (done)
               break;
            first = last;
         }
      }

      return ok;
   }

   void AnalyzeCbHits(CbHitsFlow* hf)
   {
      int ibank = hf->fCbIndex;

      if (fFlags->fPrint || fFlags->fPrintHits) {
         printf("%s: index %d, num_inputs %d, hits: %zu\n", hf->fCbBankName.c_str(), hf->fCbIndex, hf->fNumInputs, hf->fHits.size());
         if (fFlags->fPrintHits)
            PrintCbHits(hf->fHits);
      }

      //if (ibank == 0) {
      //for (size_t i=0; i<hits.size(); i++) {
      //   printf("%s: hit %zu, time %.6f sec, %d+%d, channel %2d (%d)\n", cbbank->name.c_str(), i, hits[i].time, hits[i].timestamp, hits[i].epoch, hits[i].channel, (hits[i].flags&CB_HIT_FLAG_TE));
      //}
      //}

      fCheckHitsOk &= CheckHits(hf->fCbBankName.c_str(), hf->fHits);

      if (fFlags->fCheck) {
         if (size_t(ibank) >= fCbHits.size())
            fCbHits.resize(ibank+1);
         if (size_t(ibank) >= fCbBanks.size())
            fCbBanks.resize(ibank+1);
         fCbBanks[ibank] = hf->fCbBankName;
         for (size_t i=0; i<hf->fHits.size(); i++) {
            CbHit* hit = &hf->fHits[i];
            if (hit->channel >= fCbHits[ibank].size())
               fCbHits[ibank].resize(hit->channel+1);
            fCbHits[ibank][hit->channel].push_back(*hit);
         }
      }
   }

   //TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   //{
   //   return flow;
   //}

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      CbHitsFlow* hf = flow->Find<CbHitsFlow>();
      if (hf)
         AnalyzeCbHits(hf);
      return flow;
   }
};

class CbkoModuleFactory: public TAFactory
{
public:
   Flags fFlags;

   void Usage()
   {
      printf("CbkoModuleFactory Usage:\n");
      printf("--print-cb-data : print chronobox data\n");
      printf("--check-cb-data : check chronobox data\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      //printf("Init!\n");
      //printf("Arguments:\n");
      for (unsigned i=0; i<args.size(); i++) {
         //printf("arg[%d]: [%s]\n", i, args[i].c_str());
         if (args[i] == "--print-cb-data") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--print-cb-hits") {
            fFlags.fPrintHits = true;
         }
         if (args[i] == "--check-cb-data") {
            fFlags.fCheck = true;
         }
      }
   }
   
   void Finish()
   {
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      return new CbkoModule(runinfo, &fFlags);
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
