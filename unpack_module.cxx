//
// unpack_module.cxx
//
// K.Olchanski
//

#include <stdio.h>

#include "manalyzer.h"
#include "midasio.h"

#include "Unpack.h"
#include "FeamEVB.h"
#include "AgEVB.h"
#include "AgFlow.h"

#include "ncfm.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

static std::string join(const char* sep, const std::vector<std::string> &v)
{
   std::string s;
   for (unsigned i=0; i<v.size(); i++) {
      if (i>0)
         s += sep;
      s += v[i];
   }
   return s;
}

static std::vector<std::string> split(const std::string& s, char seperator)
{
   std::vector<std::string> output;
   
   std::string::size_type prev_pos = 0, pos = 0;

   while((pos = s.find(seperator, pos)) != std::string::npos)
      {
         std::string substring( s.substr(prev_pos, pos-prev_pos) );
         output.push_back(substring);
         prev_pos = ++pos;
      }

   output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word
   return output;
}

class UnpackModule: public TARunObject
{
public:
   Ncfm*       fCfm = NULL;
   Alpha16EVB* fA16Evb = NULL;
   FeamEVB*    fFeamEvb = NULL;
   AgEVB*      fAgEvb = NULL;

   std::vector<std::string> fFeamBanks;
   std::vector<std::string> fAdcMap;
   
   UnpackModule(TARunInfo* runinfo)
      : TARunObject(runinfo)
   {
      printf("UnpackModule::ctor!\n");

      fCfm     = new Ncfm(getenv("AG_CFM"));
      fA16Evb  = NULL;
      fFeamEvb = NULL;
      fAgEvb   = NULL;
   }

   ~UnpackModule()
   {
      printf("UnpackModule::dtor!\n");
      DELETE(fCfm);
      DELETE(fA16Evb);
      DELETE(fFeamEvb);
      DELETE(fAgEvb);
   }

   void LoadFeamBanks(int runno)
   {
      fFeamBanks = fCfm->ReadFile("feam", "banks", runno);
      printf("Loaded feam banks: %s\n", join(" ", fFeamBanks).c_str());
   }

   void LoadAdcMap(int runno)
   {
      fAdcMap = fCfm->ReadFile("adc", "map", runno);
      printf("Loaded adc map: %s\n", join(", ", fAdcMap).c_str());
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("UnpackRun::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      LoadFeamBanks(runinfo->fRunNo);

      bool have_a16  = true;
      bool have_feam = true;

      if (runinfo->fRunNo == 537) {
         have_a16 = false;
      }

      if (have_a16) {
         fA16Evb  = new Alpha16EVB();
         fA16Evb->Reset();
         fA16Evb->Configure(runinfo->fRunNo);
         LoadAdcMap(runinfo->fRunNo);
         fA16Evb->fMap.Init(fAdcMap);
         fA16Evb->fMap.Print();
      }

      if (have_feam) {
         fFeamEvb = new FeamEVB(fFeamBanks.size(), 125.0*1e6, 100000/1e9);
         //fFeamEvb->fSync.fTrace = true;
      }
      
      fAgEvb = new AgEVB(100.0*1e6, 125.0*1e6, 50.0*1e-6, 100, 90, true);
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("UnpackRun::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));

      int count_feam = 0;

      if (fFeamEvb) {
         //printf("UnpackRun::EndRun: FeamEVB state:\n");
         //fFeamEvb->Print();

         while (1) {
            FeamEvent *e = fFeamEvb->GetLastEvent();
            if (!e)
               break;
            
            if (1) {
               printf("Unpacked FEAM event: ");
               e->Print();
               printf("\n");
            }
            
            if (0) {
               for (unsigned i=0; i<e->modules.size(); i++) {
                  if (!e->modules[i])
                     break;
                  FeamModuleData* m = e->modules[i];
                  printf("position %2d, cnt %4d, ts_trig: 0x%08x %14.3f usec, ts_incr %14.3f usec\n", m->fPosition, m->cnt, m->ts_trig, m->fTime*1e6, m->fTimeIncr*1e6);
                  //a->Print();
                  //printf("\n");
               }
            }
            
            if (fAgEvb) {
               count_feam += 1;
               fAgEvb->AddFeamEvent(e);
               e = NULL;
            }
            
            DELETE(e);
         }

         printf("UnpackRun::EndRun: FeamEVB final state:\n");
         fFeamEvb->Print();
      }
      
      printf("UnpackRun::EndRun: Unpacked %d last FEAM events\n", count_feam);

      // Handle leftover AgEVB events

      int count_agevent = 0;

      if (fAgEvb) {
         printf("UnpackRun::EndRun: AgEVB state:\n");
         fAgEvb->Print();
         
         while (1) {
            AgEvent *e = fAgEvb->GetLastEvent();
            if (!e)
               break;
            
            if (1) {
               printf("Unpacked AgEvent: ");
               e->Print();
               printf("\n");
            }
            
            count_agevent += 1;
            
            delete e;
         }

         printf("UnpackRun::EndRun: AgEVB final state:\n");
         fAgEvb->Print();
      }
      
      printf("UnpackRun::EndRun: Unpacked %d last AgEvent events\n", count_agevent);
   }
   
   void PauseRun(TARunInfo* runinfo)
   {
      printf("UnpackRun::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id != 1)
         return flow;

      if (fA16Evb) {
         Alpha16Event* e = UnpackAlpha16Event(fA16Evb, event);

         if (e) {
            if (1) {
               printf("Unpacked Alpha16 event: ");
               e->Print();
               printf("\n");
            }

            if (fAgEvb && e->time >= 0) {
               fAgEvb->AddAlpha16Event(e);
               e = NULL;
            }

            DELETE(e);
         }
      }

      if (fFeamEvb) {
         FeamEvent *e = UnpackFeamEvent(fFeamEvb, event, fFeamBanks);
         if (e) {
            if (1) {
               printf("Unpacked FEAM event: ");
               e->Print();
               printf("\n");
            }

            if (0) {
               for (unsigned i=0; i<e->modules.size(); i++) {
                  if (!e->modules[i])
                     break;
                  FeamModuleData* m = e->modules[i];
                  printf("position %2d, cnt %4d, ts_trig: 0x%08x %14.3f usec, ts_incr %14.3f usec\n", m->fPosition, m->cnt, m->ts_trig, m->fTime*1e6, m->fTimeIncr*1e6);
                  //a->Print();
                  //printf("\n");
               }
            }

            if (fAgEvb) {
               fAgEvb->AddFeamEvent(e);
               e = NULL;
            }

            DELETE(e);
         }
      }

      if (fAgEvb) {
         AgEvent* e = fAgEvb->Get();

         if (e) {
            printf("Have AgEvent: ");
            e->Print();
            printf("\n");

            if (e->complete && e->a16 && e->feam) {
               double ta = e->a16->timeIncr;
               double tf = e->feam->timeIncr;
               printf("  incr %f %f sec, diff %f ns, count %d\n", ta, tf, (tf-ta)*1e9, e->counter);
            }

            flow = new AgEventFlow(flow, e);
         }
      }

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("UnpackRun::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class UnpackModuleFactory: public TAFactory
{
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("UnpackModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
#if 0
         if (args[i] == "--nopads")
            fDoPads = false;
         if (args[i] == "--plot1")
            fPlotPad = atoi(args[i+1].c_str());
#endif
      }
   }

   void Finish()
   {
      printf("UnpackModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("UnpackModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new UnpackModule(runinfo);
   }
};

static TARegister tar(new UnpackModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
