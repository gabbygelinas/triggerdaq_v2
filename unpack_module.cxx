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

#if 0
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
#endif

class UnpackFlags
{
public:
   bool fPrint = false;
   bool fNoAdc = false;
   bool fNoPwb = false;
};

class UnpackModule: public TARunObject
{
public:
   UnpackFlags* fFlags = NULL;
   Ncfm*       fCfm = NULL;
   Alpha16Asm* fAdcAsm = NULL;
   FeamEVB*    fFeamEvb = NULL;
   AgEVB*      fAgEvb = NULL;

   std::vector<std::string> fFeamBanks;
   std::vector<std::string> fAdcMap;
   
   UnpackModule(TARunInfo* runinfo, UnpackFlags* flags)
      : TARunObject(runinfo)
   {
      printf("UnpackModule::ctor!\n");

      fFlags   = flags;
      fCfm     = new Ncfm(getenv("AG_CFM"));
      fAdcAsm  = NULL;
      fFeamEvb = NULL;
      fAgEvb   = NULL;
   }

   ~UnpackModule()
   {
      printf("UnpackModule::dtor!\n");
      DELETE(fCfm);
      DELETE(fAdcAsm);
      DELETE(fFeamEvb);
      DELETE(fAgEvb);
   }

   bool LoadFeamBanks(int runno)
   {
      fFeamBanks = fCfm->ReadFile("feam", "banks", runno);
      printf("Loaded feam banks: %s\n", join(" ", fFeamBanks).c_str());
      return fFeamBanks.size() > 0;
   }

   bool LoadAdcMap(int runno)
   {
      fAdcMap = fCfm->ReadFile("adc", "map", runno);
      printf("Loaded adc map: %s\n", join(", ", fAdcMap).c_str());
      return fAdcMap.size() > 0;
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("UnpackRun::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      bool have_trg = true;
      bool have_adc = true;
      bool have_pwb = true;

      if (runinfo->fRunNo < 1244)
         have_trg = false;

      if (fFlags->fNoAdc)
         have_adc = false;

      if (fFlags->fNoPwb)
         have_pwb = false;

      if (have_adc)
         have_adc &= LoadAdcMap(runinfo->fRunNo);

      if (have_pwb)
         have_pwb &= LoadFeamBanks(runinfo->fRunNo);

      if (have_adc) {
         fAdcAsm  = new Alpha16Asm();
         fAdcAsm->fMap.Init(fAdcMap);
         fAdcAsm->fMap.Print();
         if (runinfo->fRunNo < 808) {
            fAdcAsm->fTsFreq = 100e6;
         }
      }

      if (have_pwb) {
         fFeamEvb = new FeamEVB(fFeamBanks.size(), 125.0*1e6, 100000/1e9);
         //fFeamEvb->fSync.fTrace = true;
      }

      fAgEvb = new AgEVB(62.5*1e6, 125.0*1e6, 125.0*1e6, 50.0*1e-6, 100, 90, true);
      //fAgEvb->fSync.fTrace = true;

      if (!have_trg)
         fAgEvb->fSync.fModules[AGEVB_TRG_SLOT].fDead = true;
      if (!have_adc)
         fAgEvb->fSync.fModules[AGEVB_ADC_SLOT].fDead = true;
      if (!have_pwb)
         fAgEvb->fSync.fModules[AGEVB_PWB_SLOT].fDead = true;
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
               printf("Unpacked PWB event: ");
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

      if (1) {
         const time_t now = time(NULL);
         const time_t t = event->time_stamp;
         int dt = now - t;
         printf("UnpackModule: serial %d, time %d, age %d, date %s\n", event->serial_number, event->time_stamp, dt, ctime(&t));
      }

      if (1) {
         TMBank* atat_bank = event->FindBank("ATAT");

         if (atat_bank) {
            TrigEvent* e = UnpackTrigEvent(event, atat_bank);

            if (e) {
               if (1) {
                  printf("Unpacked TRG event: ");
                  e->Print();
                  printf("\n");
               }
               
               if (fAgEvb && e->time >= 0) {
                  fAgEvb->AddTrigEvent(e);
                  e = NULL;
               }
               
               DELETE(e);
            }
         }
      }

      if (fAdcAsm) {
         Alpha16Event* e = UnpackAlpha16Event(fAdcAsm, event);

         if (e) {
            if (1) {
               printf("Unpacked ADC event: ");
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
         while (1) {
            if (!e)
               break;

            if (1) {
               printf("Unpacked PWB event: ");
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

            e = fFeamEvb->Get();
         }
      }

      if (fAgEvb) {
         while (1) {
            AgEvent* e = fAgEvb->Get();
            if (!e)
               break;

            if (1) {
               printf("Unpacked AgEvent:   ");
               e->Print();
               printf("\n");
            }

            runinfo->fFlowQueue.push_back(new AgEventFlow(NULL, e));
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
   UnpackFlags fFlags;

public:
   void Init(const std::vector<std::string> &args)
   {
      printf("UnpackModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--print")
            fFlags.fPrint = true;
         if (args[i] == "--noadc")
            fFlags.fNoAdc = true;
         if (args[i] == "--nopwb")
            fFlags.fNoPwb = true;
      }
   }

   void Finish()
   {
      printf("UnpackModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("UnpackModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new UnpackModule(runinfo, &fFlags);
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
