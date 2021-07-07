//
// unpack_module.cxx
//
// K.Olchanski
//

#include <stdio.h>
#include <assert.h>

#include "manalyzer.h"
#include "midasio.h"

#include "AgAsm.h"
#include "AgFlow.h"

#include "ncfm.h"

#ifdef _TIME_ANALYSIS_
#include "AnalysisTimer.h"
#endif

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

#include "PwbAsm.h"

PwbAsm* gPwbAsm = NULL;

class UnpackFlags
{
public:
   bool fPrint = false;
   int  fAgeSec = 0;
   int  fSleepUSec = 0;
};

class UnpackModule: public TARunObject
{
public:
   UnpackFlags* fFlags = NULL;
   Ncfm*       fCfm = NULL;
   AgAsm*      fAgAsm = NULL;

   std::vector<std::string> fFeamBanks;
   std::vector<std::string> fAdcMap;

   bool fTrace = false;
   
   UnpackModule(TARunInfo* runinfo, UnpackFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("UnpackModule::ctor!\n");

      fFlags   = flags;
      std::string agcfmdb_path="agcfmdb";
      if (getenv("AGRELEASE"))
      {
         agcfmdb_path=getenv("AGRELEASE");
         agcfmdb_path+="/agana/agcfmdb";
      }
      fCfm     = new Ncfm(agcfmdb_path.c_str());
   }

   ~UnpackModule()
   {
      if (fTrace)
         printf("UnpackModule::dtor!\n");
      if (fCfm) {
         delete fCfm;
         fCfm = NULL;
      }
      if (fAgAsm) {
         delete fAgAsm;
         fAgAsm = NULL;
      }
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
      if (fTrace)
         printf("UnpackModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      int adc32_rev = 0;
      if (0) {
      } else if (runinfo->fRunNo >= 902660) {
         adc32_rev = 11;
      } else if (runinfo->fRunNo >= 900000) {
         adc32_rev = 1;
      } else if (runinfo->fRunNo >= 2724) {
         adc32_rev = 11;
      } else if (runinfo->fRunNo >= 1694) {
         adc32_rev = 1;
      }

      fAgAsm = new AgAsm();

      fAgAsm->fConfAdc32Rev = adc32_rev;
      
      fAgAsm->fAdcMap = fCfm->ReadFile("adc", "map", runinfo->fRunNo);
      printf("Loaded adc map: %s\n", join(", ", fAgAsm->fAdcMap).c_str());
      
      fAgAsm->fFeamBanks = fCfm->ReadFile("feam", "banks", runinfo->fRunNo);
      printf("Loaded feam banks: %s\n", join(", ", fAgAsm->fFeamBanks).c_str());
   }

   void PreEndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackModule::PreEndRun, run %d\n", runinfo->fRunNo);
   }
   
   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackModule::EndRun, run %d\n", runinfo->fRunNo);
   }
   
   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("UnpackModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("ResumeModule, run %d\n", runinfo->fRunNo);
   }

   uint32_t fNextSerialNumber = 0;

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      if (event->event_id != 1)
         return flow;

      if (fFlags->fAgeSec > 0) {
         const time_t now = time(NULL);
         const time_t t = event->time_stamp;
         int dt = now - t;
         //printf("UnpackModule: serial %d, time %d, age %d, date %s", event->serial_number, event->time_stamp, dt, ctime(&t));
         if (dt >= fFlags->fAgeSec) {
            printf("skipping old event, age %d!\n", dt);
            return flow;
         }
      }

#ifdef _TIME_ANALYSIS_ 
      START_TIMER 
#endif    
      assert(fAgAsm != NULL);

      if (event->serial_number != fNextSerialNumber) {
         printf("UnpackModule: event serial number jump from %d to %d, lost %d events\n", fNextSerialNumber, event->serial_number, event->serial_number - fNextSerialNumber);
         fAgAsm->Reset();
      }

      fNextSerialNumber = event->serial_number + 1;

      if (0) {
         return flow;
      }

      AgEvent* e = fAgAsm->UnpackEvent(event);

      if (fFlags->fPrint) {
         printf("Unpacked AgEvent:   ");
         e->Print();
         printf("\n");
      }

      if (fFlags->fSleepUSec) {
         ::usleep(fFlags->fSleepUSec);
      }

      if (0) {
         delete e;
         return flow;
      }

#ifdef _TIME_ANALYSIS_ 
      if (TimeModules) {
         flow = new AgAnalysisReportFlow(flow, "unpack_module(AgAsm)", timer_start);
      }
#endif 

      return new AgEventFlow(flow, e);
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("UnpackModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class UnpackModuleFactory: public TAFactory
{
public:
   UnpackFlags fFlags;

public:
   void Usage()
   {
      printf("UnpackModuleFactory flags:\n");
      printf("--print -- print something about every event\n");
      printf("--skip-old SEC -- skip events older then given age in seconds (for online use)\n");
      printf("--sleep-usec USEC -- sleep after unpacking each event (for debug use)\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("UnpackModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--print") {
            fFlags.fPrint = true;
         } else if (args[i] == "--skip-old") {
            i++;
            fFlags.fAgeSec = atoi(args[i].c_str());
         } else if (args[i] == "--sleep-usec") {
            i++;
            fFlags.fSleepUSec = atoi(args[i].c_str());
         }
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
