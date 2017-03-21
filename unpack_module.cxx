//
// unpack_module.cxx
//
// K.Olchanski
//

#include <stdio.h>

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include <vector>
#include <deque>

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "Unpack.h"

//#include "Waveform.h"
//#include "TsSync.h"
//#include "Feam.h"

#include "FeamEVB.h"
#include "AgEVB.h"
#include "Unpack.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

//static const double TSNS = 16.0;
static const double TSNS = 8.0/.99999821102751183809;

class UnpackModule: public TAModuleInterface
{
public:
   void Init(const std::vector<std::string> &args);
   void Finish();
   TARunInterface* NewRun(TARunInfo* runinfo);
};

class UnpackRun: public TARunInterface
{
public:
   UnpackModule* fModule;

   Alpha16EVB* fA16Evb;
   FeamEVB*    fFeamEvb;
   AgEVB*      fAgEvb;
   
   UnpackRun(TARunInfo* runinfo, UnpackModule* m)
      : TARunInterface(runinfo)
   {
      printf("UnpackRun::ctor!\n");
      fModule = m;

      fA16Evb  = NULL;
      fFeamEvb = NULL;
      fAgEvb   = NULL;
   }

   ~UnpackRun()
   {
      printf("UnpackRun::dtor!\n");
      DELETE(fA16Evb);
      DELETE(fFeamEvb);
      DELETE(fAgEvb);
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      fA16Evb  = new Alpha16EVB();
      fFeamEvb = new FeamEVB(MAX_FEAM, 1.0/TSNS*1e9);
      fAgEvb = new AgEVB(100.0*1e6/100.0, 125.0*1e6/100.0);

      fA16Evb->Reset();
      fA16Evb->Configure(runinfo->fRunNo);
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("EndRun, run %d\n", runinfo->fRunNo);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("PauseRun, run %d\n", runinfo->fRunNo);
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

            if (e->complete) {
               if (fAgEvb) {
                  fAgEvb->AddAlpha16Event(e);
                  e = NULL;
               }
            }

            DELETE(e);
         }
      }

      if (fFeamEvb) {
         FeamEvent *e = UnpackFeamEvent(fFeamEvb, event);
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
                  printf("module %2d, cnt %4d, ts_trig: 0x%08x %14.3f usec, ts_incr %14.3f usec\n", m->module, m->cnt, m->ts_trig, m->fTime*1e6, m->fTimeIncr*1e6);
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

            if (e->complete) {
               double ta1 = e->a16->eventTime;
               double ta2 = e->a16->prevEventTime;
               double ta = (ta1-ta2)/1e9;
               double tf = e->feam->timeIncr;
               printf("incr %f %f sec, diff %f ns, count %d\n", ta, tf, (tf-ta)*1e9, e->counter);
            }

            DELETE(e);
         }
      }

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

void UnpackModule::Init(const std::vector<std::string> &args)
{
   printf("UnpackModule::Init!\n");

   for (unsigned i=0; i<args.size(); i++) {
#if 0
      if (args[i] == "--nopads")
         fDoPads = false;
      if (args[i] == "--plot1")
         fPlotPad = atoi(args[i+1].c_str());
#endif
   }
}
   
void UnpackModule::Finish()
{
   printf("UnpackModule::Finish!\n");
}
   
TARunInterface* UnpackModule::NewRun(TARunInfo* runinfo)
{
   printf("UnpackModule::NewRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
   return new UnpackRun(runinfo, this);
}

static TARegisterModule tarm(new UnpackModule);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
