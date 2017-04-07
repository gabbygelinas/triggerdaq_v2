//
// reco_module.cxx
//
// reconstruction of TPC data
//

#include <stdio.h>
#include <iostream>

#include "TH2D.h"

#include "manalyzer.h"
#include "midasio.h"

#include "AgFlow.h"

#include "Signals.hh"
#include "TPCBase.hh"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

using std::cout;
using std::cerr;
using std::endl;

class RecoModule: public TAModuleInterface
{
public:
   void Init(const std::vector<std::string> &args);
   void Finish();
   TARunInterface* NewRun(TARunInfo* runinfo);
};

class RecoRun: public TARunInterface
{
private:
   Signals *signals = NULL;
public:
   TH2D* h_aw_padcol;

   RecoRun(TARunInfo* runinfo)
      : TARunInterface(runinfo)
   {
      printf("RecoRun::ctor!\n");
   }

   ~RecoRun()
   {
      printf("RecoRun::dtor!\n");
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("RecoRun::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      signals = new Signals(runinfo->fRunNo);
      h_aw_padcol = new TH2D("h_aw_padcol", "anode pad coincidences;anode;pad sector", TPCBase::NanodeWires, 0, TPCBase::NanodeWires, TPCBase::npadsec, 0, TPCBase::npadsec);
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("RecoRun::EndRun, run %d\n", runinfo->fRunNo);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
      // DELETE(signals);
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("RecoRun::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("RecoRun::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      // printf("RecoRun::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();

      //int force_plot = false;

      AgEvent* age = ef->fEvent;

      // use:
      //
      // age->feam --- pads data
      // age->a16  --- aw data
      //

      if(age->feam && age->a16){
         signals->Reset(age,10,16);
         int ntimes = signals->Analyze(age,1,1);
         cout << "KKKK " << ntimes << " times: " << signals->sanode.size() << '\t' << signals->spad.size() << endl;
         for(auto sa: signals->sanode){
            for(auto sp: signals->spad){
               // cout << "KKKK " << sa.i << '\t' << sp.sec << endl;
               if(abs(sa.t-sp.t) < 50)
                  h_aw_padcol->Fill(sa.i, sp.sec);
            }
         }
      }
      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("RecoRun::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

void RecoModule::Init(const std::vector<std::string> &args)
{
   printf("RecoModule::Init!\n");

   //fDoPads = true;
   //fPlotPad = -1;
   //fPlotPadCanvas = NULL;

   for (unsigned i=0; i<args.size(); i++) {
      //if (args[i] == "--nopads")
      //   fDoPads = false;
      //if (args[i] == "--plot1")
      //   fPlotPad = atoi(args[i+1].c_str());
   }
}

void RecoModule::Finish()
{
   printf("RecoModule::Finish!\n");
}

TARunInterface* RecoModule::NewRun(TARunInfo* runinfo)
{
   printf("RecoModule::NewRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
   return new RecoRun(runinfo);
}

static TARegisterModule tarm(new RecoModule);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
