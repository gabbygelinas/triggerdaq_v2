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
#include "SpacePoints.hh"
#include "PointsFinder.hh"
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
   PointsFinder *pf = NULL;
public:
   TH2D* h_aw_padcol;
   TH1D* h_timediff;
   TH1D* h_firsttimediff;
   TH2D* h_firsttime;

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
      pf = new PointsFinder(runinfo->fRunNo);
      pf->SetNoiseThreshold(100);
      pf->SetMatchPadThreshold(1);
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("signalAnalysis");
      dir->cd();
      h_aw_padcol = new TH2D("h_aw_padcol", "anode pad coincidences;anode;pad sector", TPCBase::NanodeWires, 0, TPCBase::NanodeWires, TPCBase::npadsec, 0, TPCBase::npadsec);
      h_timediff = new TH1D("h_timediff", "pad / anode time difference;t_a-t_p",4000,-2000,2000);
      h_firsttimediff = new TH1D("h_firsttimediff", "pad / anode first time difference;t_a0-t_p0",4000,-2000,2000);
      h_firsttime = new TH2D("h_firsttime", "pad / anode first time;t_a0;t_p0",1000,0,10000,500,0,8000);
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("RecoRun::EndRun, run %d\n", runinfo->fRunNo);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
      DELETE(pf);
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

      double t_pad_first = 1e6;
      double t_aw_first = 1e6;
      if(age->feam && age->a16){
         signals->Reset(age,10,16);
         pf->Reset();
         int ntimes = signals->Analyze(age,1,1);
         cout << "KKKK " << ntimes << " times: " << signals->sanode.size() << '\t' << signals->spad.size() << endl;
         int nmax = std::max(signals->sanode.size(), signals->spad.size());
         for(int i = 0; i < nmax; i++){
            cout << "KKKK " << ((i<signals->sanode.size())?(signals->sanode[i].t):-1) << '\t' << ((i<signals->spad.size())?(signals->spad[i].t):-1) << endl;
         }
         bool first = true;
         for(auto sa: signals->sanode){
            if(sa.t < t_aw_first) t_aw_first = sa.t;
            for(auto sp: signals->spad){
               if(first)
                  if(sp.t < t_pad_first) t_pad_first = sp.t;
               // cout << "KKKK " << sa.i << '\t' << sp.sec << endl;
               h_timediff->Fill(sa.t-sp.t);
               if(abs(sa.t-sp.t) < 16)
                  h_aw_padcol->Fill(sa.i, sp.sec);
            }
            first = false;
         }

         h_firsttimediff->Fill(t_aw_first-t_pad_first);
         h_firsttime->Fill(t_aw_first,t_pad_first);

         int nxy = pf->FindPointsXY(age->a16);
         cout << "PPPP XY " << nxy << endl;
         int nz = pf->FindPointsZ(age->feam);
         cout << "PPPP  Z " << nz << endl;

         auto fullSigs = signals->MatchPads();
         cout << "PPPP  " << fullSigs.size() << endl;
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
