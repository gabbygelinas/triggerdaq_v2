//
// reco_module.cxx
//
// reconstruction of TPC data
//

#include <stdio.h>
#include <iostream>

#include "TH2D.h"
#include "TPolyMarker3D.h"

#include "manalyzer.h"
#include "midasio.h"

#include "AgFlow.h"

#include "Signals.hh"
#include "SpacePoints.hh"
#include "PointsFinder.hh"
#include "TSpacePoint.hh"
#include "TPCBase.hh"

#include "TrackViewer.hh"

#include "settings.hh"
#include "TEvent.hh"
extern double gMagneticField;
extern int gVerb;

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
// private:
//    PointsFinder *pf = NULL;
public:
   TH2D* h_aw_padcol;
   TH1D* h_timediff;
   TH1D* h_atimes;
   TH1D* h_ptimes;
   TH1D* h_firsttimediff;
   TH2D* h_firsttime;

   TH2D* h_times_aw_p;
   TH2D* h_times_aw_p_match;
   TH1D* h_timediff_aw_p_match;

   TH2D *h_unmatched;
   TH1D *h_unmatched_a;
   TH1D *h_unmatched_p;
   TH2D *h_unmatched_p_v_NP;
   TH2D *h_unmatched_p_v_NA;
   TH2D *h_unmatched_a_v_NP;
   TH2D *h_unmatched_a_v_NA;

   TCanvas *cTimes;

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
      // signals = new Signals(runinfo->fRunNo);
      // pf = new PointsFinder(runinfo->fRunNo);
      // pf->SetNoiseThreshold(100);
      // pf->SetMatchPadThreshold(0.1);
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("signalAnalysis");
      dir->cd();
      TPCBase::TPCBaseInstance()->SetPrototype(true);
      int naw = TPCBase::TPCBaseInstance()->GetNumberOfAnodeWires();
      int nps = TPCBase::TPCBaseInstance()->GetNumberPadsRow();
      h_aw_padcol = new TH2D("h_aw_padcol", "anode pad coincidences;anode;pad sector", naw, 0, naw, nps, 0, nps);
      h_timediff = new TH1D("h_timediff", "pad / anode time difference;t_a-t_p",4000,-2000,2000);
      h_atimes = new TH1D("h_atimes", "anode times;t in ns",1000,0,10000);
      h_ptimes = new TH1D("h_ptimes", "pad times;t in ns",500,0,8000);
      h_firsttimediff = new TH1D("h_firsttimediff", "pad / anode first time difference;t_a0-t_p0",4000,-2000,2000);
      h_firsttime = new TH2D("h_firsttime", "pad / anode first time;t_a0;t_p0",1000,0,10000,500,0,8000);

      h_times_aw_p = new TH2D("h_times_aw_p", "pad / anode times;t_a;t_p",1000,0,10000,500,0,8000);
      h_times_aw_p_match = new TH2D("h_times_aw_p_match", "pad / anode times matched;t_a;t_p",1000,0,10000,500,0,8000);
      h_timediff_aw_p_match = new TH1D("h_timediff_aw_p_match", "pad / anode time difference;t_a-t_p",4000,-2000,2000);

      h_unmatched = new TH2D("h_unmatched", "Unmatched signals;number of unmatched anode signals;number of unmatched pad signals",1000,0,1000,1000,0,1000);
      h_unmatched_a = new TH1D("h_unmatched_a", "Unmatched anode signals ratio;number of unmatched anode signals / total number of anode signals",1000,0,1);
      h_unmatched_p = new TH1D("h_unmatched_p", "Unmatched pad signals ratio;number of unmatched pad signals / total number of pad signals",1000,0,1);
      h_unmatched_a_v_NA = new TH2D("h_unmatched_a_v_NA", "Unmatched anode signals vs NA;total number of anode signals;number of unmatched anode signals",1000,0,1000,1000,0,1000);
      h_unmatched_a_v_NP = new TH2D("h_unmatched_a_v_NP", "Unmatched anode signals vs NP;total number of pad signals;number of unmatched anode signals",1000,0,1000,1000,0,1000);
      h_unmatched_p_v_NA = new TH2D("h_unmatched_p_v_NA", "Unmatched pad signals vs NA;total number of anode signals;number of unmatched pad signals",1000,0,1000,1000,0,1000);
      h_unmatched_p_v_NP = new TH2D("h_unmatched_p_v_NP", "Unmatched pad signals vs NP;total number of pad signals;number of unmatched pad signals",1000,0,1000,1000,0,1000);

      gMagneticField=0.;
      gVerb = 2;
      TLookUpTable::LookUpTableInstance()->SetGas("arco2",0.28);
      TLookUpTable::LookUpTableInstance()->SetB(gMagneticField);

      TrackViewer::TrackViewerInstance()->StartViewer();

      TrackViewer::TrackViewerInstance()->StartDeconv();
      TrackViewer::TrackViewerInstance()->StartCoincView();
      TrackViewer::TrackViewerInstance()->StartHitsView();

      cTimes = new TCanvas("cTimes","aw and pad times");
      h_atimes->SetFillStyle(0);
      h_ptimes->SetFillStyle(0);
      h_atimes->SetLineColor(kBlue);
      h_ptimes->SetLineColor(kRed);
      h_atimes->SetMinimum(0);
      h_atimes->SetMaximum(2);
      h_atimes->Draw();
      h_ptimes->Draw("same");
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("RecoRun::EndRun, run %d\n", runinfo->fRunNo);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
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

      TEvent anEvent( event->serial_number, runinfo->fRunNo );

      // const Signals *signals = anEvent.GetSignals();

      // use:
      //
      // age->feam --- pads data
      // age->a16  --- aw data
      //


      double t_pad_first = 1e6;
      double t_aw_first = 1e6;
      if(age->feam && age->a16){
         if(age->feam->complete && age->a16->complete && !age->feam->error && !age->a16->error){
            anEvent.RecEvent( age );
            // pf->Reset();
            // pf->GetSignals()->Reset(age,10,16);
            h_atimes->Reset();
            h_ptimes->Reset();
            // int ntimes = signals->Analyze(age,1,1);
            // cout << "KKKK " << ntimes << " times: " << signals->sanode.size() << '\t' << signals->spad.size() << endl;
            int nmax = std::max(anEvent.GetSignals()->sanode.size(), anEvent.GetSignals()->spad.size());
            // for(int i = 0; i < nmax; i++){
            //    cout << "KKKK " << ((i<signals->sanode.size())?(signals->sanode[i].t):-1) << '\t' << ((i<signals->spad.size())?(signals->spad[i].t):-1) << endl;
            // }

            int unmatchedAnodeSignals, unmatchedPadSignals;
            anEvent.GetNumberOfUnmatched(unmatchedAnodeSignals,unmatchedPadSignals);
            h_unmatched->Fill(unmatchedAnodeSignals,unmatchedPadSignals);
            int NA = anEvent.GetSignals()->sanode.size();
            int NP = anEvent.GetSignals()->spad.size();
            h_unmatched_a->Fill(double(unmatchedAnodeSignals)/double(NA));
            h_unmatched_p->Fill(double(unmatchedPadSignals)/double(NP));
            h_unmatched_a_v_NA->Fill(NA,unmatchedAnodeSignals);
            h_unmatched_a_v_NP->Fill(NP,unmatchedAnodeSignals);
            h_unmatched_p_v_NA->Fill(NA,unmatchedPadSignals);
            h_unmatched_p_v_NP->Fill(NP,unmatchedPadSignals);

            bool first = true;
            cout << NA << endl;
            for(auto sa: anEvent.GetSignals()->sanode){
               if(sa.t < t_aw_first) t_aw_first = sa.t;
               h_atimes->Fill(sa.t);
               double r, phi;
               TPCBase::TPCBaseInstance()->GetAnodePosition(sa.i, r, phi, true);
               std::pair<int,int> pad = TPCBase::TPCBaseInstance()->FindPad(0, phi);
               for(auto sp: anEvent.GetSignals()->spad){
                  if(first){
                     if(sp.t < t_pad_first) t_pad_first = sp.t;
                     h_ptimes->Fill(sp.t);
                  }
                  // cout << "KKKK " << sa.i << '\t' << sp.sec << endl;
                  h_timediff->Fill(sa.t-sp.t);
                  h_times_aw_p->Fill(sa.t, sp.t);

                  if(sp.sec == pad.first){
                     h_times_aw_p_match->Fill(sa.t, sp.t);
                     h_timediff_aw_p_match->Fill(sa.t-sp.t);
                  }

                  if(abs(sa.t-sp.t) < 16)
                     h_aw_padcol->Fill(sa.i, sp.sec);
               }
               first = false;
            }
            cTimes->Update();

            h_firsttimediff->Fill(t_aw_first-t_pad_first);
            h_firsttime->Fill(t_aw_first,t_pad_first);

         }
      }
      // TrackViewer::TrackViewerInstance()->DrawPoints( pf->GetPoints() );
      // TrackViewer::TrackViewerInstance()->DrawPoints2D(anEvent.GetPointsArray() );
      TrackViewer::TrackViewerInstance()->DrawPoints(anEvent.GetPointsArray() );

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
