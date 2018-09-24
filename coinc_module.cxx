//
// coinc_module.cxx
//
// analysis of wires/pads coincidence
//
// A.Capra
//

#include <stdio.h>
#include "manalyzer.h"
#include "midasio.h"
    
#include "AgFlow.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

class CoincModule: public TARunObject
{
public:
   bool fTrace = true;
   double fCoincTime;

   TH2D* h_aw_pad_num_hits;
   TH2D* h_aw_pad_hits;
   TH2D* h_aw_pad_time;

   TH1D* h_coinc;
   TH2D* h_padrow_awamp;
   TH1D* h_padrow_awamp_px;
   // Profile histograms are used to display the mean value of Y and its error for each bin in X. 
   TProfile* h_padrow_awamp_pfx;

   TH2D* h_padrow_awamp_pc;
   TH1D* h_padrow_awamp_pc_px;
   TProfile* h_padrow_awamp_pc_pfx;

   CoincModule(TARunInfo* runinfo): TARunObject(runinfo),
                                    fCoincTime(16.)
   {
      if(fTrace)
         printf("CoincModule::ctor!\n");
   }
  
   ~CoincModule()
   {
      if(fTrace)
         printf("CoincModule::dtor!\n");
   }
  
   void BeginRun(TARunInfo* runinfo)
   {
      if(fTrace)
         printf("CoincModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("coinc");
      dir->cd(); // select correct ROOT directory

      h_aw_pad_num_hits = new TH2D("h_aw_pad_num_hits", "number of aw vs pad hits; number if hits in aw; number of hits in pads", 
                                   50, 0., 250., 50, 0., 250.);
      h_aw_pad_hits = new TH2D("h_aw_pad_hits", "hits in aw vs hits in pads; tpc wire; pad column", 
                               256, 0., 256, 32., 0., 32.);
      h_aw_pad_time = new TH2D("h_aw_pad_time", "time of hits in aw vs pads; time in aw, ns; time in pads, ns", 
                               50, 1000., 6000., 50, 1000., 6000.);

      h_coinc = new TH1D("h_coinc","Coincidence AW*Pads;Number;Events",100,0.,100.);
      h_padrow_awamp = new TH2D("h_padrow_awamp","AW amplitude vs Z;Pad Row;ADC counts",576,0.,576.,
                                300,0.,30000.);
      h_padrow_awamp_pc = new TH2D("h_padrow_awamp_pc","AW amplitude vs Z  Proportional Region;Pad Row;ADC counts",
                                   576,0.,576.,300,0.,30000.);


   }

   void EndRun(TARunInfo* runinfo)
   {
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      gDirectory->cd("coinc"); // select correct ROOT directory
    
      h_padrow_awamp_px = h_padrow_awamp->ProjectionX();
      h_padrow_awamp_px->SetMinimum(0.);
      TString ptitle(h_padrow_awamp->GetTitle());
      ptitle+="   Projection";
      h_padrow_awamp_px->SetTitle(ptitle.Data());

      h_padrow_awamp_pfx = h_padrow_awamp->ProfileX();
      h_padrow_awamp_pfx->SetMinimum(0.);
      TString pftitle(h_padrow_awamp->GetTitle());
      pftitle+="  Profile;";
      pftitle+=h_padrow_awamp->GetXaxis()->GetTitle();
      pftitle+=";Average ADC counts";
      h_padrow_awamp_pfx->SetTitle(pftitle.Data());


      h_padrow_awamp_pc_px = h_padrow_awamp_pc->ProjectionX();
      h_padrow_awamp_pc_px->SetMinimum(0.);
      ptitle=h_padrow_awamp_pc->GetTitle();
      ptitle+="   Projection";
      h_padrow_awamp_pc_px->SetTitle(ptitle.Data());

      h_padrow_awamp_pc_pfx = h_padrow_awamp_pc->ProfileX();
      h_padrow_awamp_pc_pfx->SetMinimum(0.);
      pftitle = h_padrow_awamp_pc->GetTitle();
      pftitle+="  Profile;";
      pftitle+=h_padrow_awamp_pc->GetXaxis()->GetTitle();
      pftitle+=";Average ADC counts";
      h_padrow_awamp_pc_pfx->SetTitle(pftitle.Data());

      if(fTrace)
         printf("CoincModule::EndRun, run %d\n", runinfo->fRunNo);
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if(fTrace)
         printf("CoincModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if(fTrace)
         printf("CoincModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();


      if(eawh && eph)
         {
            if( fTrace ) 
               printf("coinc event %d, time %f, anode wire hits: %d, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, 
                      (int)eawh->fAwHits.size(), (int)eph->fPadHits.size());
	
            //      h_aw_pad_num_hits->Fill(eawh->fAwHits.size(), eph->fPadHits.size());
	
            h_aw_pad_num_hits->Fill(eawh->fAwHits.size(), eph->fPadHits.size());
            double counter=0.;
            for( unsigned i=0; i<eph->fPadHits.size(); i++ )
               {
                  int col = eph->fPadHits[i].tpc_col;
                  col+=1;
                  if( col == 32 ) col = 0;
                  assert(col<32);
                  for( unsigned j=0; j<eawh->fAwHits.size(); j++ )
                     {
                        int aw = eawh->fAwHits[j].wire >= 256 ? eawh->fAwHits[j].wire-256 : eawh->fAwHits[j].wire;
                        h_aw_pad_hits->Fill(aw, col);
                        h_aw_pad_time->Fill(eawh->fAwHits[j].time, eph->fPadHits[i].time_ns);
                        int sec = aw/8;
                        if( sec != col ) continue;
	      
                        if( fabs(eawh->fAwHits[j].time-eph->fPadHits[i].time_ns) > fCoincTime ) continue;
	      
                        h_padrow_awamp->Fill( eph->fPadHits[i].tpc_row, eawh->fAwHits[j].amp );

                        if( eawh->fAwHits[j].time > 800. && eawh->fAwHits[j].time < 1200. ) // ns
                           h_padrow_awamp_pc->Fill( eph->fPadHits[i].tpc_row, eawh->fAwHits[j].amp );

                        ++counter;
                     }// loop aw
               } // loop pads
            h_coinc->Fill(counter);
         } // if aw and pad hits

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("CoincModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }

};

class CoincModuleFactory: public TAFactory
{
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("CoincModuleFactory::Init!\n");
      for (unsigned i=0; i<args.size(); i++) 
         {    }
   }

   void Finish()
   {
      printf("CoincModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("CoincModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new CoincModule(runinfo);
   }
};

static TARegister tar(new CoincModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
