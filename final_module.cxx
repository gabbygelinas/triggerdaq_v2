//
// final_module.cxx
//
// final analysis of TPC data
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
#include "AgFlow.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

class FinalModule: public TAModuleInterface
{
public:
   void Init(const std::vector<std::string> &args);
   void Finish();
   TARunInterface* NewRun(TARunInfo* runinfo);
   
   //bool fDoPads;
   //int  fPlotPad;
   //TCanvas* fPlotPadCanvas;
};

class FinalRun: public TARunInterface
{
public:
   FinalModule* fModule;

   TCanvas* fC;

   TH1D* h_num_aw_hits;
   TH1D* h_num_pad_hits;
   TH2D* h_num_aw_pad_hits;

   TH2D* h_aw_pad_hits;

   TH2D* h_aw_pad_time;

#if 0
   TProfile* hbmean_prof;
   TProfile* hbrms_prof;

   TH2D* h2led2amp;

   TH1D* hnhits;
   TH1D* hled_hit;
   TH1D* hamp_hit;
   TH1D* hamp_hit_pedestal;
#endif
   
   FinalRun(TARunInfo* runinfo, FinalModule* m)
      : TARunInterface(runinfo)
   {
      printf("FinalRun::ctor!\n");
      fModule = m;

      fC = new TCanvas();
   }

   ~FinalRun()
   {
      printf("FinalRun::dtor!\n");
      DELETE(fC);
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* dir = gDirectory->mkdir("final");
      dir->cd(); // select correct ROOT directory

      dir->mkdir("summary")->cd();

      h_num_aw_hits = new TH1D("h_num_aw_hits", "number of anode wire hits", 100, 0, 100);
      h_num_pad_hits = new TH1D("h_num_pad_hits", "number of cathode pad hits", 100, 0, 100);
      h_num_aw_pad_hits = new TH2D("h_num_aw_pad_hits", "number of aw vs pad hits", 50, 0, 100, 100, 0, 50);

      h_aw_pad_hits = new TH2D("h_aw_pad_hits", "hits in aw vs hits in pads", 2432, -0.5, 2432-0.5, 128, -0.5, 128-0.5);

      h_aw_pad_time = new TH2D("h_aw_pad_time", "time of hits in aw vs pads", 50, 0, 500, 70, 0, 700);

#if 0
      hbrms_all  = new TH1D("hbrms",  "baseline rms",  100, 0, 200);

      hbmean_prof = new TProfile("hbmean_prof", "baseline mean vs channel", nchan, -0.5, nchan-0.5);
      hbrms_prof  = new TProfile("hbrms_prof",  "baseline rms vs channel",  nchan, -0.5, nchan-0.5);

      hamp_all   = new TH1D("hamp",   "pulse height", 100, 0, 17000);
      hled_all   = new TH1D("hled",   "pulse leading edge, adc time bins", 100, 0, 900);

      h2led2amp  = new TH2D("h2led2amp", "pulse amp vs time, adc time bins", 100, 0, 900, 100, 0, 17000);

      hled_all_cut = new TH1D("hled_cut",   "pulse leading edge, adc time bins, with cuts", 100, 0, 900);
      hamp_all_cut = new TH1D("hamp_cut",   "pulse height, with cuts", 100, 0, 17000);

      hnhits = new TH1D("hnhits", "hits per channel", nchan, -0.5, nchan-0.5);
      hled_hit = new TH1D("hled_hit", "hit time, adc time bins", 100, 0, 900);
      hamp_hit = new TH1D("hamp_hit", "hit pulse height", 100, 0, 17000);
      hamp_hit_pedestal = new TH1D("hamp_hit_pedestal", "hit pulse height, zoom on pedestal", 100, 0, 300);

      // FIXME: who deletes this?
      hbmean = new TH1D*[nchan];
      hbrms  = new TH1D*[nchan];
      hwaveform = new TH1D*[nchan];
      hamp = new TH1D*[nchan];
      hled = new TH1D*[nchan];
   
      pads->mkdir("baseline_mean")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbmean%04d", i);
         sprintf(title, "chan %04d baseline mean", i);
         hbmean[i] = new TH1D(name, title, 100, 0, 17000);
      }

      pads->mkdir("baseline_rms")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hbrms%04d", i);
         sprintf(title, "chan %04d baseline rms", i);
         hbrms[i] = new TH1D(name, title, 100, 0, 200);
      }

      pads->mkdir("chan_waveform")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hwaveform%04d", i);
         sprintf(title, "chan %04d waveform", i);
         hwaveform[i] = new TH1D(name, title, nbins, -0.5, nbins-0.5);
      }

      pads->mkdir("chan_amp")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hamp%04d", i);
         sprintf(title, "chan %04d pulse height", i);
         hamp[i] = new TH1D(name, title, 100, 0, 17000);
      }

      pads->mkdir("chan_led")->cd();

      for (int i=0; i<nchan; i++) {
         char name[256];
         char title[256];
         sprintf(name, "hled%04d", i);
         sprintf(title, "chan %04d pulse leading edge, adc bins", i);
         hled[i] = new TH1D(name, title, 100, 0, 900);
      }
#endif
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

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();

      if (!eawh) {
         return flow;
      }

      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();

      if (!eph) {
         return flow;
      }

      //int force_plot = false;

      if (1) {
         printf("UUU event %d, time %f, anode wire hits: %d, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size(), (int)eph->fPadHits.size());
      }

      h_num_aw_hits->Fill(eawh->fAwHits.size());
      h_num_pad_hits->Fill(eph->fPadHits.size());
      h_num_aw_pad_hits->Fill(eph->fPadHits.size(), eawh->fAwHits.size());

      for (unsigned i=0; i<eph->fPadHits.size(); i++) {
         for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
            h_aw_pad_hits->Fill(eph->fPadHits[i].chan, eawh->fAwHits[j].chan);
            h_aw_pad_time->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);
         }
      }
         
      //hamp[ichan]->Fill(wamp);

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
#if 0
         // plot waveforms

         fC->Clear();
         fC->Divide(2,3);
         
         if (1) {
            fC->cd(1);
            TH1D* hh = new TH1D("hh", "hh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hh->SetBinContent(ibin+1, e->adcs[0]->adc[0][0][ibin]);
            }
            hh->Draw();
         }
         
         if (1) {
            fC->cd(2);
            TH1D* hhh = new TH1D("hhh", "hhh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hhh->SetBinContent(ibin+1, e->adcs[0]->adc[0][0][ibin]);
            }
            hhh->SetMinimum(-33000);
            hhh->SetMaximum(+33000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(3);
            int nbins = ww[39]->nsamples;
            TH1D* hhh = new TH1D("hhhh", "hhhh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[39]->samples[ibin]);
            }
            hhh->SetMinimum(-9000);
            hhh->SetMaximum(+9000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(4);
            int nbins = ww[iplot]->nsamples;
            TH1D* hhh = new TH1D("hhhhh", "hhhhh", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               hhh->SetBinContent(ibin+1, ww[iplot]->samples[ibin]);
            }
            hhh->SetMinimum(-9000);
            hhh->SetMaximum(+9000);
            hhh->Draw();
         }
         
         if (1) {
            fC->cd(5);
            int nbins = ww[33]->nsamples;
            TH1D* h33 = new TH1D("h33", "h33", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               h33->SetBinContent(ibin+1, ww[33]->samples[ibin]);
            }
            h33->SetMinimum(-9000);
            h33->SetMaximum(+9000);
            h33->Draw();
         }
         
         if (1) {
            fC->cd(6);
            int nbins = ww[34]->nsamples;
            TH1D* h34 = new TH1D("h34", "h34", nbins, 0, nbins);
            for (int ibin=0; ibin<nbins; ibin++) {
               h34->SetBinContent(ibin+1, ww[34]->samples[ibin]);
            }
            h34->SetMinimum(-9000);
            h34->SetMaximum(+9000);
            h34->Draw();
         }
         
         fC->Modified();
         fC->Draw();
         fC->Update();
#endif
      }

#if 0
      if (fModule->fPlotPad >= 0) {
         if (!fModule->fPlotPadCanvas)
            fModule->fPlotPadCanvas = new TCanvas("FEAM PAD", "FEAM PAD", 900, 650);

         TCanvas*c = fModule->fPlotPadCanvas;

         c->cd();

         int nbins = ww[fModule->fPlotPad]->nsamples;
         TH1D* h = new TH1D("h", "h", nbins, 0, nbins);
         for (int ibin=0; ibin<nbins; ibin++) {
            h->SetBinContent(ibin+1, ww[fModule->fPlotPad]->samples[ibin]);
         }

         h->SetMinimum(-9000);
         h->SetMaximum(+9000);
         h->Draw();
         
         c->Modified();
         c->Draw();
         c->Update();
      }
#endif

#if 0
      time_t now = time(NULL);

      if (force_plot) {
         static time_t plot_next = 0;
         if (now > plot_next) {
            //fATX->PlotA16Canvas();
            plot_next = time(NULL) + 15;
         }
      }

      static time_t t = 0;

      if (now - t > 15) {
         t = now;
         //fATX->Plot();
      }

      *flags |= TAFlag_DISPLAY;
#endif

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

void FinalModule::Init(const std::vector<std::string> &args)
{
   printf("FinalModule::Init!\n");

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
   
void FinalModule::Finish()
{
   printf("FinalModule::Finish!\n");

   //DELETE(fPlotPadCanvas);
}
   
TARunInterface* FinalModule::NewRun(TARunInfo* runinfo)
{
   printf("FinalModule::NewRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
   return new FinalRun(runinfo, this);
}

static TARegisterModule tarm(new FinalModule);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
