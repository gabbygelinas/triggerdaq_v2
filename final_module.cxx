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
   FinalModule* fModule = NULL;

   TCanvas* fC = NULL;
   TCanvas* fPH = NULL;

   TH1D* h_num_aw_hits;
   TH1D* h_num_pad_hits;
   TH2D* h_num_aw_pad_hits;

   TH1D* h_aw_time;
   TH1D* h_aw_amp;
   TH2D* h_aw_amp_time;

   TH2D* h_aw_aw_hits;
   TH2D* h_aw_aw_time;
   TH2D* h_aw_aw_amp;

   TH1D* h_pad_time;
   TH1D* h_pad_amp;
   TH2D* h_pad_amp_time;

   TH2D *h_pad_amp_pad;
   TH2D *h_pad_time_pad;

   TH2D* h_aw_pad_hits;

   TH2D* h_aw_pad_time;
   TH2D *h_aw_amp_aw;
   TH2D *h_aw_time_aw;

   TH2D* h_aw_pad_time_drift;
   TH2D* h_aw_pad_amp_pc;

   FinalRun(TARunInfo* runinfo, FinalModule* m)
      : TARunInterface(runinfo)
   {
      printf("FinalRun::ctor!\n");
      fModule = m;

      //      fC = new TCanvas();
      fPH = new TCanvas("fPH","Pulseheights",600,1000);
      fPH->Divide(1,2);
   }

   ~FinalRun()
   {
      printf("FinalRun::dtor!\n");
      DELETE(fC);
      DELETE(fPH);
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

      h_aw_time = new TH1D("h_aw_time", "aw hit time", 700, 0, 700);
      h_aw_amp = new TH1D("h_aw_amp", "aw hit pulse height", 170, 0, 17000);
      h_aw_amp_time = new TH2D("h_aw_amp_time", "aw p.h. vs time", 70, 0, 700, 50, 0, 17000);

      h_aw_aw_hits = new TH2D("h_aw_aw_hits", "hits in aw vs aw", 128, -0.5, 128-0.5, 128, -0.5, 128-0.5);
      h_aw_aw_time = new TH2D("h_aw_aw_time", "time in aw vs aw", 70, 0, 700, 70, 0, 700);
      h_aw_aw_amp  = new TH2D("h_aw_aw_amp",  "p.h. in aw vs aw", 50, 0, 17000, 50, 0, 17000);

      h_pad_time = new TH1D("h_pad_time", "pad hit time", 500, 0, 500);
      h_pad_amp = new TH1D("h_pad_amp", "pad hit pulse height", 600, 0, 60000);
      h_pad_amp_time = new TH2D("h_pad_amp_time", "pad p.h vs time", 50, 0, 500, 50, 0, 60000);
      int npads = MAX_FEAM*MAX_FEAM_PAD_COL*MAX_FEAM_PAD_ROWS;
      h_pad_amp_pad = new TH2D("h_pad_amp_pad", "pad p.h vs pad number",npads , -0.5, npads-0.5, 600, 0, 60000);
      h_pad_time_pad = new TH2D("h_pad_time_pad", "pad time vs pad number",npads , -0.5, npads-0.5, 500, 0, 500);
      fPH->cd(2);
      h_pad_amp_pad->Draw();

      h_aw_pad_hits = new TH2D("h_aw_pad_hits", "hits in aw vs hits in pads", 4*8, -0.5, 4*8-0.5, 128, -0.5, 128-0.5);

      h_aw_pad_time = new TH2D("h_aw_pad_time", "time of hits in aw vs pads", 50, 0, 500, 70, 0, 700);
      h_aw_amp_aw = new TH2D("h_aw_amp_aw", "aw p.h vs aw number", 256, -0.5, 256.-0.5, 1700, 0, 17000);
      h_aw_time_aw = new TH2D("h_aw_time_aw", "aw time vs aw number", 256, -0.5, 256.-0.5, 700, 0, 700);
      fPH->cd(1);
      h_aw_amp_aw->Draw();

      h_aw_pad_time_drift = new TH2D("h_aw_pad_time_drift", "time of hits in aw vs pads, drift region", 50, 0, 500, 70, 0, 700);

      h_aw_pad_amp_pc = new TH2D("h_aw_pad_amp_pc", "p.h. of hits in aw vs pads, pc region", 50, 0, 60000, 50, 0, 17000);
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
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();

      //int force_plot = false;

      AgEvent* age = ef->fEvent;

      if (age->a16 && age->feam) {
         static bool first = true;
         static double a16_time0 = 0;
         static double feam_time0 = 0;
         static uint32_t a16_ts0 = 0;
         static uint32_t feam_ts0 = 0;

         double a16_time  = age->a16->eventTime/1e9;
         double feam_time = age->feam->time;

         uint32_t a16_ts  = 0;
         uint32_t feam_ts = 0;

         if (1)
            a16_ts = age->a16->udpPacket[0].eventTimestamp;

         if (age->feam->modules[0])
            feam_ts = age->feam->modules[0]->ts_start;

         if (first) {
            first = false;
            a16_time0 = a16_time;
            feam_time0 = feam_time;
            a16_ts0 = a16_ts;
            feam_ts0 = feam_ts;
         }

         double atr = a16_time-a16_time0;
         double ftr = feam_time-feam_time0;
         double dtr = atr-ftr;

         uint32_t a16_tsr = a16_ts - a16_ts0;
         uint32_t feam_tsr = feam_ts - feam_ts0;

         uint32_t a16_tsr_ns = a16_tsr * 10.0;
         uint32_t feam_tsr_ns = feam_tsr * 8.0;

         int dts_ns = (a16_tsr_ns%0x80000000) - (feam_tsr_ns%0x80000000);

         printf("TTT: %d %d, %f %f, diff %f, ts 0x%08x 0x%08x, ns: %12d %12d, diff %d\n", age->a16->eventNo, age->feam->counter, atr, ftr, dtr, a16_tsr, feam_tsr, a16_tsr_ns, feam_tsr_ns, dts_ns);
      }

      if (eawh && eph) {
         if (1) {
            printf("UUU event %d, time %f, anode wire hits: %d, pad hits: %d\n", ef->fEvent->counter, ef->fEvent->time, (int)eawh->fAwHits.size(), (int)eph->fPadHits.size());
         }

         h_num_aw_hits->Fill(eawh->fAwHits.size());
         h_num_pad_hits->Fill(eph->fPadHits.size());
         h_num_aw_pad_hits->Fill(eph->fPadHits.size(), eawh->fAwHits.size());

         for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
            h_aw_time->Fill(eawh->fAwHits[j].time);
            h_aw_amp->Fill(eawh->fAwHits[j].amp);
            h_aw_amp_time->Fill(eawh->fAwHits[j].time, eawh->fAwHits[j].amp);
            h_aw_amp_aw->Fill(eawh->fAwHits[j].chan, eawh->fAwHits[j].amp);
            h_aw_time_aw->Fill(eawh->fAwHits[j].chan, eawh->fAwHits[j].time);

            for (unsigned k=0; k<eawh->fAwHits.size(); k++) {
               if (k==j)
                  continue;
               h_aw_aw_hits->Fill(eawh->fAwHits[j].chan, eawh->fAwHits[k].chan);
               h_aw_aw_time->Fill(eawh->fAwHits[j].time, eawh->fAwHits[k].time);
               h_aw_aw_amp->Fill(eawh->fAwHits[j].amp, eawh->fAwHits[k].amp);
            }
         }

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            h_pad_time->Fill(eph->fPadHits[i].time);
            h_pad_amp->Fill(eph->fPadHits[i].amp);
            h_pad_amp_time->Fill(eph->fPadHits[i].time, eph->fPadHits[i].amp);
            h_pad_amp_pad->Fill(eph->fPadHits[i].col*MAX_FEAM_PAD_ROWS + eph->fPadHits[i].row, eph->fPadHits[i].amp);
            h_pad_time_pad->Fill(eph->fPadHits[i].col*MAX_FEAM_PAD_ROWS + eph->fPadHits[i].row, eph->fPadHits[i].time);
         }

         for (unsigned i=0; i<eph->fPadHits.size(); i++) {
            for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
               int xcol = (eph->fPadHits[i].ifeam%8)*4 + eph->fPadHits[i].col;
               h_aw_pad_hits->Fill(xcol, eawh->fAwHits[j].chan);
               h_aw_pad_time->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);

               if ((eawh->fAwHits[j].time > 200) && eph->fPadHits[i].time > 200) {
                  h_aw_pad_time_drift->Fill(eph->fPadHits[i].time, eawh->fAwHits[j].time);
               }

               if ((eawh->fAwHits[j].time < 200) && eph->fPadHits[i].time < 200) {
                  h_aw_pad_amp_pc->Fill(eph->fPadHits[i].amp, eawh->fAwHits[j].amp);
               }
            }
         }
      }

      //hamp[ichan]->Fill(wamp);

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (do_plot) {
         fPH->GetPad(1)->Modified();
         fPH->GetPad(2)->Modified();
         fPH->Update();
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
