//
// MIDAS analyzer example 2: ROOT analyzer
//
// K.Olchanski
//

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "Alpha16.h"
#include "Unpack.h"
#include "AgFlow.h"

// ALPHA-g common analysis code

#define SHOW_ALPHA16 16

#include <iostream>
#include "TVirtualFFT.h"
#include "Waveform.h"

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

#if 0
Waveform* NewWaveform(const Alpha16Waveform* a, double scale)
{
   Waveform* w = new Waveform(a->size());
   for (unsigned s=0; s<a->size(); s++) {
      w->samples[s] = (*a)[s] * scale;
   }
   return w;
}
#endif

Waveform* NewWaveform(const std::vector<int>* a, double scale)
{
   Waveform* w = new Waveform(a->size());
   for (unsigned s=0; s<a->size(); s++) {
      w->samples[s] = (*a)[s] * scale;
   }
   return w;
}

struct PlotHistograms
{
   TCanvas* fCanvas;

   TH1D* fHbaseline;
   TH1D* fHbaselineRms;
   TProfile* fHbaselineRmsVsChan;
   TH1D* fHph;
   TH1D* fHle;
   TH1D* fHph_r1;
   TH1D* fHph_r2;
   TH1D* fHph_r3;
   TH1D* fHle_r1;
   TH1D* fHle_r2;
   TH1D* fHle_r3;
   TH1D* fHlex;
   TH1D* fHocc;

   TH1D* fHph1;
   TH1D* fHph2;

   TH2D* fHph3;

   TH1D* fHocc1;
   TH1D* fHocc2;

   TH1D* fHtdiff;
   TH1D* fHtdiff2;

   TProfile* fHph2occ1;
   TProfile* fHph2occ2;

   TH1D* fHleCal;
   TH1D* fHoccCal;
   TProfile* fHleVsChanCal;

   PlotHistograms(TCanvas* c) // ctor
   {
      if (!c) {
         c = new TCanvas("Histograms", "Histograms", 1100, 850);
         if (!(c->GetShowEventStatus()))
            c->ToggleEventStatus();
         if (!(c->GetShowToolBar()))
            c->ToggleToolBar();
      }

      fCanvas = c;

      fCanvas->cd();
      fCanvas->Divide(3,6);

      int max_adc = 9000;

      int i=1;

      fCanvas->cd(i++);
      fHbaseline = new TH1D("baseline", "baseline", 100, -1000, 1000);
      fHbaseline->Draw();

      fCanvas->cd(i++);
      fHbaselineRms = new TH1D("baseline_rms", "baseline_rms", 100, 0, 200);
      fHbaselineRms->Draw();

      fCanvas->cd(i++);
      fHbaselineRmsVsChan = new TProfile("baseline_rms_vs_chan", "baseline_rms_vs_chan", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHbaselineRmsVsChan->SetMinimum(0);
      fHbaselineRmsVsChan->Draw();

      fCanvas->cd(i++);
      fHph = new TH1D("pulse_height", "pulse_height", 100, 0, max_adc);
      fHph->Draw();

      fHph_r1 = new TH1D("pulse_height_r1", "pulse_height_r1", 100, 0, max_adc);
      fHph_r2 = new TH1D("pulse_height_r2", "pulse_height_r2", 100, 0, max_adc);
      fHph_r3 = new TH1D("pulse_height_r3", "pulse_height_r3", 100, 0, max_adc);
      
      fCanvas->cd(i++);
      fHle = new TH1D("pulse_time", "pulse_time", 100, 0, 1000);
      fHle->Draw();

      fHle_r1 = new TH1D("pulse_time_r1", "pulse_time_r1", 100, 0, 1000);
      fHle_r2 = new TH1D("pulse_time_r2", "pulse_time_r2", 100, 0, 1000);
      fHle_r3 = new TH1D("pulse_time_r3", "pulse_time_r3", 100, 0, 1000);

      fCanvas->cd(i++);
      fHlex = new TH1D("pulse_time_expanded", "pulse_time_expanded", 100, 100, 200);
      fHlex->Draw();

      fCanvas->cd(i++);
      fHocc = new TH1D("channel_occupancy", "channel_occupancy", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHocc->SetMinimum(0);
      fHocc->Draw();

      fCanvas->cd(i++);
      fHph1 = new TH1D("pulse_height_pc", "pulse_height_pc", 100, 0, max_adc);
      fHph1->Draw();

      fCanvas->cd(i++);
      fHph2 = new TH1D("pulse_height_drift", "pulse_height_drift", 100, 0, max_adc);
      fHph2->Draw();

      fCanvas->cd(i++);
      fHph3 = new TH2D("pulse_height_vs_drift", "pulse_height_vs_drift", 50, 0, 1000, 50, 0, max_adc);
      fHph3->Draw();

      fCanvas->cd(i++);
      fHocc1 = new TH1D("channel_occupancy_pc", "channel_occupancy_pc", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHocc1->SetMinimum(0);
      fHocc1->Draw();

      fCanvas->cd(i++);
      fHocc2 = new TH1D("channel_occupancy_drift", "channel_occupancy_drift", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHocc2->SetMinimum(0);
      fHocc2->Draw();

      fCanvas->cd(i++);
      fHtdiff = new TH1D("time_between_events", "time_between_events", 2000, 0, 2.0);
      fHtdiff2 = new TH1D("time_between_events_expanded", "time_between_events_expanded", 1000, 0, 10.0);
      fHtdiff2->Draw();

      fCanvas->cd(i++);
      fHph2occ1 = new TProfile("pulse_height_profile_pc", "pulse_height_profile_pc", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHph2occ1->SetMinimum(0);
      fHph2occ1->Draw();

      fCanvas->cd(i++);
      fHph2occ2 = new TProfile("pulse_height_profile_drift", "pulse_height_profile_drift", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHph2occ2->SetMinimum(0);
      fHph2occ2->Draw();

      fCanvas->cd(i++);
      fHleCal = new TH1D("pulse_time_cal", "pulse_time_cal", 100, 100, 120);
      fHleCal->Draw();

      fCanvas->cd(i++);
      fHoccCal = new TH1D("channel_occupancy_cal", "channel_occupancy_cal", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHoccCal->SetMinimum(0);
      fHoccCal->Draw();

      fCanvas->cd(i++);
      fHleVsChanCal = new TProfile("pulse_time_vs_chan_cal", "pulse_time_vs_chan_cal", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHleVsChanCal->Draw();

      Draw();
   }

   ~PlotHistograms() // dtor
   {
      if (fCanvas)
         delete fCanvas;
   }

   void Draw()
   {
      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};

#if 0
struct PlotHistogramsPads
{
   TCanvas* fCanvas;

   TH1D* fHbaseline;
   TH1D* fHbaselineRms;
   TH1D* fHph;
   TH1D* fHle;
   TH1D* fHlex;
   TH1D* fHocc;

   TH1D* fHocc1;
   TH1D* fHocc2;

   TH1D* fHph1;
   TH1D* fHph2;

   TH2D* fHph3;

   PlotHistogramsPads(TCanvas* c) // ctor
   {
      if (!c) {
         c = new TCanvas("HistogramsPads", "HistogramsPads", 1100, 850);
         if (!(c->GetShowEventStatus()))
            c->ToggleEventStatus();
         if (!(c->GetShowToolBar()))
            c->ToggleToolBar();
      }

      fCanvas = c;

      fCanvas->cd();
      fCanvas->Divide(3,4);

      fCanvas->cd(1);
      fHbaseline = new TH1D("pads_baseline", "baseline", 100, -1000, 1000);
      fHbaseline->Draw();

      fCanvas->cd(2);
      fHbaselineRms = new TH1D("pads_baseline_rms", "baseline_rms", 50, 0, 50);
      fHbaselineRms->Draw();

      fCanvas->cd(3);
      fHph = new TH1D("pads_pulse_height", "pulse_height", 100, 0, 500);
      fHph->Draw();

      fCanvas->cd(4);
      fHle = new TH1D("pads_pulse_time", "pulse_time", 100, 0, 1000);
      fHle->Draw();

      fCanvas->cd(5);
      fHlex = new TH1D("pads_pulse_time_expanded", "pulse_time_expanded", 100, 100, 200);
      fHlex->Draw();

      fCanvas->cd(6);
      fHocc = new TH1D("pads_channel_occupancy", "channel_occupancy", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
      fHocc->Draw();

      fCanvas->cd(7);
      fHph1 = new TH1D("pads_pulse_height_pc", "pulse_height_pc", 100, 0, 500);
      fHph1->Draw();

      fCanvas->cd(8);
      fHph2 = new TH1D("pads_pulse_height_drift", "pulse_height_drift", 100, 0, 500);
      fHph2->Draw();

      fCanvas->cd(9);
      fHph3 = new TH2D("pads_pulse_height_vs_drift", "pulse_height_vs_drift", 50, 0, 1000, 50, 0, 500);
      fHph3->Draw();

      fCanvas->cd(10);
      fHocc1 = new TH1D("pads_channel_occupancy_pc", "channel_occupancy_pc", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
      fHocc1->Draw();

      fCanvas->cd(11);
      fHocc2 = new TH1D("pads_channel_occupancy_drift", "channel_occupancy_drift", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
      fHocc2->Draw();

      Draw();
   }

   ~PlotHistogramsPads() // dtor
   {
      if (fCanvas)
         delete fCanvas;
   }

   void Draw()
   {
      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};
#endif

class A16ChanHistograms
{
public:
   std::string fNameBase;
   std::string fTitleBase;
   int fNbins = 0;

   TH1D* hbmean = NULL;
   TH1D* hbrms  = NULL;

   TH1D* hwaveform_first = NULL;
   TH1D* hwaveform_max   = NULL;
   TH1D* hwaveform_max_drift = NULL;
   TH1D* hwaveform_avg   = NULL;
   TH1D* hwaveform_avg_drift = NULL;

   int nwf = 0;
   int nwf_drift = 0;
   double fMaxWamp = 0;
   double fMaxWampDrift = 0;


public:
   A16ChanHistograms(const char* xname, const char* xtitle, TDirectory* dir, int nbins) // ctor
   {
      printf("Create name [%s] title [%s] with %d bins\n", xname, xtitle, nbins);
      
      TDirectory* dir_first = dir->GetDirectory("achan_waveform_first");
      if(!dir_first) dir_first = dir->mkdir("achan_waveform_first");
      TDirectory* dir_max = dir->GetDirectory("achan_waveform_max");
      if(!dir_max) dir_max = dir->mkdir("achan_waveform_max");
      TDirectory* dir_max_drift = dir->GetDirectory("achan_waveform_max_drift");
      if(!dir_max_drift) dir_max_drift = dir->mkdir("achan_waveform_max_drift");
      TDirectory* dir_avg = dir->GetDirectory("achan_waveform_avg");
      if(!dir_avg) dir_avg = dir->mkdir("achan_waveform_avg");
      TDirectory* dir_avg_drift = dir->GetDirectory("achan_waveform_avg_drift");
      if(!dir_avg_drift) dir_avg_drift = dir->mkdir("achan_waveform_avg_drift");

      fNameBase = xname;
      fTitleBase = xtitle;
      fNbins = nbins;

      char name[256];
      char title[256];

      sprintf(name, "hawf_first_%s", xname);
      sprintf(title, "%s first waveform", xtitle);

      dir_first->cd();
      hwaveform_first = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hawf_max_%s", xname);
      sprintf(title, "%s biggest waveform", xtitle);
      dir_max->cd();
      hwaveform_max = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hawf_max_drift_%s", xname);
      sprintf(title, "%s biggest waveform, drift region", xtitle);
      dir_max_drift->cd();
      hwaveform_max_drift = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hawf_avg_%s", xname);
      sprintf(title, "%s average waveform", xtitle);
      dir_avg->cd();
      hwaveform_avg = new TH1D(name, title, nbins, -0.5, nbins-0.5);

      sprintf(name, "hawf_avg_drift_%s", xname);
      sprintf(title, "%s average waveform, drift region", xtitle);
      dir_avg_drift->cd();
      hwaveform_avg_drift = new TH1D(name, title, nbins, -0.5, nbins-0.5);
   }

   ~A16ChanHistograms() // dtor
   {

   }

};

struct PlotA16
{
   TCanvas* fCanvas = NULL;
   std::vector<TH1D*> fH;

   int fModule = 0;
   int fSection = 0;

   PlotA16(TCanvas* c, int module, int section) // ctor
   {
      if (!c) {
         char title[256];
         sprintf(title, "ADC mod%d section %d", module, section);
         c = new TCanvas(title, title, 900, 650);
         //c = new TCanvas("ALPHA16 ADC", "ALPHA16 ADC", 900, 650);
         if (!(c->GetShowEventStatus()))
            c->ToggleEventStatus();
         if (!(c->GetShowToolBar()))
            c->ToggleToolBar();
      }

      fCanvas = c;

      fCanvas->cd();
      fCanvas->Divide(4,4);

      fModule = module;
      fSection = section;

      for (int i=0; i<16; i++) {
         fH.push_back(NULL);
      }
   }

   ~PlotA16()
   {
      if (fCanvas) {
         delete fCanvas;
         fCanvas = NULL;
      }
   }

   void DrawChannel(int i, const Alpha16Channel* c)
   {
      assert(i>=0 && i<16);

      fCanvas->cd(i+1);

      int color = 1;
      
      if (!fH[i]) {
         char name[256];
         sprintf(name, "mod%dch%d", c->adc_module, c->adc_chan);
         
         fH[i] = new TH1D(name, name, c->adc_samples.size(), 0, c->adc_samples.size());
         
         fH[i]->Draw();
         fH[i]->SetMinimum(-(1<<15));
         fH[i]->SetMaximum(1<<15);
         //fH[i]->SetMinimum(-2000);
         //fH[i]->SetMaximum(2000);
         fH[i]->GetYaxis()->SetLabelSize(0.10);
         //fH[i]->SetLineColor(color);
      }

      fH[i]->SetLineColor(color);

      for (unsigned s=0; s<c->adc_samples.size(); s++) {
         int v = c->adc_samples[s];
         //printf("bin %d, v %d\n", i, v);
         fH[i]->SetBinContent(s+1, v);
      }
   }

   void ClearChannel(int i)
   {
      if (fH[i]) {
         fH[i]->SetLineColor(2);
         //fH[i]->Clear();
      }
   }
   
   void Draw(const Alpha16Event* e)
   {
      // colors:
      // 0 = white
      // 1 = black
      // 2 = red
      // 3 = green
      // 4 = blue

      for (int i=0; i<16; i++) {
         ClearChannel(i);
      }

      for (unsigned i=0; i<e->hits.size(); i++) {
         Alpha16Channel* c = e->hits[i];
         if (c->adc_module != fModule)
            continue;
         if (fSection == 0) { // 16 onboard ADCs
            if (c->adc_chan < 0 || c->adc_chan >= 16)
               continue;
            DrawChannel(c->adc_chan, c);
         } else if (fSection == 1) {
            if (c->adc_chan < 16 || c->adc_chan >= 32)
               continue;
            DrawChannel(c->adc_chan-16, c);
         } else if (fSection == 2) {
            if (c->adc_chan < 32 || c->adc_chan >= 48)
               continue;
            DrawChannel(c->adc_chan-32, c);
         }
      }

      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};

#if 0
struct PlotWire
{
   TCanvas* fCanvas;
   TH1D* hh[3];

   PlotWire(TCanvas* c) // ctor
   {
      if (!c) {
         c = new TCanvas("ALPHA-g One-Wire TPC prototype", "ALPHA-g One-Wire TPC prototype", 900, 500);
      }

      fCanvas = c;

      fCanvas->cd();
      fCanvas->Divide(1, 1);

      for (int i=0; i<3; i++) {
         hh[i] = NULL;
      }
   }

   void Fill(const Alpha16Event* e)
   {
#if 0
      int diff = t2 - t1;
      printf("slot %d, jitter %3d: %3d %3d\n", slot, diff, t1, t2);

      if (hg[slot])
         hg[slot]->Fill(diff);
#endif
   }

   void Draw(const Alpha16Event* e)
   {
      // colors:
      // 0 = white
      // 1 = black
      // 2 = red
      // 3 = green
      // 4 = blue

      if (1) {
         fCanvas->cd(1);

         int color = 1;
         for (int i=0; i<3; i++, color++) {
            int ch = i;

            if (!hh[i]) {
               char name[256];
               sprintf(name, "hh%d", i);
               hh[i] = new TH1D(name, name, e->waveform[ch].size(), 0, e->waveform[ch].size());

               if (color == 1)
                  hh[i]->Draw();
               else
                  hh[i]->Draw("same");
               hh[i]->SetMinimum(-(1<<15));
               hh[i]->SetMaximum(1<<15);
               hh[i]->SetLineColor(color);
            }

            for (unsigned s=0; s<e->waveform[ch].size(); s++)
               hh[i]->SetBinContent(s+1, e->waveform[ch][s]);
         }
      }

#if 0
      if (1) {
         fCanvas->cd(1+slot*3+2);

         int ii = slot;
         if (!hg[ii]) {
            char name[256];
            sprintf(name, "hg%d", ii);
            hg[ii] = new TH1D(name, name, 21, -10, 10);
         }

         hg[ii]->Draw();
      }
#endif
   }

   void Update()
   {
      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};
#endif

#if 0
struct FileWire
{
   FILE *fp1;
   FILE *fp2;
   FILE *fp3;

   int chan1;
   int chan2;
   int chan3;

   int max_write;
   int count;

   FileWire(int x_max_write) // ctor
   {
      chan1 = 0;
      chan2 = 1;
      chan3 = 2;
      count = 0;
      max_write = x_max_write;
   }

   void Open(int runno)
   {
      char buf[2560];
      const char* dir = "/pool8tb/agdaq/wire/text";

      sprintf(buf, "%s/run%05dchan%02d.txt", dir, runno, chan1);

      fp1 = fopen(buf, "w");
      assert(fp1);

      sprintf(buf, "%s/run%05dchan%02d.txt", dir, runno, chan2);

      fp2 = fopen(buf, "w");
      assert(fp2);

      sprintf(buf, "%s/run%05dchan%02d.txt", dir, runno, chan3);

      fp3 = fopen(buf, "w");
      assert(fp3);

      count = 0;
   }

   void Close()
   {
      fclose(fp1);
      fclose(fp2);
      fclose(fp3);

      fp1 = NULL;
      fp2 = NULL;
      fp3 = NULL;

      printf("Write %d events\n", count);
   }

   void PrintWave(FILE* fp, int chan, const Alpha16Waveform* w)
   {
      fprintf(fp, "%d", chan);
      for (unsigned s=0; s<w->size(); s++)
         fprintf(fp, ",%d", (int)((*w)[s]));
      fprintf(fp, "\n");
   }

   void Fill(const Alpha16Event* e)
   {
      PrintWave(fp1, chan1, &e->waveform[chan1]);
      PrintWave(fp2, chan2, &e->waveform[chan2]);
      PrintWave(fp3, chan3, &e->waveform[chan3]);
      count++;

      if (max_write)
         if (count > max_write) {
            Close();
            exit(1);
         }
   }

};
#endif

#if 0
class AlphaTpcX
{
public:
   int fRunNo;

   Alpha16EVB* fEvb;

   //PlotA16* fPlotA16[SHOW_ALPHA16];

   //Alpha16Event* fLastEvent;

   //PlotHistograms* fH;
   //PlotHistogramsPads* fP;

   int fCountEarlyBad;
   int fCountGood;
   int fCountBad;

   TDirectory *dnoise = NULL;
   TDirectory *dwf = NULL;
   TCanvas *c;
   TH1D *hwf0, *hwf1, *hfftsum0, *hfftsum1, *hbase0, *hbase1, *hRMS0, *hRMS1;
   std::vector<A16ChanHistograms*> fHC;

   int entries = 0;

   AlphaTpcX()
   {
      //fH = new PlotHistograms(NULL);
      //fP = NULL; // new PlotHistogramsPads(NULL);

      fEvb = new Alpha16EVB();

      dnoise = gDirectory->mkdir("noise");
      dwf = gDirectory->mkdir("waveforms");
      dnoise->cd();
      c = new TCanvas("cnoise","noise analysis",800,1200);
      c->Divide(1,6);
      hwf0 = new TH1D("hwf0","waveform 0",701,0,701);
      c->cd(1);
      hwf0->Draw();
      c->GetPad(2)->SetLogy();

      hwf1 = new TH1D("hwf1","waveform 1",701,0,701);
      c->cd(3);
      hwf1->Draw();
      c->GetPad(4)->SetLogy();

      hfftsum0 = new TH1D("hfftsum0","FFT sum",701,0,701);
      c->cd(5);
      hfftsum0->Draw();
      c->GetPad(5)->SetLogy();

      hfftsum1 = new TH1D("hfftsum1","FFT sum",701,0,701);
      c->cd(6);
      hfftsum1->Draw();
      c->GetPad(6)->SetLogy();

      hbase0 = new TH1D("hbase0", "baseline 0", 2000, -1000, 1000);
      hbase1 = new TH1D("hbase1", "baseline 1", 2000, -1000, 1000);

      hRMS0 = new TH1D("hRMS0", "RMS 0", 1000, 0, 100);
      hRMS1 = new TH1D("hRMS1", "RMS 1", 1000, 0, 100);

      //for (int i=0; i<SHOW_ALPHA16; i++)
      //fPlotA16[i] = NULL;

      //fLastEvent = NULL;

      fCountEarlyBad = 0;
      fCountGood = 0;
      fCountBad = 0;
   }

   ~AlphaTpcX()
   {
      if (fEvb) {
         delete fEvb;
         fEvb = NULL;
      }

#if 0
      for (int i=0; i<SHOW_ALPHA16; i++) {
         if (fPlotA16[i]) {
            delete fPlotA16[i];
            fPlotA16[i] = NULL;
         }
      }
#endif

      //if (fLastEvent) {
      //   delete fLastEvent;
      //   fLastEvent = NULL;
      //}

#if 0
      if (fH) {
         delete fH;
         fH = NULL;
      }

      if (fP) {
         delete fP;
         fP = NULL;
      }
#endif
   }

#if 0
   void CreateA16Canvas()
   {
      for (int i=0; i<SHOW_ALPHA16; i++) {
         char title[256];
         sprintf(title, "ADC %2d", i+1);
         TCanvas *c = new TCanvas(title, title, 900, 400);
         fPlotA16[i] = new PlotA16(c, i*NUM_CHAN_ALPHA16);
      }
   }
#endif

   void BeginRun(int run)
   {
      fRunNo = run;
      fEvb->Reset();
      fEvb->Configure(run);
   }

   void EndRun()
   {
      printf("AlphaTpcX::EndRun: early bad events: %d, good events: %d, bad events: %d, total %d events\n", fCountEarlyBad, fCountGood, fCountBad, fCountEarlyBad + fCountGood + fCountBad);
      hfftsum0->Scale(1./double(hbase0->GetEntries()*sqrt(701.)));
      hfftsum1->Scale(1./double(hbase1->GetEntries()*sqrt(701.)));
      for(auto *hc: fHC){
         if(hc->nwf) hc->hwaveform_avg->Scale(1./double(hc->nwf));
         if(hc->nwf_drift) hc->hwaveform_avg_drift->Scale(1./double(hc->nwf_drift));
      }
   }

   int Event(Alpha16Event* e, std::vector<AgAwHit>* hits)
   {
      bool verbose = false;
      
      //printf("new event: "); e->Print();
      //e->error = false;

      if (e->error || !e->complete) {
         if (fCountGood)
            fCountBad++;
         else
            fCountEarlyBad++;
         return 0;
      }

      fCountGood++;

      if (0) {
         //printf("event: "); e->Print();
         printf("a16: event %4d, ts %14.3f usec, ts_incr %14.3f usec\n", e->counter, e->time*1e6, e->timeIncr*1e6);
         return 0;
      }

      printf("Have Alpha16 event: "); e->Print(); printf("\n");

      //if (fLastEvent) {
      //   delete fLastEvent;
      //   fLastEvent = NULL;
      //}
      //
      //fLastEvent = e;

      bool force_plot = false;

      //printf("analyzing: "); e->Print();

      // analyze the event

      int nhits = 0;

#if 0
      for (int i=0; i<SHOW_ALPHA16*NUM_CHAN_ALPHA16; i++)
         if (e->udpPresent[i]) {

            //if (i != 76)
            //   continue;

            Waveform* w = NewWaveform(&e->waveform[i], 1.0/4.0);

#if 0
            if(i == 96 || i == 127)
               AnalyzeNoise(w, i);
#endif

            double b, brms;
            b = baseline(w, 0, 100, NULL, &brms);

            double wmin = min(w);
            double wmax = max(w);

            if (0 && i>=0 && i<32) {
               fP->fHbaseline->Fill(b);
               fP->fHbaselineRms->Fill(brms);

               int pph = wmax - b;

               if (pph > 30) {

                  fP->fHph->Fill(pph);

                  if (pph > 100)
                     force_plot = true;

                  int ple = led(w, b, 1.0, pph/2.0);

                  if (1 || pph > 100) {
                     fP->fHle->Fill(ple);
                     fP->fHlex->Fill(ple);
                     fP->fHocc->Fill(i);

                     if (ple > 150 && ple < 180) {
                        fP->fHocc1->Fill(i);
                        fP->fHph1->Fill(pph);
                     }

                     if (ple > 180 && ple < 580) {
                        fP->fHocc2->Fill(i);
                        fP->fHph2->Fill(pph);
                     }

                     fP->fHph3->Fill(ple, pph);
                  }
               }

               delete w;
               continue;
            }

            if (0 && i < 48) {
               delete w;
               continue;
            }

            fH->fHbaseline->Fill(b);
            fH->fHbaselineRms->Fill(brms);
            fH->fHbaselineRmsVsChan->Fill(i, brms);

            double ph = b - wmin;

            ////// Plot waveforms
            if (i >= fHC.size()) {
               for (unsigned j=fHC.size(); j<=i; j++)
                  fHC.push_back(NULL);
            }

            if (fHC[i] == NULL){
               char xname[256];
               char xtitle[256];
               sprintf(xname, "hawf_%03d", i);
               sprintf(xtitle, "Waveform anode %03d", i);
               fHC[i] = new A16ChanHistograms(xname, xtitle, dwf, w->nsamples);
            }

            // save first waveform

            bool doPrint = false;
            if (fHC[i]->hwaveform_first->GetEntries() == 0) {
               if (doPrint)
                  printf("saving first waveform %d\n", i);
               for (int j=0; j< w->nsamples; j++)
                  fHC[i]->hwaveform_first->SetBinContent(j+1, w->samples[j]);
            }

            // save biggest waveform

            if (ph > fHC[i]->fMaxWamp) {
               fHC[i]->fMaxWamp = ph;
               if (doPrint)
                  printf("saving biggest waveform %d\n", i);
               for (int j=0; j< w->nsamples; j++)
                  fHC[i]->hwaveform_max->SetBinContent(j+1, w->samples[j]);
            }

            // add to average waveform

            for (int j=0; j< w->nsamples; j++)
               fHC[i]->hwaveform_avg->AddBinContent(j+1, w->samples[j]);
            fHC[i]->nwf++;

            ////////

            //if (ph > 120) {
            if (ph > 250) {
               fH->fHph->Fill(ph);

               int le = led(w, b, -1.0, ph/2.0);

               fH->fHle->Fill(le);
               fH->fHlex->Fill(le);
               fH->fHocc->Fill(i);

               if (le > 150 && le < 180) {
                  fH->fHocc1->Fill(i);
                  fH->fHph1->Fill(ph);

                  if (ph < 7000)
                     fH->fHph2occ1->Fill(i, ph);
               }

               if (le > 180 && le < 580) {
                  fH->fHocc2->Fill(i);
                  fH->fHph2->Fill(ph);

                  if (ph < 7000)
                     fH->fHph2occ2->Fill(i, ph);
               }

               fH->fHph3->Fill(le, ph);

               int calStart = 160;
               int calEnd = 180;

               if (le > calStart && le < calEnd) {
                  fH->fHleCal->Fill(le);
                  fH->fHoccCal->Fill(i);
                  fH->fHleVsChanCal->Fill(i, le);
               }

               if (ph > 4000) {
                  nhits++;
                  //printf("samples %d %d, ", e->waveform[i].size(), w->nsamples);
                  if (verbose)
                     printf("chan %4d: baseline %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %4d\n", i, b, brms, wmin, wmax, ph, le);
               }

               bool hit = false;

               if (le > 150 && le < 580 && ph > 600) {
                  hit = true;
               }

               if (hit) {
                  AgAwHit h;
                  h.chan = i;
                  h.time = le;
                  h.amp = ph;
                  hits->push_back(h);
               }
                           // save biggest drift region waveform

               if (le > 180){
                  if(ph > fHC[i]->fMaxWampDrift) {
                     fHC[i]->fMaxWampDrift = ph;
                     if (doPrint)
                        printf("saving biggest drift waveform %d\n", i);
                     for (int j=0; j< w->nsamples; j++)
                        fHC[i]->hwaveform_max_drift->SetBinContent(j+1, w->samples[j]);
                  }
                  // add to average waveform

                  for (int j=0; j< w->nsamples; j++)
                     fHC[i]->hwaveform_avg_drift->AddBinContent(j+1, w->samples[j]);
                  fHC[i]->nwf_drift++;
               }


            }

            delete w;
         }
#endif

#if 0
      double tdiff = e->timeIncr;
      //printf("time diff %f sec\n", tdiff);

      fH->fHtdiff->Fill(tdiff);
      fH->fHtdiff2->Fill(tdiff);
#endif

      if (verbose)
         printf("Found %d hits\n", nhits);

      if (nhits >= 5)
         force_plot = true;

      return force_plot;
   }

#if 0
   void Plot()
   {
      time_t tstart = time(NULL);

      if (fH) {
         printf("plot H start\n");
         fH->Draw();
         printf("plot H done\n");
      }

      if (fP) {
         printf("plot P start\n");
         fP->Draw();
         printf("plot P done\n");
      }

      time_t tend = time(NULL);
      int elapsed = tend-tstart;

      printf("plotting: done, %d sec!\n", elapsed);
   }
#endif

#if 0
   void PlotA16Canvas(Alpha16Event* e)
   {
      time_t tstart = time(NULL);

      printf("plotting:   "); e->Print();

      for (int i=0; i<SHOW_ALPHA16; i++) {
         if (fPlotA16[i]) {
            printf("plot %d start\n", i);
            fPlotA16[i]->Draw(e);
            printf("plot %d done\n", i);
         }
      }

      time_t tend = time(NULL);
      int elapsed = tend-tstart;

      printf("plotting: done, %d sec!\n", elapsed);
   }
#endif

   void AnalyzeNoise(Waveform *w, short i){
      TH1D *h, *hRMS, *hbase, *hfftsum;
      short index;
      switch(i){
      case 96: h = hwf0; hbase = hbase0; hRMS = hRMS0; hfftsum = hfftsum0; index = 2; break;
      case 127: h = hwf1; hbase = hbase1; hRMS = hRMS1; hfftsum = hfftsum1; index = 4; break;
      default: std::cerr << "Something's wrong!" << std::endl;
      }
      for(int b = 0; b < w->nsamples; b++)
         h->SetBinContent(b+1, w->samples[b]);

      double b, brms;
      b = baseline(w, 0, w->nsamples, NULL, &brms);
      hbase->Fill(b);
      hRMS->Fill(brms);

      TH1 *hfft = 0;
      TVirtualFFT::SetTransform(0);
      hfft = h->FFT(hfft, "MAG");
      c->cd(index);
      hfft->Draw();
      hfftsum->Add(hfft);
      entries++;
      for(int p = 1; p <= 6; p++)
         c->GetPad(p)->Modified();
      c->Update();
   }
};
#endif

class A16Flags
{
public:
   bool fPlotWF = false;
   bool fDoPlotAll = false;
};

class A16Module: public TARunObject
{
public:
   A16Flags* fFlags = NULL;
   int fCounter = 0;
   //AlphaTpcX* fATX = NULL;
   bool fTrace = false;

   PlotHistograms* fH;
   std::vector<A16ChanHistograms*> fHC;
   std::vector<PlotA16*> fPlotA16;

   TDirectory *dnoise = NULL;
   TDirectory *dwf = NULL;

public:
   A16Module(TARunInfo* runinfo, A16Flags* f)
      : TARunObject(runinfo)
   {
      fTrace = false;
      if (fTrace)
         printf("A16Module::ctor!\n");
      fFlags = f;

      runinfo->fRoot->fOutputFile->cd();
      TDirectory* aw = gDirectory->mkdir("aw");
      aw->cd(); // select correct ROOT directory

      fH = new PlotHistograms(NULL);

      dnoise = gDirectory->mkdir("noise");
      dwf = gDirectory->mkdir("waveforms");

      //fATX = new AlphaTpcX();
#if 0
      if (fFlags->fPlotWF)
         fATX->CreateA16Canvas();
#endif
   }

   ~A16Module()
   {
      if (fTrace)
         printf("A16Module::dtor!\n");
      if (fH) {
         delete fH;
         fH = NULL;
      }
      //if (fATX)
      //   delete fATX;
      //fATX = NULL;
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("A16Module::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      fCounter = 0;
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      //fATX->BeginRun(runinfo->fRunNo);

      fPlotA16.push_back(new PlotA16(NULL, 1, 0));
      fPlotA16.push_back(new PlotA16(NULL, 2, 0));
      fPlotA16.push_back(new PlotA16(NULL, 3, 0));
      fPlotA16.push_back(new PlotA16(NULL, 4, 0));
      fPlotA16.push_back(new PlotA16(NULL, 5, 0));
      fPlotA16.push_back(new PlotA16(NULL, 6, 0));
      fPlotA16.push_back(new PlotA16(NULL, 6, 1));
      fPlotA16.push_back(new PlotA16(NULL, 6, 2));
      fPlotA16.push_back(new PlotA16(NULL, 7, 0));
      fPlotA16.push_back(new PlotA16(NULL, 8, 0));

      fPlotA16.push_back(new PlotA16(NULL, 9, 0));
      fPlotA16.push_back(new PlotA16(NULL, 10, 0));
      fPlotA16.push_back(new PlotA16(NULL, 11, 0));
      fPlotA16.push_back(new PlotA16(NULL, 12, 0));
      fPlotA16.push_back(new PlotA16(NULL, 13, 0));
      fPlotA16.push_back(new PlotA16(NULL, 14, 0));
      fPlotA16.push_back(new PlotA16(NULL, 15, 0));
      fPlotA16.push_back(new PlotA16(NULL, 16, 0));
      fPlotA16.push_back(new PlotA16(NULL, 17, 0));
      fPlotA16.push_back(new PlotA16(NULL, 18, 0));
   }

   void EndRun(TARunInfo* runinfo)
   {
      printf("A16Module::EndRun, run %d, events %d\n", runinfo->fRunNo, fCounter);
      time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
      //fATX->EndRun();
      //char fname[1024];
      //sprintf(fname, "output%05d.pdf", runinfo->fRunNo);
      //fATX->fH->fCanvas->SaveAs(fname);
   }

   void PauseRun(TARunInfo* runinfo)
   {
      printf("A16Module::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("A16Module::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void AnalyzeHit(const Alpha16Channel* hit, std::vector<AgAwHit>* flow_hits)
   {
      int i = 0;
      int r = 1;

      if (hit->adc_module == 1) {
         i = 16*0 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 2) {
         i = 16*1 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 3) {
         i = 16*2 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 4) {
         i = 16*3 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 5) {
         i = 16*4 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 6) {
         if (hit->adc_chan < 16) {
            i = 16*5 + hit->adc_chan;
            r = 1;
         } else {
            // 32ch ADC
            i = 16*8 + hit->adc_chan - 16;
            r = 2;
         }
      } else if (hit->adc_module == 7) {
         i = 16*6 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 8) {
         i = 16*7 + hit->adc_chan;
         r = 1;
      } else if (hit->adc_module == 11) {
         i = 16*10 + hit->adc_chan;
         r = 3;
      } else if (hit->adc_module == 12) {
         i = 16*11 + hit->adc_chan;
         r = 3;
      } else if (hit->adc_module == 9) {
         i = 16*12 + hit->adc_chan;
         r = 3;
      } else if (hit->adc_module == 14) {
         i = 16*13 + hit->adc_chan;
         r = 3;
      } else if (hit->adc_module == 15) {
         i = 16*14 + hit->adc_chan;
         r = 3;
      } else if (hit->adc_module == 18) {
         i = 16*15 + hit->adc_chan;
         r = 3;
      } else {
         i = -1;
      }
      
      printf("hit: bank [%s], adc_module %d, adc_chan %d, preamp_pos %d, preamp_wire %d, tpc_wire %d. nbins %d+%d, seqno %d, region %d\n",
             hit->bank.c_str(),
             hit->adc_module,
             hit->adc_chan,
             hit->preamp_pos,
             hit->preamp_wire,
             hit->tpc_wire,
             hit->first_bin,
             (int)hit->adc_samples.size(),
             i,
             r);

      if (i < 0)
         return;

      Waveform* w = NewWaveform(&hit->adc_samples, 1.0);
      
      ////// Plot waveforms
      
      if (i >= fHC.size()) {
         for (unsigned j=fHC.size(); j<=i; j++)
            fHC.push_back(NULL);
      }

      if (fHC[i] == NULL){
         char xname[256];
         char xtitle[256];
         sprintf(xname, "hawf_%03d", i);
         sprintf(xtitle, "Waveform seqno %03d", i);
         fHC[i] = new A16ChanHistograms(xname, xtitle, dwf, w->nsamples);
      }

#if 0
      if(i == 96 || i == 127)
         AnalyzeNoise(w, i);
#endif

      double b, brms;
      b = baseline(w, 0, 100, NULL, &brms);
      
      double wmin = min(w);
      double wmax = max(w);

      fH->fHbaseline->Fill(b);
      fH->fHbaselineRms->Fill(brms);
      fH->fHbaselineRmsVsChan->Fill(i, brms);
      
      double ph = b - wmin;

      // save first waveform
      
      bool doPrint = false;
      if (fHC[i]->hwaveform_first->GetEntries() == 0) {
         if (doPrint)
            printf("saving first waveform %d\n", i);
         for (int j=0; j< w->nsamples; j++)
            fHC[i]->hwaveform_first->SetBinContent(j+1, w->samples[j]);
      }
      
      // save biggest waveform
      
      if (ph > fHC[i]->fMaxWamp) {
         fHC[i]->fMaxWamp = ph;
         if (doPrint)
            printf("saving biggest waveform %d\n", i);
         for (int j=0; j< w->nsamples; j++)
            fHC[i]->hwaveform_max->SetBinContent(j+1, w->samples[j]);
      }
      
      // add to average waveform
      
      for (int j=0; j< w->nsamples; j++)
         fHC[i]->hwaveform_avg->AddBinContent(j+1, w->samples[j]);
      fHC[i]->nwf++;

      double ph_hit_thr = 250;

      if (r == 1)
         ph_hit_thr = 1000;
      else if (r == 2)
         ph_hit_thr = 500;
      
      if (ph > ph_hit_thr) {
         fH->fHph->Fill(ph);

         if (r==1)
            fH->fHph_r1->Fill(ph);
         else if (r==2)
            fH->fHph_r2->Fill(ph);
         else if (r==3)
            fH->fHph_r3->Fill(ph);
         
         int le = led(w, b, -1.0, ph/2.0);
         
         fH->fHle->Fill(le);

         if (r==1)
            fH->fHle_r1->Fill(le);
         else if (r==2)
            fH->fHle_r2->Fill(le);
         else if (r==3)
            fH->fHle_r3->Fill(le);

         fH->fHlex->Fill(le);
         fH->fHocc->Fill(i);
         
         if (le > 150 && le < 180) {
            fH->fHocc1->Fill(i);
            fH->fHph1->Fill(ph);
            
            if (ph < 7000)
               fH->fHph2occ1->Fill(i, ph);
         }
         
         if (le > 180 && le < 580) {
            fH->fHocc2->Fill(i);
            fH->fHph2->Fill(ph);
            
            if (ph < 7000) {
               fH->fHph2occ2->Fill(i, ph);
            }
         }
         
         fH->fHph3->Fill(le, ph);
         
         int calStart = 160;
         int calEnd = 180;
         
         if (le > calStart && le < calEnd) {
            fH->fHleCal->Fill(le);
            fH->fHoccCal->Fill(i);
            fH->fHleVsChanCal->Fill(i, le);
         }
         
         //if (ph > 4000) {
         //nhits++;
         //printf("samples %d %d, ", e->waveform[i].size(), w->nsamples);
         if (1) {
            printf("chan %4d: baseline %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %4d\n", i, b, brms, wmin, wmax, ph, le);
         }
         //}
         
         bool have_hit = false;
         
         if (le > 150 && le < 580 && ph > 600) {
            have_hit = true;
         }
         
         if (have_hit) {
            AgAwHit h;
            h.chan = i;
            h.time = le;
            h.amp = ph;
            flow_hits->push_back(h);
         }

         // save biggest drift region waveform
         
         if (le > 180) {
            if(ph > fHC[i]->fMaxWampDrift) {
               fHC[i]->fMaxWampDrift = ph;
               if (doPrint)
                  printf("saving biggest drift waveform %d\n", i);
               for (int j=0; j< w->nsamples; j++)
                  fHC[i]->hwaveform_max_drift->SetBinContent(j+1, w->samples[j]);
            }

            // add to average waveform
            
            for (int j=0; j< w->nsamples; j++)
               fHC[i]->hwaveform_avg_drift->AddBinContent(j+1, w->samples[j]);
            fHC[i]->nwf_drift++;
         }
      }
      
      delete w;
   }
   
   TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      if (fTrace)
         printf("A16Module::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      Alpha16Event* e = ef->fEvent->a16;

      if (!e) {
         return flow;
      }

      printf("---> in a16module: ");
      e->Print();
      printf("\n");

      for (unsigned i=0; i<fPlotA16.size(); i++) {
         fPlotA16[i]->Draw(e);
      }

      AgAwHitsFlow* flow_hits = new AgAwHitsFlow(flow);
      flow = flow_hits;

      for (unsigned i=0; i<e->hits.size(); i++) {
         AnalyzeHit(e->hits[i], &flow_hits->fAwHits);
      }

      *flags |= TAFlag_DISPLAY;

      //int force_plot = fATX->Event(e, &hits->fAwHits);
      //if (fFlags->fDoPlotAll)
      //force_plot = true;

      time_t now = time(NULL);

#if 0
      if (force_plot) {
         static time_t plot_next = 0;
         if (now > plot_next) {
            fATX->PlotA16Canvas(e);
            plot_next = time(NULL) + 15;
         }
      }
#endif

#if 0
      static time_t t = 0;

      if (now - t > 15) {
         t = now;
         fATX->Plot();
      }
#endif

      fCounter++;

#if 0
      delete e;
#endif

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class A16ModuleFactory: public TAFactory
{
public:
   A16Flags fFlags;
   
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("A16ModuleFactory::Init!\n");
      
      fFlags.fPlotWF = false;
      
      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--wf")
            fFlags.fPlotWF = true;
         if (args[i] == "--wfall") {
            fFlags.fDoPlotAll = true;
            fFlags.fPlotWF = true;
         }
      }
      
      TARootHelper::fgDir->cd(); // select correct ROOT directory
   }
   
   void Finish()
   {
      printf("A16ModuleFactory::Finish!\n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("A16ModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new A16Module(runinfo, &fFlags);
   }

public:
};

static TARegister tar(new A16ModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
