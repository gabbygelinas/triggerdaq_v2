// ALPHA-g common analysis code

#define SHOW_ALPHA16 8

#include "Waveform.h"

Waveform* NewWaveform(const Alpha16Waveform* a, double scale)
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
   TH1D* fHlex;
   TH1D* fHocc;

   TH1D* fHph1;
   TH1D* fHph2;

   TH2D* fHph3;

   TH1D* fHocc1;
   TH1D* fHocc2;

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
      fHbaselineRms = new TH1D("baseline_rms", "baseline_rms", 100, 0, 100);
      fHbaselineRms->Draw();

      fCanvas->cd(i++);
      fHbaselineRmsVsChan = new TProfile("baseline_rms_vs_chan", "baseline_rms_vs_chan", SHOW_ALPHA16*NUM_CHAN_ALPHA16, -0.5, SHOW_ALPHA16*NUM_CHAN_ALPHA16-0.5);
      fHbaselineRmsVsChan->SetMinimum(0);
      fHbaselineRmsVsChan->Draw();

      fCanvas->cd(i++);
      fHph = new TH1D("pulse_height", "pulse_height", 100, 0, max_adc);
      fHph->Draw();

      fCanvas->cd(i++);
      fHle = new TH1D("pulse_time", "pulse_time", 100, 0, 1000);
      fHle->Draw();

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
      // empty slot

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

static int x = 0;

struct PlotA16
{
   TCanvas* fCanvas;
   TH1D* fH[16];
   int fFirstChan;

   PlotA16(TCanvas* c, int first_chan) // ctor
   {
      if (!c) {
         c = new TCanvas("ALPHA16 ADC", "ALPHA16 ADC", 900, 650);
         if (!(c->GetShowEventStatus()))
            c->ToggleEventStatus();
         if (!(c->GetShowToolBar()))
            c->ToggleToolBar();
      }

      fCanvas = c;

      fCanvas->cd();
      fCanvas->Divide(4,4);

      fFirstChan = first_chan;

      for (int i=0; i<16; i++)
         fH[i] = NULL;
   }

   ~PlotA16()
   {
      if (fCanvas) {
         delete fCanvas;
         fCanvas = NULL;
      }
   }

   void Draw(const Alpha16Event* adc)
   {
      // colors:
      // 0 = white
      // 1 = black
      // 2 = red
      // 3 = green
      // 4 = blue

      for (int i=0; i<16; i++) {
         fCanvas->cd(i+1);
         
         int color = 1;

         int ichan = fFirstChan + i;

         if (!adc->udpPresent[ichan])
            continue;
         
         if (!fH[i]) {
            char name[256];
            sprintf(name, "a16ch%02d_x%d", i, x++);

            fH[i] = new TH1D(name, name, adc->waveform[ichan].size(), 0, adc->waveform[ichan].size());

            fH[i]->Draw();
            fH[i]->SetMinimum(-(1<<15));
            fH[i]->SetMaximum(1<<15);
            //fH[i]->SetMinimum(-2000);
            //fH[i]->SetMaximum(2000);
            fH[i]->GetYaxis()->SetLabelSize(0.10);
            fH[i]->SetLineColor(color);
         }

         for (unsigned s=0; s<adc->waveform[ichan].size(); s++) {
            int v = adc->waveform[ichan][s];
            //printf("bin %d, v %d\n", i, v);
            fH[i]->SetBinContent(s+1, v);
         }
      }

      fCanvas->Modified();
      fCanvas->Draw();
      fCanvas->Update();
   }
};

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

class AlphaTpcX
{
public:
   int fRunNo;
   
   Alpha16EVB* fEvb;
   
   PlotA16* fPlotA16[SHOW_ALPHA16];
   
   Alpha16Event* fLastEvent;

   PlotHistograms* fH;
   PlotHistogramsPads* fP;

   int fCountEarlyBad;
   int fCountGood;
   int fCountBad;

   AlphaTpcX()
   {
      fH = new PlotHistograms(NULL);
      fP = NULL; // new PlotHistogramsPads(NULL);

      fEvb = new Alpha16EVB();

      for (int i=0; i<SHOW_ALPHA16; i++)
         fPlotA16[i] = NULL;

      fLastEvent = NULL;

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

      for (int i=0; i<SHOW_ALPHA16; i++) {
         if (fPlotA16[i]) {
            delete fPlotA16[i];
            fPlotA16[i] = NULL;
         }
      }

      if (fLastEvent) {
         delete fLastEvent;
         fLastEvent = NULL;
      }

      if (fH) {
         delete fH;
         fH = NULL;
      }
         
      if (fP) {
         delete fP;
         fP = NULL;
      }
   }

   void CreateA16Canvas()
   {
      for (int i=0; i<SHOW_ALPHA16; i++) {
         char title[256];
         sprintf(title, "ADC %2d", i+1);
         TCanvas *c = new TCanvas(title, title, 900, 400);
         fPlotA16[i] = new PlotA16(c, i*NUM_CHAN_ALPHA16);
      }
   }

   void BeginRun(int run)
   {
      fRunNo = run;
      fEvb->Reset();
      fEvb->Configure(run);
   }

   void EndRun()
   {
      printf("EndRun: early bad events: %d, good events: %d, bad events: %d, total %d events\n", fCountEarlyBad, fCountGood, fCountBad, fCountEarlyBad + fCountGood + fCountBad);
   }

   int Event(Alpha16Event* e)
   {
      printf("new event: "); e->Print();
      
      if (e->error || !e->complete) {
         if (fCountGood)
            fCountBad++;
         else
            fCountEarlyBad++;
         return 0;
      }

      fCountGood++;
      
      printf("event: "); e->Print();
      
      if (fLastEvent) {
         delete fLastEvent;
         fLastEvent = NULL;
      }
      
      fLastEvent = e;
      
      bool force_plot = false;
      
      printf("analyzing: "); e->Print();
      
      // analyze the event

      int nhits = 0;

      for (int i=0; i<SHOW_ALPHA16*NUM_CHAN_ALPHA16; i++)
         if (e->udpPresent[i]) {

            //if (i != 76)
            //   continue;
            
            Waveform* w = NewWaveform(&e->waveform[i], 1.0/4.0);

            double b, brms;
            b = baseline(w, 100, NULL, &brms);

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

               if (le > 100 && le < 120) {
                  fH->fHleCal->Fill(le);
                  fH->fHoccCal->Fill(i);
                  fH->fHleVsChanCal->Fill(i, le);
               }

               if (ph > 4000) {
                  nhits++;
                  //printf("samples %d %d, ", e->waveform[i].size(), w->nsamples);
                  printf("chan %4d: baseline %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %4d\n", i, b, brms, wmin, wmax, ph, le);
               }
            }

            delete w;
         }

      printf("Found %d hits\n", nhits);

      if (nhits >= 5)
         force_plot = true;

      return force_plot;
   }
      
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

   void PlotA16Canvas()
   {
      if (!fLastEvent)
         return;

      time_t tstart = time(NULL);
      
      printf("plotting:   "); fLastEvent->Print();

      for (int i=0; i<SHOW_ALPHA16; i++) {
         if (fPlotA16[i]) {
            printf("plot %d start\n", i);
            fPlotA16[i]->Draw(fLastEvent);
            printf("plot %d done\n", i);
         }
      }

      time_t tend = time(NULL);
      int elapsed = tend-tstart;

      printf("plotting: done, %d sec!\n", elapsed);
   }
};

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
