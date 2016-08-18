#include <stdio.h>
#include <iostream>

#include "TRootanaDisplay.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "TInterestingEventManager.hxx"

#include "Alpha16.h"
#include "Unpack.h"

#define SHOW_ALPHA16 6

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

   TH1D* fHocc1;
   TH1D* fHocc2;

   TH1D* fHph1;
   TH1D* fHph2;

   TH2D* fHph3;

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
      fCanvas->Divide(3,4);

      //int max_adc = 8000;
      int max_adc = 1000;

      int i=1;

      fCanvas->cd(i++);
      fHbaseline = new TH1D("baseline", "baseline", 100, -1000, 1000);
      fHbaseline->Draw();

      fCanvas->cd(i++);
      fHbaselineRms = new TH1D("baseline_rms", "baseline_rms", 50, 0, 50);
      fHbaselineRms->Draw();

      fCanvas->cd(i++);
      fHbaselineRmsVsChan = new TProfile("baseline_rms_vs_chan", "baseline_rms_vs_chan", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
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
      fHocc = new TH1D("channel_occupancy", "channel_occupancy", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
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
      fHocc1 = new TH1D("channel_occupancy_pc", "channel_occupancy_pc", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
      fHocc1->Draw();

      fCanvas->cd(i++);
      fHocc2 = new TH1D("channel_occupancy_drift", "channel_occupancy_drift", SHOW_ALPHA16*NUM_CHAN_ALPHA16, 0, SHOW_ALPHA16*NUM_CHAN_ALPHA16-1);
      fHocc2->Draw();

      Draw();
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
            //fH[i]->SetMinimum(-(1<<15));
            //fH[i]->SetMaximum(1<<15);
            fH[i]->SetMinimum(-2000);
            fH[i]->SetMaximum(2000);
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

#if 0
class Alpha16Canvas: public TCanvasHandleBase {
public:
   
   int fRunNo;
   
   std::string fBankName;
   
   Alpha16EVB* fEvb;
   
   PlotA16* fPlotA16;
   PlotWire* fPlotWire;
   FileWire* fFileWire;
   
   Alpha16Event* fLastEvent;
   
   Alpha16Canvas(const char* tab_title, const char* bankname, PlotA16* plotA16, PlotWire* plotWire, int fileWrite) // ctor
      : TCanvasHandleBase(tab_title)
   {
      fBankName = bankname;
      fEvb = new Alpha16EVB(NUM_CHAN_ALPHA16, 512);

      fPlotA16 = plotA16;
      fPlotWire = plotWire;
      fFileWire = NULL;

      fLastEvent = NULL;

      if (fileWrite) {
         fFileWire = new FileWire(fileWrite);
      }
   };
   
   ~Alpha16Canvas() // dtor
   {
      printf("Alpha16Canvas dtor!\n");
   };
   
   void BeginRun(int transition,int run,int time)
   {
      printf("Alpha16Canvas::BeginRun() for bank %s\n", fBankName.c_str());
      
      fRunNo = run;

      if (fFileWire)
         fFileWire->Open(run);
      
      fEvb->Reset();
   };
   
   void EndRun(int transition,int run,int time)
   {
      if (fFileWire)
         fFileWire->Close();
   };
   
   void ResetCanvasHistograms()
   {
      printf("Alpha16Canvas::ResetCanvasHistograms()\n");
   };
   
   void UpdateCanvasHistograms(TDataContainer& dataContainer)
   {
      //printf("Alpha16Canvas::UpdateCanvasHistograms()\n");
      
      const TMidasEvent& me = dataContainer.GetMidasEvent();
      
      void *ptr;
      int bklen, bktype;
      int status;
      
      status = me.FindBank(fBankName.c_str(), &bklen, &bktype, &ptr);
      if (status == 1) {
         //printf("ALPHA16 pointer: %p, module %d, status %d, len %d, type %d\n", ptr, module, status, bklen, bktype);
      
         // print header
         
         int packetType = Alpha16Packet::PacketType(ptr, bklen);
         int packetVersion = Alpha16Packet::PacketVersion(ptr, bklen);
         
         if (0) {
            printf("Header:\n");
            printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
            printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
         }
         
         if (packetType == 1 && packetVersion == 1) {
            uint32_t ts = Alpha16Packet::PacketTimestamp(ptr, bklen);
            Alpha16Event *e = fEvb->FindEvent(0, ts);
            if (e)
               fEvb->AddBank(e, 0, ptr, bklen);
         } else {
            printf("unknown packet type %d, version %d\n", packetType, packetVersion);
            return;
         }
      } else {
         int imodule = 0;
         for (int i=0; i<NUM_CHAN_ALPHA16; i++) {
            char bname[5];
            sprintf(bname, "%2s%02d", fBankName.c_str(), i);

            status = me.FindBank(bname, &bklen, &bktype, &ptr);
            if (status == 1) {
               //printf("ALPHA16 bname %s, pointer: %p, status %d, len %d, type %d\n", bname, ptr, status, bklen, bktype);
      
               // print header
               
               int packetType = Alpha16Packet::PacketType(ptr, bklen);
               int packetVersion = Alpha16Packet::PacketVersion(ptr, bklen);
               
               if (0) {
                  printf("Header:\n");
                  printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
                  printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
               }
               
               if (packetType == 1 && packetVersion == 1) {
                  uint32_t ts = Alpha16Packet::PacketTimestamp(ptr, bklen);
                  Alpha16Event *e = fEvb->FindEvent(imodule, ts);
                  if (e)
                     fEvb->AddBank(e, imodule, ptr, bklen);
               } else {
                  printf("unknown packet type %d, version %d\n", packetType, packetVersion);
                  return;
               }
            }
         }
      }

    Alpha16Event* e = fEvb->GetNextEvent();

    if (!e)
      return;

    if (!e->complete)
      return;
    
    printf("event: "); e->Print();

    if (fLastEvent) {
      delete fLastEvent;
      fLastEvent = NULL;
    }
    
    fLastEvent = e;

    bool force_plot = false;

    //printf("next event: "); e->Print();

    // analyze the event

    if (fPlotWire)
       fPlotWire->Fill(e);

    if (fFileWire)
       fFileWire->Fill(e);
   };

  void PlotCanvas(TDataContainer& dataContainer, TRootEmbeddedCanvas *embedCanvas)
  {
     printf("Alpha16Canvas::PlotCanvas() for bank %s\n", fBankName.c_str());

    if (!fLastEvent)
      return;

    printf("plotting:   "); fLastEvent->Print();

    if (fPlotA16)
       fPlotA16->Draw(fLastEvent);

    if (fPlotWire) {
       fPlotWire->Draw(fLastEvent);
       fPlotWire->Update();
    }
  }
};
#endif

class AlphaTpcCanvas: public TCanvasHandleBase {
public:
   
   int fRunNo;
   
   Alpha16EVB* fEvb;
   
   PlotA16* fPlotA16[SHOW_ALPHA16];
   
   Alpha16Event* fLastEvent;

   PlotHistograms* fH;
   PlotHistogramsPads* fP;
   
   AlphaTpcCanvas() // ctor
      : TCanvasHandleBase("ALPHA TPC")
   {
      fH = new PlotHistograms(NULL);
      fP = new PlotHistogramsPads(NULL);

      fEvb = new Alpha16EVB(NUM_CHAN_ALPHA16, 512);

      for (int i=0; i<SHOW_ALPHA16; i++) {
         char title[256];
         sprintf(title, "ADC %2d", i+1);
         TCanvas *c = new TCanvas(title, title, 900, 400);
         fPlotA16[i] = new PlotA16(c, i*NUM_CHAN_ALPHA16);
      }
         
      fLastEvent = NULL;
   };
   
   ~AlphaTpcCanvas() // dtor
   {
      printf("AlphaTpcCanvas dtor!\n");
   };
   
   void BeginRun(int transition,int run,int time)
   {
      printf("AlphaTpcCanvas::BeginRun()\n");
      
      fRunNo = run;
      
      fEvb->Reset();
   };
   
   void EndRun(int transition,int run,int time)
   {
   };
   
   void ResetCanvasHistograms()
   {
      printf("AlphaTpcCanvas::ResetCanvasHistograms()\n");
   };
   
   void UpdateCanvasHistograms(TDataContainer& dataContainer)
   {
      //printf("AlphaTpcCanvas::UpdateCanvasHistograms()\n");

      Alpha16Event* e = UnpackAlpha16Event(fEvb, &dataContainer.GetMidasEvent());

#if 0      
      const TMidasEvent& me = dataContainer.GetMidasEvent();

      Alpha16Event* e = fEvb->NewEvent();
      
      for (int imodule = 0; imodule < SHOW_ALPHA16; imodule++) {
         for (int i=0; i<NUM_CHAN_ALPHA16; i++) {
            void *ptr;
            int bklen, bktype;
            int status;

            char bname[5];
            sprintf(bname, "A%01d%02d", 1+imodule, i);

            status = me.FindBank(bname, &bklen, &bktype, &ptr);
            if (status == 1) {
               //printf("ALPHA16 bname %s, pointer: %p, status %d, len %d, type %d\n", bname, ptr, status, bklen, bktype);
      
               // print header
               
               int packetType = Alpha16Packet::PacketType(ptr, bklen);
               int packetVersion = Alpha16Packet::PacketVersion(ptr, bklen);
               
               if (0) {
                  printf("Header:\n");
                  printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
                  printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
               }
               
               if (packetType == 1 && packetVersion == 1) {
                  fEvb->AddBank(e, imodule, ptr, bklen);
               } else {
                  printf("unknown packet type %d, version %d\n", packetType, packetVersion);
               }
            }
         }
      }

      fEvb->CheckEvent(e);
#endif
      printf("new event: "); e->Print();
      
      if (e->error)
         return;
      
      if (!e->complete)
         return;
      
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
                     
                     if (ple > 150 && ple < 170) {
                        fP->fHocc1->Fill(i);
                        fP->fHph1->Fill(pph);
                     }
                     
                     if (ple > 170 && ple < 580) {
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

            if (ph > 120) {
               fH->fHph->Fill(ph);

               int le = led(w, b, -1.0, ph/2.0);

               fH->fHle->Fill(le);
               fH->fHlex->Fill(le);
               fH->fHocc->Fill(i);

               if (le > 150 && le < 170) {
                  fH->fHocc1->Fill(i);
                  fH->fHph1->Fill(ph);
               }

               if (le > 170 && le < 580) {
                  fH->fHocc2->Fill(i);
                  fH->fHph2->Fill(ph);
               }

               fH->fHph3->Fill(le, ph);

               if (ph > 4000) {
                  nhits++;
                  //printf("samples %d %d, ", e->waveform[i].size(), w->nsamples);
                  printf("chan %4d: baseline %8.1f, rms %4.1f, range %8.1f %6.1f, pulse %6.1f, le %4d\n", i, b, brms, wmin, wmax, ph, le);
               }
            }

            delete w;
         }

      printf("Found %d hits\n", nhits);

      //if (nhits >= 2)
      //force_plot = true;

      TInterestingEventManager::instance()->Reset();
      if (force_plot) {
         printf("INTERESTING!\n");
         TInterestingEventManager::instance()->SetInteresting();
      }
   };

   void PlotCanvas(TDataContainer& dataContainer, TRootEmbeddedCanvas *embedCanvas)
   {
      printf("AlphaTpcCanvas::PlotCanvas()\n");
      
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
};

int doWrite = 0;

class MyTestLoop: public TRootanaDisplay { 

public:
	
   // An analysis manager.  Define and fill histograms in 
   // analysis manager.
   //TAnaManager *anaManager;
   
   MyTestLoop()
   {
      SetOutputFilename("example_output");
      DisableRootOutput(false);
      // Number of events to skip before plotting one.
      //SetNumberSkipEvent(10);
      // Choose to use functionality to update after X seconds
      SetOnlineUpdatingBasedSeconds();
      // Uncomment this to enable the 'interesting event' functionality.
      //iem_t::instance()->Enable();
   }
   
   void AddAllCanvases()
   {
      // Set up tabbed canvases

      AddSingleCanvas(new AlphaTpcCanvas());

#if 0      
      PlotA16* plotA16 = new PlotA16(NULL, 0);
      PlotWire* plotWire = new PlotWire(NULL);
      //AddSingleCanvas(new Alpha16Canvas("WIRE", plotA16, plotWire));
      AddSingleCanvas(new Alpha16Canvas("WIRE", "AZ", plotA16, plotWire, doWrite));
#endif

#if 0      
      TCanvas *c1 = new TCanvas("ALPHA16 ADC1", "ALPHA16 ADC1", 900, 650);
      PlotA16* plotA16_ADC1 = new PlotA16(c1, 0*NUM_CHAN_ALPHA16);
      AddSingleCanvas(new Alpha16Canvas("ADC1", "A1", plotA16_ADC1, NULL, 0));
      
      TCanvas *c2 = new TCanvas("ALPHA16 ADC2", "ALPHA16 ADC2", 900, 650);
      PlotA16* plotA16_ADC2 = new PlotA16(c2, 1*NUM_CHAN_ALPHA16);
      AddSingleCanvas(new Alpha16Canvas("ADC2", "A2", plotA16_ADC2, NULL, 0));
      
      TCanvas *c3 = new TCanvas("ALPHA16 ADC3", "ALPHA16 ADC3", 900, 650);
      PlotA16* plotA16_ADC3 = new PlotA16(c3, 2*NUM_CHAN_ALPHA16);
      AddSingleCanvas(new Alpha16Canvas("ADC3", "A3", plotA16_ADC3, NULL, 0));
      
      TCanvas *c4 = new TCanvas("ALPHA16 ADC4", "ALPHA16 ADC4", 900, 650);
      PlotA16* plotA16_ADC4 = new PlotA16(c4, 3*NUM_CHAN_ALPHA16);
      AddSingleCanvas(new Alpha16Canvas("ADC4", "A4", plotA16_ADC4, NULL, 0));
      
      TCanvas *c5 = new TCanvas("ALPHA16 ADC5", "ALPHA16 ADC5", 900, 650);
      PlotA16* plotA16_ADC5 = new PlotA16(c5, 4*NUM_CHAN_ALPHA16);
      AddSingleCanvas(new Alpha16Canvas("ADC5", "A5", plotA16_ADC5, NULL, 0));
#endif

      SetDisplayName("Example Display");
   };

   virtual ~MyTestLoop() {};

   void BeginRun(int transition,int run,int time)
   {
      std::cout << "User BOR method" << std::endl;
   }

   void EndRun(int transition,int run,int time)
   {
      std::cout << "User EOR method" << std::endl;
   }

   void ResetHistograms()
   {
   }

   void UpdateHistograms(TDataContainer& dataContainer)
   {
   }
   
   void PlotCanvas(TDataContainer& dataContainer)
   {
   }
}; 

int main(int argc, char *argv[])
{
   for (int i=0; i<argc; i++) {
      if (strcmp(argv[i], "--wire") == 0) {
         doWrite = atoi(argv[i+1]);
         printf("Will write %d events for 1-wire TPC prototype\n", doWrite);
         argc = i;
         break;
         //argv[i] = "";
         //argv[i+1] = "";
         i++;
      }
   }

   TInterestingEventManager::instance()->Enable();
   
   MyTestLoop::CreateSingleton<MyTestLoop>();  
   return MyTestLoop::Get().ExecuteLoop(argc, argv);
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
