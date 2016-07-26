#include <stdio.h>
#include <iostream>

#include "TRootanaDisplay.hxx"
#include "TH1D.h"

#include "TInterestingEventManager.hxx"

#include "Alpha16.h"

static int x = 0;

struct PlotA16
{
   TCanvas* fCanvas;
   TH1D* fH[16];

   PlotA16(TCanvas* c) // ctor
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

         if (!adc->chan[i])
            continue;
         
         if (!fH[i]) {
            char name[256];
            sprintf(name, "a16ch%02d_x%d", i, x++);
            fH[i] = new TH1D(name, name, adc->chan[i]->w->nsamples, 0, adc->chan[i]->w->nsamples);

            fH[i]->Draw();
            fH[i]->SetMinimum(-(1<<15));
            fH[i]->SetMaximum(1<<15);
            fH[i]->SetLineColor(color);
         }

         for (int s=0; s<adc->chan[i]->w->nsamples; s++)
            fH[i]->SetBinContent(s+1, adc->chan[i]->w->samples[s]);
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
               hh[i] = new TH1D(name, name, e->chan[ch]->w->nsamples, 0, e->chan[ch]->w->nsamples);

               if (color == 1)
                  hh[i]->Draw();
               else
                  hh[i]->Draw("same");
               hh[i]->SetMinimum(-(1<<15));
               hh[i]->SetMaximum(1<<15);
               hh[i]->SetLineColor(color);
            }
               
            for (int s=0; s<e->chan[ch]->w->nsamples; s++)
               hh[i]->SetBinContent(s+1, e->chan[ch]->w->samples[s]);
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

   void PrintWave(FILE* fp, int chan, const Waveform* w)
   {
      fprintf(fp, "%d", chan);
      for (int s=0; s<w->nsamples; s++)
         fprintf(fp, ",%d", (int)w->samples[s]);
      fprintf(fp, "\n");
   }

   void Fill(const Alpha16Event* e)
   {
      PrintWave(fp1, chan1, e->chan[chan1]->w);
      PrintWave(fp2, chan2, e->chan[chan2]->w);
      PrintWave(fp3, chan3, e->chan[chan3]->w);
      count++;

      if (max_write)
         if (count > max_write) {
            Close();
            exit(1);
         }
   }

};

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
      fEvb = new Alpha16EVB;
      fPlotA16 = plotA16;
      fPlotWire = plotWire;
      fLastEvent = NULL;
      fFileWire = NULL;
      
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
      int bklen,bktype;
      int status;
      
      status = me.FindBank(fBankName.c_str(), &bklen, &bktype, &ptr);
      if (status == 1) {
         //printf("ALPHA16 pointer: %p, module %d, status %d, len %d, type %d\n", ptr, module, status, bklen, bktype);
      
         // print header
         
         int packetType = getUint8(ptr, 0);
         int packetVersion = getUint8(ptr, 1);
         
         if (0) {
            printf("Header:\n");
            printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
            printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
         }
         
         if (packetType == 1 && packetVersion == 1) {
            Alpha16Packet* p = new Alpha16Packet(ptr, bklen*4);
            //p->Print();
            //printf("trig %d, ts 0x%08x, chan %d\n", p->acceptedTrigger, p->eventTimestamp, p->channelId);
            fEvb->AddPacket(p);
         } else {
            printf("unknown packet type %d, version %d\n", packetType, packetVersion);
            return;
         }
      } else {
         for (int i=0; i<16; i++) {
            char bname[5];
            sprintf(bname, "%2s%02d", fBankName.c_str(), i);

            status = me.FindBank(bname, &bklen, &bktype, &ptr);
            if (status == 1) {
               //printf("ALPHA16 bname %s, pointer: %p, status %d, len %d, type %d\n", bname, ptr, status, bklen, bktype);
      
               // print header
               
               int packetType = getUint8(ptr, 0);
               int packetVersion = getUint8(ptr, 1);
               
               if (0) {
                  printf("Header:\n");
                  printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
                  printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
               }
               
               if (packetType == 1 && packetVersion == 1) {
                  Alpha16Packet* p = new Alpha16Packet(ptr, bklen*4);
                  //p->Print();
                  //printf("trig %d, ts 0x%08x, chan %d\n", p->acceptedTrigger, p->eventTimestamp, p->channelId);
                  fEvb->AddPacket(p);
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

int doWrite = 0;

class MyTestLoop: public TRootanaDisplay { 

public:
	
  // An analysis manager.  Define and fill histograms in 
  // analysis manager.
  //TAnaManager *anaManager;

  MyTestLoop() {
    SetOutputFilename("example_output");
    DisableRootOutput(false);
    // Number of events to skip before plotting one.
    //SetNumberSkipEvent(10);
    // Choose to use functionality to update after X seconds
    SetOnlineUpdatingBasedSeconds();
    // Uncomment this to enable the 'interesting event' functionality.
    //iem_t::instance()->Enable();
  }

  void AddAllCanvases(){

    // Set up tabbed canvases

    PlotA16* plotA16 = new PlotA16(NULL);
    PlotWire* plotWire = new PlotWire(NULL);
    //AddSingleCanvas(new Alpha16Canvas("WIRE", plotA16, plotWire));
    AddSingleCanvas(new Alpha16Canvas("WIRE", "AZ", plotA16, plotWire, doWrite));

    TCanvas *c1 = new TCanvas("ALPHA16 ADC1", "ALPHA16 ADC1", 900, 650);
    PlotA16* plotA16_ADC1 = new PlotA16(c1);
    AddSingleCanvas(new Alpha16Canvas("ADC1", "A1", plotA16_ADC1, NULL, 0));

    TCanvas *c2 = new TCanvas("ALPHA16 ADC2", "ALPHA16 ADC2", 900, 650);
    PlotA16* plotA16_ADC2 = new PlotA16(c2);
    AddSingleCanvas(new Alpha16Canvas("ADC2", "A2", plotA16_ADC2, NULL, 0));

    TCanvas *c3 = new TCanvas("ALPHA16 ADC3", "ALPHA16 ADC3", 900, 650);
    PlotA16* plotA16_ADC3 = new PlotA16(c3);
    AddSingleCanvas(new Alpha16Canvas("ADC3", "A3", plotA16_ADC3, NULL, 0));

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
