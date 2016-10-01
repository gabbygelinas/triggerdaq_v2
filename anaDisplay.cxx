#include <stdio.h>
#include <iostream>

#include "TRootanaDisplay.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "TInterestingEventManager.hxx"

#include "Alpha16.h"
#include "Unpack.h"

#include "anaCommon.cxx"

class AlphaTpcCanvas: public TCanvasHandleBase {
public:

   AlphaTpcX* fATX;
   
   AlphaTpcCanvas() // ctor
      : TCanvasHandleBase("ALPHA TPC")
   {
      fATX = new AlphaTpcX();
   };
   
   ~AlphaTpcCanvas() // dtor
   {
      printf("AlphaTpcCanvas dtor!\n");
      if (fATX)
         delete fATX;
   };
   
   void BeginRun(int transition,int run,int time)
   {
      printf("AlphaTpcCanvas::BeginRun()\n");
      fATX->BeginRun(run);
   };
   
   void EndRun(int transition,int run,int time)
   {
      fATX->EndRun();
   };
   
   void ResetCanvasHistograms()
   {
      printf("AlphaTpcCanvas::ResetCanvasHistograms()\n");
   };
   
   void UpdateCanvasHistograms(TDataContainer& dataContainer)
   {
      //printf("AlphaTpcCanvas::UpdateCanvasHistograms()\n");

      Alpha16Event* e = UnpackAlpha16Event(fATX->fEvb, &dataContainer.GetMidasEvent());

      int force_plot = fATX->Event(e);

      TInterestingEventManager::instance()->Reset();
      if (force_plot) {
         printf("INTERESTING!\n");
         TInterestingEventManager::instance()->SetInteresting();
      }
   };

   void PlotCanvas(TDataContainer& dataContainer, TRootEmbeddedCanvas *embedCanvas)
   {
      printf("AlphaTpcCanvas::PlotCanvas()\n");
      fATX->Plot();
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
      SetSecondsBeforeUpdating(15);
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
