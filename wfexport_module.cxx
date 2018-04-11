//
// MIDAS analyzer example 2: ROOT analyzer
//
// K.Olchanski
//

#include <stdio.h>

#include "manalyzer.h"
#include "midasio.h"

#include <assert.h> // assert()

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "AgFlow.h"

#define ADC_MIN -33000
#define ADC_MAX  33000

// adjusted for 12-bit range 0xFFF = 4095
#define ADC_BINS 410
#define ADC_RANGE 4100
#define ADC_RANGE_PED 200

#define ADC_MIN_ADC -2048
#define ADC_OVERFLOW 4099

#define ADC_BINS_PULSER 100
#define ADC_RANGE_PULSER 300

#define ADC_RANGE_RMS 50

#define ADC_RMS_FPN_MIN 0
#define ADC_RMS_FPN_MAX 3.25

#define ADC_PULSER_TIME 450

#define NUM_SEQSCA (3*80+79)

#define NUM_TIME_BINS 512
#define MAX_TIME_BINS 512
#define MAX_TIME_NS 8200

class WfExportFlags
{
public:
   bool fExportWaveforms = false;

public:
   WfExportFlags() // ctor
   {
   }

   ~WfExportFlags() // dtor
   {
   }
};

class WfExportModule: public TARunObject
{
public:
   WfExportFlags* fFlags = NULL;

   //TDirectory* hdir_summary;
   //TDirectory* hdir_feam;
   //TDirectory* hdir_pads;
   //FeamHistograms fHF[MAX_FEAM];
   //std::vector<ChanHistograms*> fHC;

   bool fTrace = false;

   WfExportModule(TARunInfo* runinfo, WfExportFlags* f)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("WfExportModule::ctor!\n");

      fFlags = f;
   }

   ~WfExportModule()
   {
      if (fTrace)
         printf("WfExportModule::dtor!\n");
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("WfExportModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("ResumeRun, run %d\n", runinfo->fRunNo);
   }

   void ExportPads(TARunInfo* runinfo, const FeamEvent* pads)
   {
      if (1) {
         printf("ExportPads: ");
         pads->Print();
         printf("\n");
      }

      if (pads->error)
         return;

      for (unsigned ihit=0; ihit < pads->hits.size(); ihit++) {
         const FeamChannel* c = pads->hits[ihit];
         if (!c)
            continue;

         int ipwb = c->imodule;
         int pwb_column = c->pwb_column;
         int pwb_ring = c->pwb_ring;
         int isca = c->sca;
         int ichan = c->sca_readout;
         int scachan = c->sca_chan;
         int col = c->pad_col;
         int row = c->pad_row;

         //printf("hit %d: ipos %2d, isca %d, sca_readout %2d, sca_chan %2d, tpc col %2d, row %2d\n", ihit, ipos, isca, ichan, scachan, col, row);

         int first_bin = c->first_bin;
         unsigned nbins = c->adc_samples.size();

         assert(isca >= 0 && isca < 4);
         assert(ichan >=0 && ichan < 80);

         int seqsca = isca*80 + ichan;

         bool scachan_is_pad = (scachan > 0);
         bool scachan_is_fpn = (scachan >= -4) && (scachan <= -1);

         char xdir[256];
         char xhdr[256];
         char xname[256];
         char xtitle[256];

         sprintf(xdir, "pwb%02d_c%dr%d", ipwb, pwb_column, pwb_ring);
         sprintf(xhdr, "pwb %02d, col %d, ring %d", ipwb, pwb_column, pwb_ring);
         
         if (scachan_is_pad) {
            sprintf(xname, "%s_%03d_sca%d_ri%02d_scachan%02d_col%02d_row%02d", xdir, seqsca, isca, ichan, scachan, col, row);
            sprintf(xtitle, "%s, sca %d, readout index %d, sca chan %d, col %d, row %d", xhdr, isca, ichan, scachan, col, row);
         } else if (scachan_is_fpn) {
            sprintf(xname, "%s_%03d_sca%d_ri%02d_fpn%d", xdir, seqsca, isca, ichan, -scachan);
            sprintf(xtitle, "%s, sca %d, readout index %d, fpn %d", xhdr, isca, ichan, -scachan);
         } else {
            sprintf(xname, "%s_%03d_sca%d_ri%02d", xdir, seqsca, isca, ichan);
            sprintf(xtitle, "%s, sca %d, readout index %d", xhdr, isca, ichan);
         }
         
         if (fFlags->fExportWaveforms) {
            TDirectory* dir = runinfo->fRoot->fgDir;
            dir->cd();

            const char* xdir1 = "pwb_waveforms";
            
            std::string dirname = "";
            dirname += xdir1;
            dirname += "/";
            dirname += xdir;

            if (!dir->cd(dirname.c_str())) {
               dir = dir->mkdir(xdir1);
               dir->cd();
               dir = dir->mkdir(xdir);
               dir->cd();
            }
            
            dir = dir->CurrentDirectory();
            
            std::string wname = std::string(xname) + "_wf";
            std::string wtitle = std::string(xtitle) + " current waveform";
            
            TH1D* hwf = (TH1D*)dir->FindObject(wname.c_str());
            if (!hwf) {
               //printf("%s: ", wname.c_str()); dir->pwd();
               hwf = new TH1D(wname.c_str(), wtitle.c_str(), nbins, 0, nbins);
            }
            
            //printf("%s: first bin %d: ", wname.c_str(), first_bin);
            for (unsigned i=first_bin; i<nbins; i++) {
               //printf(" %4d", c->adc_samples[i]);
               hwf->SetBinContent(i+1, c->adc_samples[i]);
            }
            //printf("\n");
         }
      }
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      const AgEventFlow* ef = flow->Find<AgEventFlow>();
      //AgPadHitsFlow* padhitsf = flow->Find<AgPadHitsFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      const AgEvent* e = ef->fEvent;
      const FeamEvent* pads = e->feam;

      if (pads) {
         ExportPads(runinfo, pads);
      }

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class WfExportModuleFactory: public TAFactory
{
public:
   WfExportFlags fFlags;
   
public:
   void Init(const std::vector<std::string> &args)
   {
      printf("WfExportModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--wfexport")
            fFlags.fExportWaveforms = true;
      }
   }

   void Finish()
   {
      printf("WfExportModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("WfExportModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new WfExportModule(runinfo, &fFlags);
   }
};

static TARegister tar(new WfExportModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
