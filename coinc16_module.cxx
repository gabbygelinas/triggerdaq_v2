// Gabby Gelinas
// Module for producing coincidence plots with 4 paddles (8 scintillators)
#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "manalyzer.h"
#include "midasio.h"

#include "dltdc.h"

#include "DlFlow.h"

#include "ncfm.h"

#include <deque>

#include <TStyle.h>

#ifdef HAVE_ROOT
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#endif


class DlTdcFlags
{
public:
   bool fEnabled = false;
   bool fTriggered = false;
   bool fDebug = false;
   bool fPrint = false;
   bool fTWC = false;
};



class DlTdcCoincMap16
{
public:
   int fChanA = 32;
   int fChanB = 33;
   int fChanT = 34;

   std::vector<int> fMap;

   std::vector<int> fPair1;
   std::vector<int> fPair2;

public:
   bool Init(int runno);
};

bool DlTdcCoincMap16::Init(int runno)
{
   printf("DlTdcMap16 for run %d!\n", runno);

   fMap.resize(32+1);

   // TDC connector P1, A-side
   
   fMap[1]  =  0;
   fMap[2]  =  1;
   fMap[3]  = 10;
   fMap[4]  = 11;
   fMap[5]  =  2;
   fMap[6]  =  3;
   fMap[7]  =  8;
   fMap[8]  =  9;
   
   fMap[9]  = 15;
   fMap[10] = 14;
   fMap[11] =  7;
   fMap[12] =  6;
   fMap[13] = 13;
   fMap[14] = 12;
   fMap[15] =  5;
   fMap[16] =  4;

   // TDC connector P2, B-side
   
   fMap[17] =  0 + 16;
   fMap[18] =  1 + 16;
   fMap[19] = 10 + 16;
   fMap[20] = 11 + 16;
   fMap[21] =  2 + 16;
   fMap[22] =  3 + 16;
   fMap[23] =  8 + 16;
   fMap[24] =  9 + 16;
   
   fMap[25] = 15 + 16;
   fMap[26] = 14 + 16;
   fMap[27] =  7 + 16;
   fMap[28] =  6 + 16;
   fMap[29] = 13 + 16;
   fMap[30] = 12 + 16;
   fMap[31] =  5 + 16;
   fMap[32] =  4 + 16;
   
   fPair1.resize(16+1);
   fPair2.resize(16+1);

   // A-side
   
   fPair1[1] =  1; fPair2[1] =  9;
   fPair1[2] =  2; fPair2[2] = 10;
   fPair1[3] =  3; fPair2[3] = 11;
   fPair1[4] =  4; fPair2[4] = 12;
   fPair1[5] =  5; fPair2[5] = 13;
   fPair1[6] =  6; fPair2[6] = 14;
   fPair1[7] =  7; fPair2[7] = 15;
   fPair1[8] =  8; fPair2[8] = 16;

   // B-side

   fPair1[8+1] =  1+16; fPair2[8+1] =  9+16;
   fPair1[8+2] =  2+16; fPair2[8+2] = 10+16;
   fPair1[8+3] =  3+16; fPair2[8+3] = 11+16;
   fPair1[8+4] =  4+16; fPair2[8+4] = 12+16;
   fPair1[8+5] =  5+16; fPair2[8+5] = 13+16;
   fPair1[8+6] =  6+16; fPair2[8+6] = 14+16;
   fPair1[8+7] =  7+16; fPair2[8+7] = 15+16;
   fPair1[8+8] =  8+16; fPair2[8+8] = 16+16;

   return true;
};

class DlTdc16CoincModule: public TARunObject
{ 
public:
   DlTdcFlags* fFlags = NULL;
   //printf("into DlTdc16CoincModule class"); //testGG

   //#ifdef HAVE_ROOT
   // Declare the core file names for each plot (names as they appear on the side bar of the root browser)
   /* Declare all vectors of histograms. They follow the convension of <(int_a, int_b), *TH#D>.
   int_a and int_b run from 1 to 8 (inclusive) and denote the pair number of the scintillator 
   of interest. The pair number also gives the data channel number of one end of the scintillator.
   The other end of a scintillator is given by pair# + 8.

   Plots are grouped by the type of plot and historgrams for all possible coincidences for that 
   plot are contained within that vector.

   Notation:
   - The follow the convension plotType_pairOfInterest_pairsInCoincidence with "pa" denoting pair a.
   - cha denotes channel a (one end of pair a).
   - cut denotes a width cut
   - twc denotes time walk corrected data
   */

   // Make a vector of pairs (the lower numbered channel that we use to identify the scintillator) where each pair is one coincidence you want to examine. 
   //We do this so you don't get an overwhelming number of coincidence plots
   //std::vector<std::pair<int, int>> possibleCoinc = {std::make_pair(1, 7),std::make_pair(1, 8),std::make_pair(2, 7),std::make_pair(2, 8)};
   //printf("defined possibleCoinc"); //testGG
   std::vector<std::pair<int, int>> possibleCoinc = {std::make_pair(3, 19),std::make_pair(4, 19),std::make_pair(3, 20),std::make_pair(4, 20)};
   


   // Time difference along a scintillator
   std::map<std::pair<int,int>,TH1D*> timeDiff_pa_papb_cut; 
   std::map<std::pair<int,int>,TH1D*> timeDiff_pa_papb_cut_twc;
   std::map<std::pair<int,int>,TH1D*> timeDiff_pa_papb;

   // Time of flight
   std::map<std::pair<int,int>,TH1D*> tof_papb;
   std::map<std::pair<int,int>,TH1D*> tof_papb_cut;
   std::map<std::pair<int,int>,TH1D*> tof_papb_cut_twc;

   // Absolute time (time at one end minus the average time at the ends of the bar used for coincidence)
   std::map<std::pair<int,int>,TH1D*> absTime_chapb; 
   std::map<std::pair<int,int>,TH1D*> absTime_chapb_cut_twc;

   // 2D time difference
   std::map<std::pair<int,int>,TH2D*> timeDiff_pa_pb_papb;
   std::map<std::pair<int,int>,TH2D*> timeDiff_pa_pb_papb_cut;
   std::map<std::pair<int,int>,TH2D*> timeDiff_pa_pb_papb_cut_twc;

   // Widths with coincidence condition
   std::map<std::pair<int,int>,TH2D*> widths_pa_papb;

   // TOF vs width
   std::map<std::pair<int,int>,TH2D*> tofwidth_cha_papb;
   std::map<std::pair<int,int>,TH2D*> tofwidth_cha_papb_cut;
   std::map<std::pair<int,int>,TH2D*> tofwidth_cha_papb_cut_twc;

   // Absolute time vs width
   std::map<std::pair<int,int>,TH2D*> absTimeWidth_cha_papb;
   std::map<std::pair<int,int>,TH2D*> absTimeWidth_cha_papb_twc;

   //printf("initialized maps for storing plots"); //testGG

   char name_timeDiff_pa_papb[256];
   char name_timeDiff_pa_papb_cut[256];
   char name_timeDiff_pa_papb_cut_twc[256];

   char name_tof_papb[256];
   char name_tof_papb_cut[256];
   char name_tof_papb_cut_twc[256];

   char name_absTime_chapb[256];
   char name_absTime_chapb_cut_twc[256];

   char name_timeDiff_pa_pb_papb[256];
   char name_timeDiff_pa_pb_papb_cut[256];
   char name_timeDiff_pa_pb_papb_cut_twc[256];

   char name_widths_pa_papb[256];
   char name_widths_paHigh_papb[256];

   char name_tofwidth_cha_papb[256];
   char name_tofwidth_cha_papb_cut[256];
   char name_tofwidth_cha_papb_cut_twc[256];

   char name_absTime_cha_papb[256];
   char name_absTime_cha_papb_cut_twc[256];

   char name_absTimeWidth_cha_papb[256];
   char name_absTimeWidth_cha_papb_twc[256];


   // Declare titles for each plot
   char title_timeDiff_pa_papb[256];
   char title_timeDiff_pa_papb_cut[256];
   char title_timeDiff_pa_papb_cut_twc[256];

   char title_tof_papb[256];
   char title_tof_papb_cut[256];
   char title_tof_papb_cut_twc[256];

   char title_absTime_chapb[256];
   char title_absTime_chapb_cut_twc[256];

   char title_timeDiff_pa_pb_papb[256];
   char title_timeDiff_pa_pb_papb_cut[256];
   char title_timeDiff_pa_pb_papb_cut_twc[256];

   char title_widths_pa_papb[256];
   char title_widths_paHigh_papb[256];

   char title_tofwidth_cha_papb[256];
   char title_tofwidth_cha_papb_cut[256];
   char title_tofwidth_cha_papb_cut_twc[256];

   char title_absTimeWidth_cha_papb[256];
   char title_absTimeWidth_cha_papb_twc[256];

   double fPrevEventTimeSec = 0;

   //printf("initalized titles and names"); //testGG

   //testing
   double timeDiff_pair;
   double timeDiff_pa;
   double timeDiff_pb;

   //DlTdcMap8 *fMap16 = NULL;
   DlTdcCoincMap16 *fMap16 = NULL;

   bool fTrace = false;

   Ncfm* fCfm = NULL;

   double width_cut = 10.0; // Width to cut below in ns

   double ww_twc = 9.0;

   int numScints = 8; // total number of scintillators hooked up per arm
   std::vector<int> pairIdentifiers={1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24}; //integers used to identify the pair numbers. 1-8 correspond to e1-8, 17-24 correspond to p1-8.
   // The pair identifier +8 gives the opposite end of the scintillator of interest

   // Store the channel reading at the ends of each scintillator (channel pair) together in a map, identified by the pair number
   std::map<int, std::pair<double, double>> scintPairs; // <lower numbered pair, higher numbered pair>

   //printf("entering Dltdc16CoincModule \n"); //testGG
   DlTdc16CoincModule(TARunInfo* runinfo, DlTdcFlags* flags)
      : TARunObject(runinfo)
   {
      //printf("entering Dltdc16CoincModule\n"); //testGG
      if (fTrace)
         printf("DlTdc8Module::ctor!\n");

      printf("passed fTrace\n"); //testGG
      fModuleName = "dltdc16_module";
      fFlags   = flags;

      fCfm = new Ncfm("dlcfmdb");

      fMap16 = new DlTdcCoincMap16();
   }

   ~DlTdc16CoincModule()
   {
      if (fTrace)
         printf("DlTdc8Module::dtor!\n");

      if (fMap16) {
         delete fMap16;
         fMap16 = NULL;
      }

      if (fCfm) {
         delete fCfm;
         fCfm = NULL;
      }
   }

   void BeginRun(TARunInfo* runinfo)
   {
      printf("MODULE - coinc, begin run \n");
      if (fTrace)
         printf("DlTdcModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());

      if (!fFlags->fEnabled)
         return;

      bool conf_ok = fMap16->Init(runinfo->fRunNo);
      if (!conf_ok) {
         printf("Cannot load TDC map for run %d\n", runinfo->fRunNo);
         exit(123);
      }

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("coinc16");
      dir->cd(); // select correct ROOT directory

      //printf("about to create histograms that we can put stuff in later \n");
      std::vector<int> scintsExamined;
      // Initialize histograms for plots with times along scintillators, widths, and single channels
      for (int coincNum=0; coincNum<possibleCoinc.size(); coincNum++){
         //printf("in the loop \n"); // testingGG
         // Make the plot names and titles for each plot for the coincidence pairing of interest

         // Make a vector of each entry in this coincidence pairing so I can put it in a for loop 
         int coincPair1=possibleCoinc[coincNum].first; // Pair as in "pair" meaning scintillator. So what scintillator
         int coincPair2=possibleCoinc[coincNum].second;
         std::pair<int,int> thisCoincPairs (coincPair1,coincPair2); // the two scintillators we are looking at right now
         std::pair<int,int> thisCoincPairs_higherFirst (coincPair2,coincPair1);

         //printf("got the pairs \n"); // testingGG

         std::vector<int> thisCoincScint={coincPair1, coincPair2}; //Tracks which scintillators are in this coincidence
         //printf("tracking scintillators \n"); // testingGG

         // Doesn't specify a single scintillator or channel of interest so this can happen outside of the second for loop
         sprintf(name_tof_papb, "coinc_tof_cpair_%02d_%02d", coincPair1,coincPair2);
         sprintf(name_tof_papb_cut, "coinc_tof_cpair_%02d_%02d_cut", coincPair1,coincPair2);
         sprintf(name_tof_papb_cut_twc, "coinc_tof_cpair_%02d_%02d_cut_twc", coincPair1,coincPair2);

         sprintf(title_tof_papb, "TOF (ns) coinc with pairs %02d %02d", coincPair1,coincPair2);
         sprintf(title_tof_papb_cut, "TOF (ns) coinc with pairs %02d %02d cut", coincPair1,coincPair2);
         sprintf(title_tof_papb_cut_twc, "TOF (ns) coinc with pairs %02d %02d cut twc", coincPair1,coincPair2);

         //printf("made some names and titles \n"); // testingGG

         //printf("this Coinc Scint size %d \n", thisCoincScint.size()); // testingGG

         sprintf(name_timeDiff_pa_pb_papb, "coinc_tpair%d_%d_ns_cpair_%02d_%02d", coincPair1,coincPair2, coincPair1,coincPair2);
         sprintf(name_timeDiff_pa_pb_papb_cut, "coinc_tpair%d_%d_ns_cpair_%02d_%02d_cut", coincPair1,coincPair2, coincPair1,coincPair2);
         sprintf(name_timeDiff_pa_pb_papb_cut_twc, "coinc_tpair%d_%d_ns_cpair_%02d_%02d_cut_twc", coincPair1,coincPair2, coincPair1,coincPair2);

         sprintf(title_timeDiff_pa_pb_papb, "Time difference along pair %d vs pair %d (ns) with coinc between pairs %02d %02d", coincPair1,coincPair2, coincPair1,coincPair2);
         sprintf(title_timeDiff_pa_pb_papb_cut, "Time difference along pair %d vs pair %d (ns) with coinc between pairs %02d %02d cut", coincPair1,coincPair2,coincPair1,coincPair2);
         sprintf(title_timeDiff_pa_pb_papb_cut_twc, "Time difference along pair %d vs pair %d (ns) with coinc between pairs %02d %02d cut twc", coincPair1,coincPair2, coincPair1,coincPair2);

         sprintf(name_timeDiff_pa_papb, "coinc_tpair%d_ns_cpair_%02d_%02d", coincPair1, coincPair1,coincPair2);
         sprintf(name_timeDiff_pa_papb_cut, "coinc_tpair%d_ns_cpair_%02d_%02d_cut", coincPair1, coincPair1,coincPair2);
         sprintf(name_timeDiff_pa_papb_cut_twc, "coinc_tpair%d_ns_cpair_%02d_%02d_cut_twc", coincPair1, coincPair1,coincPair2);//

         sprintf(title_timeDiff_pa_papb, "Time difference along pair %d (ns) with coinc between pairs %02d %02d",  coincPair1, coincPair1,coincPair2);
         sprintf(title_timeDiff_pa_papb_cut, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut", coincPair1, coincPair1,coincPair2);
         sprintf(title_timeDiff_pa_papb_cut_twc, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut twc", coincPair1, coincPair1,coincPair2);

         //std::pair<int,int> scintPair_alongInterest1 (coincPair1,coincPair1+numScints);


         timeDiff_pa_papb[thisCoincPairs] = new TH1D(name_timeDiff_pa_papb, title_timeDiff_pa_papb, 200, -10, 10);
         //timeDiff_pa_papb[thisCoincPairs] = new TH1D(TString::s("coinc_tpair%d_%d_ns_cpair_%02d_%02d", coincPair1,coincPair2, coincPair1,coincPair2), title_timeDiff_pa_papb, 200, -10, 10);
         timeDiff_pa_papb_cut[thisCoincPairs] = new TH1D(name_timeDiff_pa_papb_cut, title_timeDiff_pa_papb_cut, 200, -10, 10);
         timeDiff_pa_papb_cut_twc[thisCoincPairs] = new TH1D(name_timeDiff_pa_papb_cut_twc, title_timeDiff_pa_papb_cut_twc, 200, -10, 10);

         sprintf(name_timeDiff_pa_papb, "coinc_tpair%d_ns_cpair_%02d_%02d", coincPair2, coincPair1,coincPair2);
         sprintf(name_timeDiff_pa_papb_cut, "coinc_tpair%d_ns_cpair_%02d_%02d_cut", coincPair2, coincPair1,coincPair2);
         sprintf(name_timeDiff_pa_papb_cut_twc, "coinc_tpair%d_ns_cpair_%02d_%02d_cut_twc", coincPair2, coincPair1,coincPair2);

         sprintf(title_timeDiff_pa_papb, "Time difference along pair %d (ns) with coinc between pairs %02d %02d", coincPair2, coincPair1,coincPair2);
         sprintf(title_timeDiff_pa_papb_cut, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut", coincPair2, coincPair1,coincPair2);
         sprintf(title_timeDiff_pa_papb_cut_twc, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut twc", coincPair2, coincPair1,coincPair2);

         //std::pair<int,int> scintPair_alongInterest2 (coincPair2,coincPair2+numScints);

         timeDiff_pa_papb[thisCoincPairs_higherFirst] = new TH1D(name_timeDiff_pa_papb, title_timeDiff_pa_papb, 200, -10, 10);
         timeDiff_pa_papb_cut[thisCoincPairs_higherFirst] = new TH1D(name_timeDiff_pa_papb_cut, title_timeDiff_pa_papb_cut, 200, -10, 10);
         timeDiff_pa_papb_cut_twc[thisCoincPairs_higherFirst] = new TH1D(name_timeDiff_pa_papb_cut_twc, title_timeDiff_pa_papb_cut_twc, 200, -10, 10);

         timeDiff_pa_pb_papb[thisCoincPairs] = new TH2D(name_timeDiff_pa_pb_papb, title_timeDiff_pa_pb_papb, 200, -10, 10, 200, -10, 10);
         timeDiff_pa_pb_papb_cut[thisCoincPairs] = new TH2D(name_timeDiff_pa_pb_papb_cut, title_timeDiff_pa_pb_papb_cut, 200, -10, 10, 200, -10, 10);
         timeDiff_pa_pb_papb_cut_twc[thisCoincPairs] = new TH2D(name_timeDiff_pa_pb_papb_cut_twc, title_timeDiff_pa_pb_papb_cut_twc, 200, -10, 10, 200, -10, 10);
            
         tof_papb[thisCoincPairs] = new TH1D(name_tof_papb, title_tof_papb, 200, -10, 10);
         tof_papb_cut[thisCoincPairs] = new TH1D(name_tof_papb_cut, title_tof_papb_cut, 200, -10, 10);
         tof_papb_cut_twc[thisCoincPairs] = new TH1D(name_tof_papb_cut_twc, title_tof_papb_cut_twc, 200, -10, 10);


         // Make the rest of the names and titles then apply them to histograms we will create
         for (int i=0; i<(thisCoincScint.size()); i++){
            //printf("in the second loop \n"); // testingGG
            int scintOfInterest = thisCoincScint[i];
            int otherScint;
            if (i==0){
               otherScint = thisCoincScint[1];
            }
            else {
               otherScint = thisCoincScint[0];
            }

            // Initialize a pair variable for storing the information of the scintillators we are examining for accessing our map element
            std::pair<int,int> scintPair (scintOfInterest,otherScint);

            //printf("making plot names \n"); // testingGG

            // Make plot names for plots that focus on the lower numbered channel of a scintillator
            sprintf(name_timeDiff_pa_papb, "coinc_tpair%d_ns_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_timeDiff_pa_papb_cut, "coinc_tpair%d_ns_cpair_%02d_%02d_cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_timeDiff_pa_papb_cut_twc, "coinc_tpair%d_ns_cpair_%02d_%02d_cut_twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_widths_pa_papb, "coinc_width%d_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_absTime_cha_papb, "coinc_absTime_ch%d_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_absTime_cha_papb_cut_twc, "coinc_absTime_ch%d_cpair_%02d_%02d_cut_twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_absTimeWidth_cha_papb, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d", scintOfInterest, scintOfInterest, coincPair1,coincPair2);
            sprintf(name_absTimeWidth_cha_papb_twc, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d_twc", scintOfInterest, scintOfInterest, coincPair1,coincPair2);

            //printf("making plot titles \n"); // testingGG

            // Make plot titles for plots that focus on the lower numbered channel of a scintillator
            sprintf(title_timeDiff_pa_papb, "Time difference along pair %d (ns) with coinc between pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_timeDiff_pa_papb_cut, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_timeDiff_pa_papb_cut_twc, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(title_widths_pa_papb, "Width on channel %d, coinc between pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);

            sprintf(title_tofwidth_cha_papb, "TOF vs Width on channel%d with coinc between pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut_twc, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(title_absTime_chapb, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d", scintOfInterest, otherScint, otherScint+numScints,coincPair1,coincPair2);
            sprintf(title_absTime_chapb_cut_twc, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d cut twc", scintOfInterest, otherScint, otherScint+numScints,coincPair1,coincPair2);

            sprintf(title_absTimeWidth_cha_papb, "TOF vs width%d coinc btwn pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_absTimeWidth_cha_papb_twc, "TOF vs width%d coinc btwn pairs %02d %02d twc", scintOfInterest, coincPair1,coincPair2);

            //std::cout << "Up here scintPair: " << scintPair.first << " " << scintPair.second << std::endl;

            widths_pa_papb[scintPair] = new TH2D(name_widths_pa_papb, title_widths_pa_papb, 200, -10, 10, 200, -10, 10);

            tofwidth_cha_papb[scintPair] = new TH2D(name_tofwidth_cha_papb, title_tofwidth_cha_papb, 200, -10, 10, 200, -10, 10);
            tofwidth_cha_papb_cut[scintPair] = new TH2D(name_tofwidth_cha_papb_cut, title_tofwidth_cha_papb_cut, 200, -10, 10, 200, -10, 10);
            tofwidth_cha_papb_cut_twc[scintPair] = new TH2D(name_tofwidth_cha_papb_cut_twc, title_tofwidth_cha_papb_cut_twc, 200, -10, 10, 200, -10, 10);

            absTime_chapb[scintPair] = new TH1D(name_absTime_cha_papb, title_absTime_chapb, 200, -10, 10);
            absTime_chapb_cut_twc[scintPair] = new TH1D(name_absTime_cha_papb_cut_twc, title_absTime_chapb_cut_twc, 200, -10, 10);

            absTimeWidth_cha_papb[scintPair] = new TH2D(name_absTimeWidth_cha_papb, title_absTimeWidth_cha_papb, 200, -10, 10, 200, -10, 10);
            absTimeWidth_cha_papb_twc[scintPair] = new TH2D(name_absTimeWidth_cha_papb_twc, title_absTimeWidth_cha_papb_twc, 200, -10, 10, 200, -10, 10);


            // Make plot names for plots that focus on the higher number channel on a scintillator
            sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", scintOfInterest+numScints, coincPair1,coincPair2); //new
            sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", scintOfInterest+numScints, coincPair1,coincPair2); //new
            sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", scintOfInterest+numScints, coincPair1,coincPair2);//new

            sprintf(name_absTime_cha_papb, "coinc_absTime_ch%d_cpair_%02d_%02d", scintOfInterest+numScints, coincPair1,coincPair2); //new
            sprintf(name_absTime_cha_papb_cut_twc, "coinc_absTime_ch%d_cpair_%02d_%02d_cut_twc", scintOfInterest+numScints, coincPair1,coincPair2); //new

            sprintf(name_absTimeWidth_cha_papb, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d", scintOfInterest+numScints, scintOfInterest+numScints, coincPair1,coincPair2); //new
            sprintf(name_absTimeWidth_cha_papb_twc, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d_twc", scintOfInterest+numScints, i+numScints, coincPair1,coincPair2); //new

            sprintf(name_widths_pa_papb, "coinc_width%d_cpair_%02d_%02d", scintOfInterest+numScints, coincPair1,coincPair2);


            // Make plot titles for plots that focus on the higher number channel on a scintillator
            sprintf(title_tofwidth_cha_papb, "TOF vs Width on channel%d with coinc between pairs %02d %02d", scintOfInterest+numScints, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut", scintOfInterest+numScints, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut_twc, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut twc", scintOfInterest+numScints, coincPair1,coincPair2);

            sprintf(title_absTime_chapb, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d", scintOfInterest+numScints, otherScint, otherScint+numScints,coincPair1,coincPair2);
            sprintf(title_absTime_chapb_cut_twc, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d cut twc", scintOfInterest+numScints, otherScint, otherScint+numScints,coincPair1,coincPair2);

            sprintf(title_absTimeWidth_cha_papb, "TOF vs width%d coinc btwn pairs %02d %02d", scintOfInterest+numScints, coincPair1,coincPair2);
            sprintf(title_absTimeWidth_cha_papb_twc, "TOF vs width%d coinc btwn pairs %02d %02d twc", scintOfInterest+numScints, coincPair1,coincPair2);

            sprintf(title_widths_pa_papb, "Width on channel %d, coinc between pairs %02d %02d", scintOfInterest+numScints, coincPair1,coincPair2);



            // Initalize histograms for plots that focus on the higher number channel on a scintillator
            std::pair<int,int> scintPair_interest8 (scintOfInterest+numScints,otherScint);
            //std::cout << "Up here scintPair_interest8: " << scintPair_interest8.first << " " << scintPair_interest8.second << std::endl;

            tofwidth_cha_papb[scintPair_interest8] = new TH2D(name_tofwidth_cha_papb, title_tofwidth_cha_papb, 200, -10, 10, 200, -10, 10);
            tofwidth_cha_papb_cut[scintPair_interest8] = new TH2D(name_tofwidth_cha_papb_cut, title_tofwidth_cha_papb_cut, 200, -10, 10, 200, -10, 10);
            tofwidth_cha_papb_cut_twc[scintPair_interest8] = new TH2D(name_tofwidth_cha_papb_cut_twc, title_tofwidth_cha_papb_cut_twc, 200, -10, 10, 200, -10, 10);

            absTime_chapb[scintPair_interest8] = new TH1D(name_absTime_cha_papb, title_absTime_chapb, 200, -10, 10);
            absTime_chapb_cut_twc[scintPair_interest8] = new TH1D(name_absTime_cha_papb_cut_twc, title_absTime_chapb_cut_twc, 200, -10, 10);

            absTimeWidth_cha_papb[scintPair_interest8] = new TH2D(name_absTimeWidth_cha_papb, title_absTimeWidth_cha_papb, 200, -10, 10, 200, -10, 10);
            absTimeWidth_cha_papb_twc[scintPair_interest8] = new TH2D(name_absTimeWidth_cha_papb_twc, title_absTimeWidth_cha_papb_twc, 200, -10, 10, 200, -10, 10);

            widths_pa_papb[scintPair_interest8] = new TH2D(name_widths_pa_papb, title_widths_pa_papb, 200, -10, 10, 200, -10, 10);

         }
      }

      for (int i=1; i<(numScints+1); i++){
         scintPairs.insert({i, std::make_pair(fMap16->fMap[i],fMap16->fMap[i+numScints])});
      }

   }


void PreEndRun(TARunInfo* runinfo)
   {  
      printf("MODULE - coinc, pre end run\n");
      if (fTrace)
         printf("DlTdcModule::PreEndRun, run %d\n", runinfo->fRunNo);
   }
   
   void EndRun(TARunInfo* runinfo)
   {  
      printf("MODULE - coinc, end run \n");
      if (fTrace)
         printf("DlTdcModule::EndRun, run %d\n", runinfo->fRunNo);

      if (!fFlags->fEnabled)
         return;
   }

void PauseRun(TARunInfo* runinfo)
   {
      printf("MODULE - coinc, pause run \n");
      if (fTrace)
         printf("DlTdcModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      printf("MODULE - coinc, resume run \n");
      if (fTrace)
         printf("DlTdcModule::ResumeRun, run %d\n", runinfo->fRunNo);

   }


// Define function to fill histograms
void do_quadcoinc(int scint_a, int scint_b, const DlTdcEvent& t, std::vector<double> ww_ns,TARunInfo* runinfo)
//void do_quadcoinc(int scint_a, int scint_b, const DlTdcEvent& t, std::vector<double> ww_ns)
   {
   //printf("entered do_quadcoinc \n"); //testgg
   //printf("MODULE - coinc, doing quad coinc GABBY\n");  
    /* Fill histograms for the appropriate quad coincidences
    scint_a - one scintillator pair number (1-8) you want to examine, int
    scint_b - the other scintillator pair number (1-8) you want to examine, int

    NOTE: a and b are different from A and B used to denote the two data cables
    */

    // Time difference along each scintillator (channel a+numScints) - (channel a)
    double t_scinta_ns = subtract_ns(t.GetCh(scintPairs[scint_a].second).fLe, t.GetCh(scintPairs[scint_a].first).fLe);
    double t_scintb_ns = subtract_ns(t.GetCh(scintPairs[scint_b].second).fLe, t.GetCh(scintPairs[scint_b].first).fLe);
   
    //printf("made t_scinta_ns and b \n"); //testgg

    if (t_scinta_ns > -9999 && t_scintb_ns > -9999){
         //printf("signals exist \n"); //testgg
 
         // Make a vector with int entries of scint_a, scint_b, and those two +numScints. We will use this list to fill the above calculation vectors.
         std::vector<int> channelNumbers;
         channelNumbers.push_back(scint_a);
         channelNumbers.push_back(scint_a+numScints);
         channelNumbers.push_back(scint_b);
         channelNumbers.push_back(scint_b+numScints);

         //printf("channel numbers saved \n"); //testgg
 
         // testing
         //printf("\n scint_a=%d \n", scint_a);
         //printf("Channel numbers 0, %d,num1=%d, num2=%d, number 3, %d", channelNumbers[0], channelNumbers[1],channelNumbers[2],channelNumbers[3]);

         for (int i=0; i<=3; ++i){
            for (int j=0; j<=3; ++j){
               //printf("i=%d, j=%d \n", i, j); //testgg
               // testing: Making it into this loop correctly
               int chan_i = channelNumbers[i];
               //int chan_i8 = channelNumbers[i] + 8;
               int chan_j = channelNumbers[j];
               int chan_j8 = channelNumbers[j] + 8;
               int chan_i8 = channelNumbers[i] + 8;

               // Pair of the channels we are examining for this coincidence
               std::pair<int,int> chanPair_ji (channelNumbers[j],channelNumbers[i]);
               std::pair<int,int> chanPair_ij (channelNumbers[i],channelNumbers[j]);
               //printf("chanPairs made \n"); //testgg
               //printf(" At the top: chanPair_ij=%d,%d \n",channelNumbers[i],channelNumbers[j]);

               std::pair<int,int> chanPair_i8j (channelNumbers[i]+numScints,channelNumbers[j]);
               std::pair<int,int> chanPair_j8i (channelNumbers[j]+numScints,channelNumbers[i]);
               //printf("chanPairs +numScints made \n"); //testgg
               //std::cout<< "chanPair_ij:" << chanPair_ij.first << " " << chanPair_ij.second<<std::endl;
               //std::cout<< "chanPair_ji:" <<chanPair_ji.first  << " " << chanPair_ji.second<<std::endl;
               //std::cout<< "chanPair_i8j:" <<chanPair_i8j.first  << " " << chanPair_i8j.second<<std::endl;


                // The graphs in this condition get filled
                if ((i != j) && ((chan_i<9)||(16<chan_i && chan_i<25)) && ((chan_j<9)||(16<chan_j && 16<chan_j && chan_j<25)) && (chan_i < chan_j) && (chan_i+8 != chan_j)){
                  //printf("inside loop for chanel numbers that are the different (two scintillator) \n"); //testgg
                  //printf("chan_i = %d, chan_j=%d \n", chan_i, chan_j); //testgg
                  // Doing subtractions between different scintillators. Can do non-physical and absoulte time here (do absolute by adding 8 to one)
                  // The last condition is to make sure that we follow Konstantin's convension of doing higher numbered channel minus lower
                  // Define channel 0 as the lower number channel and channel 1 and the higher numbered channel (channel 0) + 8, so the opposite channel on the same scintillator
                  // All of these are set up as a-b
                  //printf("Different scintillators \n");
                  double timeDiff_a0_b0 = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].first).fLe, t.GetCh(scintPairs[channelNumbers[j]].first).fLe);
                  double timeDiff_a0_b1 = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].first).fLe, t.GetCh(scintPairs[channelNumbers[j]].second).fLe);
                  double timeDiff_a1_b0 = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].second).fLe, t.GetCh(scintPairs[channelNumbers[j]].first).fLe);
                  double timeDiff_a1_b1 = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].second).fLe, t.GetCh(scintPairs[channelNumbers[j]].second).fLe);
                  //printf("Calculates first block \n");
                  //printf("different scint time diff calcs done \n"); //testgg

                  double tof_a0a1b0b1 = -(timeDiff_a0_b0 + timeDiff_a1_b1);
                  double tof_a0a1b0b1_twc = tof_a0a1b0b1 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[i]]) + ww_twc / sqrt(ww_ns[channelNumbers[i]+numScints]) - ww_twc / sqrt(ww_ns[channelNumbers[j]]) - ww_twc / sqrt(ww_ns[channelNumbers[j]+numScints])); 
                  //printf("Calculates second block \n");
                  //printf("different scint tof calcs done \n"); //testgg

                  double absTime_a0 = 0.5*(timeDiff_a0_b0 + timeDiff_a0_b1);
                  double absTime_a1 = 0.5*(timeDiff_a1_b0 + timeDiff_a1_b1);
                  double absTime_b0 = -0.5*(timeDiff_a0_b0 + timeDiff_a1_b0);
                  double absTime_b1 = -0.5*(timeDiff_a0_b1 + timeDiff_a1_b1);
                  //printf("Calculates third block \n");
                  //printf("different scint abs time calcs 1 done \n"); //testgg

                  double absTime_a0_twc = absTime_a0 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[j]]) + ww_twc / sqrt(ww_ns[channelNumbers[j]+numScints]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[i]]));
                  double absTime_a1_twc = absTime_a1 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[j]]) + ww_twc / sqrt(ww_ns[channelNumbers[j]+numScints]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[i]+numScints]));
                  double absTime_b0_twc = absTime_b0 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[i]]) + ww_twc / sqrt(ww_ns[channelNumbers[i]+numScints]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[j]]));
                  double absTime_b1_twc = absTime_b1 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[i]]) + ww_twc / sqrt(ww_ns[channelNumbers[i]+numScints]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[j]+numScints]));
                  //printf("Calculates fourth block \n");
                  //printf("different scint abs time 2 calcs done \n"); //testgg
                  //std::cout << "timeDiff_pa_pb_papb[chanPair_ij] " << timeDiff_pa_pb_papb[chanPair_ij] << std::endl;
                  //std::cout<< "chanPair_ij:" << chanPair_ij.first << " " << chanPair_ij.second<<std::endl;

                  timeDiff_pa = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].second).fLe, t.GetCh(scintPairs[channelNumbers[i]].first).fLe);
                  timeDiff_pb = subtract_ns(t.GetCh(scintPairs[channelNumbers[j]].second).fLe, t.GetCh(scintPairs[channelNumbers[j]].first).fLe);
                  //printf("DIFFERENT t.GetCh(scintPairs[channelNumbers[i]].second).fLe: %f \n", t.GetCh(scintPairs[channelNumbers[i]].second).fLe);
                  //printf("DIFFERENT timeDiff_pa: %.16f \n", timeDiff_pa);
                    
                  double timeDiff_pa_twc = timeDiff_pa + (ww_twc / sqrt(ww_ns[channelNumbers[j]]) - ww_twc / sqrt(ww_ns[channelNumbers[i]]));
                  double timeDiff_pb_twc = timeDiff_pa + (ww_twc / sqrt(ww_ns[channelNumbers[i]]) - ww_twc / sqrt(ww_ns[channelNumbers[j]])); //check
                  
                  //std::cout << "timeDiff_pa_papb[chanPair_ij] " << timeDiff_pa_papb[chanPair_ij] << std::endl;
                  //std::cout << "timeDiff_pa_papb[chanPair_ji] " << timeDiff_pa_papb[chanPair_ji] << std::endl;

                  timeDiff_pa_papb[chanPair_ij]->Fill(timeDiff_pa);
                  timeDiff_pa_papb[chanPair_ji]->Fill(timeDiff_pb);

                  timeDiff_pa_pb_papb[chanPair_ij]->Fill(timeDiff_pa,timeDiff_pb);
                  //printf("done timeDiff_pa_pb_papb \n");

                  // Fill histograms with no width cut
                  tof_papb[chanPair_ij]->Fill(tof_a0a1b0b1);

                  absTime_chapb[chanPair_ij]->Fill(absTime_a0);
                  absTime_chapb[chanPair_i8j]->Fill(absTime_a1);
                  absTime_chapb[chanPair_ji]->Fill(absTime_b0);
                  absTime_chapb[chanPair_j8i]->Fill(absTime_b1);

                  widths_pa_papb[chanPair_ij]->Fill(ww_ns[channelNumbers[i]],ww_ns[channelNumbers[j]]);

                  tofwidth_cha_papb[chanPair_ij]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]]);
                  tofwidth_cha_papb[chanPair_i8j]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]+numScints]);
                  tofwidth_cha_papb[chanPair_ji]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]]);
                  tofwidth_cha_papb[chanPair_j8i]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]+numScints]);
                  // end of filling histograms with no width cut



                  // Fill histograms with a width cut
                  if (ww_ns[channelNumbers[i]] > width_cut && ww_ns[channelNumbers[i]+numScints] > width_cut && ww_ns[channelNumbers[j]] > width_cut && ww_ns[channelNumbers[j]+numScints] > width_cut) {
                     //printf("Enters different scintillators width cut \n");

                     //std::cout << "chanPair_ij " << chanPair_ij.first << " " << chanPair_ij.second << std::endl;
                     //std::cout << "tof_a0a1b0b1 " << tof_a0a1b0b1 << std::endl;
                     //std::cout << "tof_a0a1b0b1_twc " << tof_a0a1b0b1_twc << std::endl;
                     //std::cout << "tof_papb_cut " << tof_papb_cut[chanPair_ij] << std::endl;
                     //std::cout << "tof_papb_cut_twc " << tof_papb_cut_twc[chanPair_ij] << std::endl;

                     timeDiff_pa_papb_cut[chanPair_ij]->Fill(timeDiff_pa);
                     timeDiff_pa_papb_cut_twc[chanPair_ij]->Fill(timeDiff_pa_twc);
                     timeDiff_pa_papb_cut[chanPair_ji]->Fill(timeDiff_pb);
                     timeDiff_pa_papb_cut_twc[chanPair_ji]->Fill(timeDiff_pb_twc);

               

                     tof_papb_cut[chanPair_ij]->Fill(tof_a0a1b0b1);
                     tof_papb_cut_twc[chanPair_ij]->Fill(tof_a0a1b0b1_twc);
                     //printf("First block width calculations done \n");

                     absTime_chapb_cut_twc[chanPair_ij]->Fill(absTime_a0_twc);
                     absTime_chapb_cut_twc[chanPair_i8j]->Fill(absTime_a1_twc);
                     absTime_chapb_cut_twc[chanPair_ji]->Fill(absTime_b0_twc);
                     absTime_chapb_cut_twc[chanPair_j8i]->Fill(absTime_b1_twc);
                     //printf("Second block width calculations done \n");

                     //std::cout << "timeDiff_pa_pb_papb_cut " << timeDiff_pa_pb_papb_cut[chanPair_ij] << std::endl;
                     timeDiff_pa_pb_papb_cut[chanPair_ij]->Fill(t_scinta_ns,t_scintb_ns);
                     timeDiff_pa_pb_papb_cut_twc[chanPair_ij]->Fill(timeDiff_pa_twc,timeDiff_pb_twc);
                     //printf("Third block width calculations done \n");

                     tofwidth_cha_papb_cut[chanPair_ij]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]]);
                     tofwidth_cha_papb_cut[chanPair_i8j]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]+numScints]);
                     tofwidth_cha_papb_cut[chanPair_ji]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]]);
                     tofwidth_cha_papb_cut[chanPair_j8i]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]+numScints]);

                     tofwidth_cha_papb_cut_twc[chanPair_ij]->Fill(tof_a0a1b0b1_twc, ww_ns[channelNumbers[i]]);
                     tofwidth_cha_papb_cut_twc[chanPair_i8j]->Fill(tof_a0a1b0b1_twc, ww_ns[channelNumbers[i]+numScints]);
                     tofwidth_cha_papb_cut_twc[chanPair_ji]->Fill(tof_a0a1b0b1_twc, ww_ns[channelNumbers[j]]);
                     tofwidth_cha_papb_cut_twc[chanPair_j8i]->Fill(tof_a0a1b0b1_twc, ww_ns[channelNumbers[j]+numScints]);

                     absTimeWidth_cha_papb_twc[chanPair_ij] -> Fill(absTime_a0_twc, ww_ns[channelNumbers[i]]);
                     absTimeWidth_cha_papb_twc[chanPair_i8j]-> Fill(absTime_a1_twc, ww_ns[channelNumbers[i]+numScints]);
                     absTimeWidth_cha_papb_twc[chanPair_ji] -> Fill(absTime_b0_twc, ww_ns[channelNumbers[j]]);
                     absTimeWidth_cha_papb_twc[chanPair_j8i]-> Fill(absTime_b1_twc, ww_ns[channelNumbers[j]+numScints]); 
                     //printf("Fourth block width calculations done \n");
                    }
                }
            }
        }
    }
}


// Call function for each possible quad coincidence in this function
void AnalyzeTdcEvent(const DlTdcEvent& t, TARunInfo* runinfo)
   {  
      printf("New AnalyzeTDCEvent \n"); //testgg
      if (fFlags->fDebug) {
         std::string s = "";
         for (size_t i=1; i<fMap16->fMap.size(); i++) {
            int tdc_ch = fMap16->fMap[i];
            if (t.HaveCh(tdc_ch))
               s += "H";
            else
               s += ".";
         }

         printf("EVENT %s, ABT %d%d%d\n", s.c_str(), t.HaveCh(fMap16->fChanA), t.HaveCh(fMap16->fChanB), t.HaveCh(fMap16->fChanT));
         printf("fMap", fMap16->fMap[1]);
      }

      ///////// check for triggered event ///////////
         
      if (fFlags->fTriggered && !t.HaveCh(fMap16->fChanT)) {
         return;
      }

///////// SET WIDTH CUTOFFS (NS)  /////////
      
      //double ww_cut_ns = 10;

      ///////// COMPUTE WIDTH AND PULSE HEIGHT ///////////

      std::vector<double> ww_ns(fMap16->fMap.size() + 1);
      //printf("made width vector \n"); //testgg

      for (size_t i=1; i<fMap16->fMap.size(); i++) {
         int tdc_ch = fMap16->fMap[i];

         double w_ns = -9999;

         if (t.HaveCh(tdc_ch)) {
            w_ns = t.GetCh(tdc_ch).fWidthNs;

            if (w_ns < 0.01) {
               printf("WWW: BAD WIDTH chan %02d, tdc_ch %02d: %f ns!\n", (int)i, (int)tdc_ch, w_ns);
               w_ns = -9999;
            }
         }

         ww_ns[i] = w_ns;
      }
      //printf("filled width vector \n"); //testgg

      // SET UP TIME WALK CORRECT VALUES //
      //double ww_twc = 9.0;

   // ALL MY COINCIDENCE FUNCTIONS
   //do_quadcoinc(4, 8, t, ww_ns);
   do_quadcoinc(3, 19, t, ww_ns, runinfo);
   do_quadcoinc(4, 19, t, ww_ns, runinfo);
   do_quadcoinc(3, 20, t, ww_ns, runinfo);
   do_quadcoinc(4, 20, t, ww_ns, runinfo);
   //printf("did a coincidence analysis \n"); //testgg

   }; 

TAFlowEvent* Analyze(TARunInfo* runinfo, TMEvent* event, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DlTdcModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
      return flow;
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      if (!fFlags->fEnabled)
         return flow;

      DlTdcEventFlow* f = flow->Find<DlTdcEventFlow>();
      if (f) {
         AnalyzeTdcEvent(*f->fDlTdcEvent, runinfo);
      }
      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("MODULE - coinc, analyze special event \n");
      if (fTrace)
         printf("DlTdcModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
      printf("passed fTrace in AnalyzeSpecialEvent\n"); //testGG
   }
}; 

class DlTdc16CoincModuleFactory: public TAFactory
{
public:
   DlTdcFlags fFlags;

public:
   void Usage()
   {//ignore
      printf("DlTdc16CoincModuleFactory flags:\n");
      printf("--dltdc8 -- enable analysis of 4 paddle data\n");
      printf("--dltdc8-triggered -- analyze only events with hit in channel T (coincidence of A and B)\n");
      printf("--dltdc8-debug -- print detailed information\n");
      printf("--dltdc8-print -- print events\n");
      printf("--dltdc8-twc -- get twc W parameter");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("MODULE - coinc16, init \n");
      printf("DlTdc16CoincModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--coinc16") {
            fFlags.fEnabled = true;
         }
         if (args[i] == "--coinc16-triggered") {
            fFlags.fTriggered = true;
         }
         if (args[i] == "--coinc16-debug") {
            fFlags.fDebug = true;
         }
         if (args[i] == "--coinc16-print") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--coinc16-twc") {
            fFlags.fTWC = true;
         }
      }
   }

   void Finish()
   {
      printf("DlTdc16CoincModuleFactory::Finish!\n");
      // printf("GABBY GABBY GABBY GABBY GABBY GABBY GABBY GABBY GABBY GABBY \n");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("DlTdc16CoincModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //return new DlTdc8Module(runinfo, &fFlags);
      return new DlTdc16CoincModule(runinfo, &fFlags);
   }
};

static TARegister tar(new DlTdc16CoincModuleFactory);
