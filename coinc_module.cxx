// Gabby Gelinas, July 17, 2024
// Module for producing coincidence plots with 8 paddles (16 scintillators)

// ASK - I INCLUDED ALL OF THE INCLUDES FROM dltdc8_module.cxx BUT IDK IF THEY ARE ALL NEEDED/WHAT EACH IS FOR
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

// copied means copied from dltdc8

// copied
class DlTdcFlags
{
public:
   bool fEnabled = false;
   bool fTriggered = false;
   bool fDebug = false;
   bool fPrint = false;
   bool fTWC = false;
};


//class DlTdcMap8
class DlTdcCoincMap8
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

//bool DlTdcMap8::Init(int runno) //TEST - Change this and the class it's part of to be DlTdcCoincMap8? I thought since this is a mapping thing
// it may need to be consistent across files
bool DlTdcCoincMap8::Init(int runno)
{
   if (runno < 156) {
      printf("DlTdcMap8 old 1458 map for run %d!\n", runno);

      fMap.resize(16+1);
      
      fMap[1] = 16; // chan1
      fMap[2] = 17; // chan2
      fMap[3] = 22; // chan3
      fMap[4] = 23; // chan4
      fMap[5] = 26; // chan5
      fMap[6] = 27; // chan6
      fMap[7] = 30; // chan7
      fMap[8] = 31; // chan8
      
      fMap[9]  = 18;
      fMap[10] = 19;
      fMap[11] = 20;
      fMap[12] = 21;
      fMap[13] = 24;
      fMap[14] = 25;
      fMap[15] = 28;
      fMap[16] = 29;
      
      fPair1.resize(8+1);
      fPair2.resize(8+1);
      
      fPair1[1] =  1; fPair2[1] =  4; // chan14
      fPair1[2] =  2; fPair2[2] =  3; // chan23
      fPair1[3] =  5; fPair2[3] =  8; // chan58
      fPair1[4] =  6; fPair2[4] =  7; // chan67
      fPair1[5] =  9; fPair2[5] = 12;
      fPair1[6] = 10; fPair2[6] = 11;
      fPair1[7] = 13; fPair2[7] = 16;
      fPair1[8] = 14; fPair2[8] = 15;
   } else {
      printf("DlTdcMap8 for run %d!\n", runno);

      fMap.resize(16+1);
      
      fMap[1]  =  0; // chan1
      fMap[2]  =  1; // chan2
      fMap[3]  = 10; // chan3
      fMap[4]  = 11; // chan4
      fMap[5]  =  2 + 16; // chan5
      fMap[6]  =  3 + 16; // chan6
      fMap[7]  =  8 + 16; // chan7
      fMap[8]  =  9 + 16; // chan8
      
      fMap[9]  = 15;
      fMap[10] = 14;
      fMap[11] =  7;
      fMap[12] =  6;
      fMap[13] = 13 + 16;
      fMap[14] = 12 + 16;
      fMap[15] =  5 + 16;
      fMap[16] =  4 + 16;
      
      fPair1.resize(8+1);
      fPair2.resize(8+1);
      
      fPair1[1] =  1; fPair2[1] =  9; // chan14
      fPair1[2] =  2; fPair2[2] = 10; // chan23
      fPair1[3] =  3; fPair2[3] = 11; // chan58
      fPair1[4] =  4; fPair2[4] = 12; // chan67
      fPair1[5] =  5; fPair2[5] = 13;
      fPair1[6] =  6; fPair2[6] = 14;
      fPair1[7] =  7; fPair2[7] = 15;
      fPair1[8] =  8; fPair2[8] = 16;
   }

   return true;
};

//class DlTdc8Module: public TARunObject //change name
class DlTdc8CoincModule: public TARunObject
{ 
public:
   DlTdcFlags* fFlags = NULL;

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

   //TEST_Sept17
   // Make a vector of PAIRS (the lower numbered channel) where each pair is one coincidence you want to examine. 
   //We do this so you don't get an overwhelming number of coincidence plots
   std::vector<std::pair<int, int>> possibleCoinc = {std::make_pair(1, 5), std::make_pair(2, 6)};

   //std::vector<std::pair<int, int>> possibleCoinc;
   int tester = possibleCoinc[0].first;

   //Make the coincidence pairs, add as many as you need here
   //std::pair<int, int> p1, p2;
   //std::pair<int, int> p1 = std::make_pair(1, 5);
   //std::pair<int, int> p2 = std::make_pair(2, 6);

   //Add coincidence pairs to the vector, add as many as you need here
   //possibleCoinc.push_back(p1);
   //possibleCoinc.push_back(p2);


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

   char name_tofwidth_cha_papb[256];
   char name_tofwidth_cha_papb_cut[256];
   char name_tofwidth_cha_papb_cut_twc[256];

   char name_absTime_cha_papb[256];
   char name_absTime_cha_papb_twc[256];

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

   char title_tofwidth_cha_papb[256];
   char title_tofwidth_cha_papb_cut[256];
   char title_tofwidth_cha_papb_cut_twc[256];

   char title_absTimeWidth_cha_papb[256];
   char title_absTimeWidth_cha_papb_twc[256];

   double fPrevEventTimeSec = 0;

   //testing
   double timeDiff_pair;

   //DlTdcMap8 *fMap8 = NULL;
   DlTdcCoincMap8 *fMap8 = NULL;

   bool fTrace = false;

   Ncfm* fCfm = NULL;

   double width_cut = 10.0; // Width to cut below in ns

   double ww_twc = 9.0;

   // Store the channel reading at the ends of each scintillator (channel pair) together in a map, identified by the pair number
   std::map<int, std::pair<double, double>> scintPairs; // <pair number, (lower# channel time measurement, higher# channel time measurement)>

   
   // DlTdc8Module(TARunInfo* runinfo, DlTdcFlags* flags)
   DlTdc8CoincModule(TARunInfo* runinfo, DlTdcFlags* flags)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("DlTdc8Module::ctor!\n");

      fModuleName = "dltdc8_module";
      fFlags   = flags;

      fCfm = new Ncfm("dlcfmdb");

      //fMap8 = new DlTdcMap8();
      fMap8 = new DlTdcCoincMap8();
   }

   //~DlTdc8Module()
   ~DlTdc8CoincModule()
   {
      if (fTrace)
         printf("DlTdc8Module::dtor!\n");

      if (fMap8) {
         delete fMap8;
         fMap8 = NULL;
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

      bool conf_ok = fMap8->Init(runinfo->fRunNo);
      if (!conf_ok) {
         printf("Cannot load TDC map for run %d\n", runinfo->fRunNo);
         exit(123);
      }

      runinfo->fRoot->fOutputFile->cd(); // select correct ROOT directory
      TDirectory* dir = gDirectory->mkdir("coinc");
      dir->cd(); // select correct ROOT directory

      // Initialize histograms
      int i;
      int j;

      // Initialize histograms for plots with times along scintillators, widths, and single channels
      // This will make lots of enteries in these maps that won't be used like pairs with the same number (1,1) or ones like (1,9) when we will only use (9,1)
      // Filtering so we only fill the enteries of interest happens later.

      //TESTING_Sept17
      // testing a way to only make the histograms we care about for this set up of coincidences
      for (int coincNum=0; possibleCoinc.size(); coincNum++){
         // Make the plot names and titles for each plot for the coincidence pairing of interest

         // Make a vector of each entry in this coincidence pairing so I can put it in a for loop 
         int coincPair1=possibleCoinc[coincNum].first; // Pair as in "pair" meaning scintillator. So what scintillator
         int coincPair2=possibleCoinc[coincNum].second;

         std::vector<int> thisCoincScint={coincPair1, coincPair2}; //Tracks which scintillators are in this coincidence

         // Doesn't specify a single scintillator or channel of interest so this can happen outside of the second for loop
         sprintf(name_tof_papb, "coinc_tof_cpair_%02d_%02d", coincPair1,coincPair2);
         sprintf(name_tof_papb_cut, "coinc_tof_cpair_%02d_%02d_cut", coincPair1,coincPair2);
         sprintf(name_tof_papb_cut_twc, "coinc_tof_cpair_%02d_%02d_cut_twc", coincPair1,coincPair2);

         sprintf(title_tof_papb, "TOF (ns) coinc with pairs %02d %02d", coincPair1,coincPair2);
         sprintf(title_tof_papb_cut, "TOF (ns) coinc with pairs %02d %02d cut", coincPair1,coincPair2);
         sprintf(title_tof_papb_cut_twc, "TOF (ns) coinc with pairs %02d %02d cut twc", coincPair1,coincPair2);


         for (int i=0; thisCoincScint.size(); i++){
            int scintOfInterest = thisCoincScint[i];
            if (i==0){
               int otherScint = thisCoincScint[1];
            }
            else {
               int otherScint = thisCoincScint[0];
            }

            // Make plot names
            sprintf(name_timeDiff_pa_papb, "coinc_tpair%d_ns_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_timeDiff_pa_papb_cut, "coinc_tpair%d_ns_cpair_%02d_%02d_cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_timeDiff_pa_papb_cut_twc, "coinc_tpair%d_ns_cpair_%02d_%02d_cut_twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_tof_papb, "coinc_tof_cpair_%02d_%02d", coincPair1,coincPair2);
            sprintf(name_tof_papb_cut, "coinc_tof_cpair_%02d_%02d_cut", coincPair1,coincPair2);
            sprintf(name_tof_papb_cut_twc, "coinc_tof_cpair_%02d_%02d_cut_twc", coincPair1,coincPair2);

            sprintf(name_widths_pa_papb, "coinc_width%d_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_widths_pa_papb, "coinc_width%d_cpair_%02d_%02d", scintOfInterest+8, coincPair1,coincPair2); //new

            sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", scintOfInterest+8, coincPair1,coincPair2); //new
            sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", scintOfInterest+8, coincPair1,coincPair2); //new
            sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", scintOfInterest+8, coincPair1,coincPair2);//new

            sprintf(name_absTime_cha_papb, "coinc_absTime_ch%d_cpair_%02d_%02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(name_absTime_cha_papb_twc, "coinc_absTime_ch%d_cpair_%02d_%02d_twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(name_absTime_cha_papb, "coinc_absTime_ch%d_cpair_%02d_%02d", scintOfInterest+8, coincPair1,coincPair2); //new
            sprintf(name_absTime_cha_papb_twc, "coinc_absTime_ch%d_cpair_%02d_%02d_twc", scintOfInterest+8, coincPair1,coincPair2); //new
            // All plus 8s need to come after the histogram is made

            sprintf(name_absTimeWidth_cha_papb, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d", scintOfInterest, scintOfInterest, coincPair1,coincPair2);
            sprintf(name_absTimeWidth_cha_papb_twc, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d_twc", scintOfInterest, scintOfInterest, coincPair1,coincPair2);

            sprintf(name_absTimeWidth_cha_papb, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d", scintOfInterest+8, scintOfInterest+8, coincPair1,coincPair2); //new
            sprintf(name_absTimeWidth_cha_papb_twc, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d_twc", scintOfInterest+8, i+8, coincPair1,coincPair2); //new


            // Make plot titles
            sprintf(title_timeDiff_pa_papb, "Time difference along pair %d (ns) with coinc between pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_timeDiff_pa_papb_cut, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_timeDiff_pa_papb_cut_twc, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(title_tof_papb, "TOF with coincidence between pairs %02d %02d", coincPair1,coincPair2);
            sprintf(title_tof_papb_cut, "TOF with coincidence between pairs %02d %02d cut", coincPair1,coincPair2);
            sprintf(title_tof_papb_cut_twc, "TOF with coincidence between pairs %02d %02d cut twc", coincPair1,coincPair2);

            sprintf(title_widths_pa_papb, "Width on channel %d, coinc between pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_widths_pa_papb, "Width on channel %d, coinc between pairs %02d %02d", scintOfInterest+8, coincPair1,coincPair2); //new

            sprintf(title_tofwidth_cha_papb, "TOF vs Width on channel%d with coinc between pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut_twc, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(title_tofwidth_cha_papb, "TOF vs Width on channel%d with coinc between pairs %02d %02d", scintOfInterest+8, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut", scintOfInterest+8, coincPair1,coincPair2);
            sprintf(title_tofwidth_cha_papb_cut_twc, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut twc", scintOfInterest+8, coincPair1,coincPair2);

            sprintf(title_absTime_chapb, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d", scintOfInterest, otherScint, otherScint+8,coincPair1,coincPair2);
            sprintf(title_absTime_chapb_cut_twc, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d cut twc", scintOfInterest, otherScint, otherScint+8,coincPair1,coincPair2);

            sprintf(title_absTime_chapb, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d", coincPair1+8, otherScint, otherScint+8,coincPair1,coincPair2);
            sprintf(title_absTime_chapb_cut_twc, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d cut twc", coincPair1+8, otherScint, otherScint+8,coincPair1,coincPair2);

            sprintf(title_absTimeWidth_cha_papb, "TOF vs width%d coinc btwn pairs %02d %02d", scintOfInterest, coincPair1,coincPair2);
            sprintf(title_absTimeWidth_cha_papb_twc, "TOF vs width%d coinc btwn pairs %02d %02d twc", scintOfInterest, coincPair1,coincPair2);

            sprintf(title_absTimeWidth_cha_papb, "TOF vs width%d coinc btwn pairs %02d %02d", scintOfInterest+8, coincPair1,coincPair2);
            sprintf(title_absTimeWidth_cha_papb_twc, "TOF vs width%d coinc btwn pairs %02d %02d twc", scintOfInterest+8, coincPair1,coincPair2);


            // Initialize histograms





            // Make plot names for plots that focus on the higher number channel on a scintillator



            // Make plot titles for plots that focus on the higher number channel on a scintillator


            // Initalize histograms for plots that focus on the higher number channel on a scintillator

         }
      }



      #if 0
      // Known way that works but produces SO many unnecessary plots
      for (i=1; i<17; i++){
          for(j=1; j<17; j++){
            // tesing
            //printf("Making names and titles for histograms \n");

            // Create the names of each plot
            if ((i < 9) && (j<9)) { // Only make file names once for the plots that are pair specific and not channel specific
               // This section worked when it was all out of an if statement, it just produced so many unneeded plots
               sprintf(name_timeDiff_pa_papb, "coinc_tpair%d_ns_cpair_%02d_%02d", i, i,j);
               sprintf(name_timeDiff_pa_papb_cut, "coinc_tpair%d_ns_cpair_%02d_%02d_cut", i, i,j);
               sprintf(name_timeDiff_pa_papb_cut_twc, "coinc_tpair%d_ns_cpair_%02d_%02d_cut_twc", i, i,j);

               sprintf(name_tof_papb, "coinc_tof_cpair_%02d_%02d", i, j);
               sprintf(name_tof_papb_cut, "coinc_tof_cpair_%02d_%02d_cut", i, j);
               sprintf(name_tof_papb_cut_twc, "coinc_tof_cpair_%02d_%02d_cut_twc", i, j);

               sprintf(name_timeDiff_pa_pb_papb, "coinc_tpair%d_tpair%d_cpair_%02d_%02d", i, j, i, j);
               sprintf(name_timeDiff_pa_pb_papb_cut, "coinc_tpair%d_tpair%d_cpair_%02d_%02d_cut", i, j, i, j);
               sprintf(name_timeDiff_pa_pb_papb_cut_twc, "coinc_tpair%d_tpair%d_cpair_%02d_%02d_cut_twc", i, j, i, j);

               sprintf(name_widths_pa_papb, "coinc_width%d_cpair_%02d_%02d", i, i, j);
               sprintf(name_widths_pa_papb, "coinc_width%d_cpair_%02d_%02d", i+8, i, j); //new

               sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", i, i, j);
               sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", i, i, j);
               sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", i, i, j);

               sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", i+8, i, j); //new
               sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", i+8, i, j); //new
               sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", i+8, i, j);//new

               sprintf(name_absTimeWidth_cha_papb, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d", i, i, i, j);
               sprintf(name_absTimeWidth_cha_papb_twc, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d_twc", i, i, i, j);

               sprintf(name_absTimeWidth_cha_papb, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d", i+8, i+8, i, j); //new
               sprintf(name_absTimeWidth_cha_papb_twc, "coinc_absTime_ch%d_width%d_cpair_%02d_%02d_twc", i+8, i+8, i, j); //new
            }

            if (j<9){
               // Need this one to go for every i but only the pair number (1-8) j since we should always be along a scintillator for j
               sprintf(name_absTime_chapb, "coinc_t%d_%d_%d_cpair_%02d_%02d", i, j, j+8,i,j);
               sprintf(name_absTime_chapb_cut_twc, "coinc_t%d_%d_%d_cpair_%02d_%02d_cut_twc", i, j, j+8,i,j);
            }

            
            // Create the titles for each plot
              sprintf(title_timeDiff_pa_papb, "Time difference along pair %d (ns) with coinc between pairs %02d %02d", i, i,j);
              sprintf(title_timeDiff_pa_papb_cut, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut", i, i,j);
              sprintf(title_timeDiff_pa_papb_cut_twc, "Time difference along pair %d (ns) with coinc between pairs %02d %02d cut twc", i, i,j);

              sprintf(title_tof_papb, "TOF (ns) coinc with pairs %02d %02d", i, j);
              sprintf(title_tof_papb_cut, "TOF (ns) coinc with pairs %02d %02d cut", i, j);
              sprintf(title_tof_papb_cut_twc, "TOF (ns) coinc with pairs %02d %02d cut twc", i, j);

              sprintf(title_absTime_chapb, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d", i, j,i,j);
              sprintf(title_absTime_chapb_cut_twc, "t%d - 0.5(t%d+t%d) (ns) coinc btwn pairs %02d %02d cut twc", i, j, j,i,j);

              sprintf(title_timeDiff_pa_pb_papb, "Time differnce along pair%d and pair%d coinc (ns)", i, j);
              sprintf(title_timeDiff_pa_pb_papb_cut, "Time differnce along pair%d and pair%d coinc (ns) cut", i, j);
              sprintf(title_timeDiff_pa_pb_papb_cut_twc, "Time differnce along pair%d and pair%d coinc (ns) cut twc", i, j);

              sprintf(title_tofwidth_cha_papb, "TOF vs Width on channel%d with coinc between pairs %02d %02d", i, i, j);
              sprintf(title_tofwidth_cha_papb_cut, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut", i, i, j);
              sprintf(title_tofwidth_cha_papb_cut_twc, "TOF vs Width on channel%d with coinc between pairs %02d %02d cut twc", i, i, j);

              sprintf(title_absTimeWidth_cha_papb, "TOF vs width%d coinc btwn pairs %02d %02d", i, i, j);
              sprintf(title_absTimeWidth_cha_papb_twc, "TOF vs width%d coinc btwn pairs %02d %02d twc", i, i, j);

            // tesing
            //printf("Done making names and titles for histograms \n");

            // Initialize a pair variable for storing the information of the scintillators we are examining for accessing our map element
            std::pair<int,int> scintPair (i,j);

              // Initalize 1D histograms
              timeDiff_pa_papb[scintPair] = new TH1D(name_timeDiff_pa_papb, title_timeDiff_pa_papb, 200, -10, 10);
              timeDiff_pa_papb_cut[scintPair] = new TH1D(name_timeDiff_pa_papb_cut, title_timeDiff_pa_papb_cut, 200, -10, 10);
              timeDiff_pa_papb_cut_twc[scintPair] = new TH1D(name_timeDiff_pa_papb_cut_twc, title_timeDiff_pa_papb_cut_twc, 200, -10, 10);

              tof_papb[scintPair] = new TH1D(name_tof_papb, title_tof_papb, 200, -10, 10);
              tof_papb_cut[scintPair] = new TH1D(name_tof_papb_cut, title_tof_papb_cut, 200, -10, 10);
              tof_papb_cut_twc[scintPair] = new TH1D(name_tof_papb_cut_twc, title_tof_papb_cut_twc, 200, -10, 10);

              absTime_chapb[scintPair] = new TH1D(name_absTime_chapb, title_absTime_chapb, 200, -10, 10);
              absTime_chapb_cut_twc[scintPair] = new TH1D(name_absTime_chapb_cut_twc, title_absTime_chapb_cut_twc, 200, -10, 10);


              // Initalize 2D histograms
              timeDiff_pa_pb_papb[scintPair] = new TH2D(name_timeDiff_pa_pb_papb, title_timeDiff_pa_pb_papb, 200, -10, 10, 200, -10, 10);
              timeDiff_pa_pb_papb_cut[scintPair] = new TH2D(name_timeDiff_pa_pb_papb_cut, title_timeDiff_pa_pb_papb_cut, 200, -10, 10, 200, -10, 10);
              timeDiff_pa_pb_papb_cut_twc[scintPair] = new TH2D(name_timeDiff_pa_pb_papb_cut_twc, title_timeDiff_pa_pb_papb_cut_twc, 200, -10, 10, 200, -10, 10);

              widths_pa_papb[scintPair] = new TH2D(name_widths_pa_papb, title_widths_pa_papb, 200, -10, 10, 200, -10, 10);

              tofwidth_cha_papb[scintPair] = new TH2D(name_tofwidth_cha_papb, title_tofwidth_cha_papb, 200, -10, 10, 200, -10, 10);
              tofwidth_cha_papb_cut[scintPair] = new TH2D(name_tofwidth_cha_papb_cut, title_tofwidth_cha_papb_cut, 200, -10, 10, 200, -10, 10);
              tofwidth_cha_papb_cut_twc[scintPair] = new TH2D(name_tofwidth_cha_papb_cut_twc, title_tofwidth_cha_papb_cut_twc, 200, -10, 10, 200, -10, 10);

              absTimeWidth_cha_papb[scintPair] = new TH2D(name_absTimeWidth_cha_papb, title_absTimeWidth_cha_papb, 200, -10, 10, 200, -10, 10);
              absTimeWidth_cha_papb_twc[scintPair] = new TH2D(name_absTimeWidth_cha_papb_twc, title_absTimeWidth_cha_papb_twc, 200, -10, 10, 200, -10, 10);
          }
      }
      #endif

      // Initialize histograms for plots involving single channels for the higher numbered channels in the pair
      for (i=9; i<17; i++){
          for (j=1; j<9; j++){
              sprintf(name_absTime_chapb, "coinc_t%d_%d_%d_cpair_%02d_%02d", i, j, j+8,i,j);
              sprintf(name_absTime_chapb_cut_twc, "coinc_t%d_%d_%d_cpair_%02d_%02d_cut_twc", i, j, j+8,i,j);

              sprintf(name_tofwidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", i, i, j);
              sprintf(name_tofwidth_cha_papb_cut, "coinc_tofwidth%d_cpair_%02d_%02d_cut", i, i, j);
              sprintf(name_tofwidth_cha_papb_cut_twc, "coinc_tofwidth%d_cpair_%02d_%02d_cut_twc", i, i, j);

              sprintf(name_absTimeWidth_cha_papb, "coinc_tofwidth%d_cpair_%02d_%02d", i, i, j);
              sprintf(name_absTimeWidth_cha_papb_twc, "coinc_tofwidth%d_cpair_%02d_%02d_twc", i, i, j);


          }
      }

      // Store the channel reading at the ends of each scintillator (channel pair) together in a map, identified by the pair number
      //std::map<int, std::pair<double, double>> scintPairs; // <pair number, (lower# channel time measurement, higher# channel time measurement)>

      // TO DO: ASK how the fMap works with calling to channels
      for (int i=1; i<9; i++){
         scintPairs.insert({i, std::make_pair(fMap8->fMap[i],fMap8->fMap[i+8])});
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
void do_quadcoinc(int scint_a, int scint_b, const DlTdcEvent& t, std::vector<double> ww_ns)
   {
   //printf("MODULE - coinc, doing quad coinc GABBY\n");  
    /* Fill histograms for the appropriate quad coincidences
    scint_a - one scintillator pair number (1-8) you want to examine, int
    scint_b - the other scintillator pair number (1-8) you want to examine, int

    NOTE: a and b are different from A and B used to denote the two data cables
    */

    // Time difference along each scintillator (channel a+8) - (channel a)
    //printf("scint_a %d and scint_b %d \n", scint_a, scint_b);
    //int scint_a8 = scint_a+8;
    //printf("scint_a+8=%d \n", scint_a8);
    //printf("test 1, %f \n", t.GetCh(scintPairs[scint_a].second).fLe);
    double t_scinta_ns = subtract_ns(t.GetCh(scintPairs[scint_a].second).fLe, t.GetCh(scintPairs[scint_a].first).fLe);
    double t_scintb_ns = subtract_ns(t.GetCh(scintPairs[scint_b].second).fLe, t.GetCh(scintPairs[scint_b].first).fLe);


    if (t_scinta_ns > -9999 && t_scintb_ns > -9999){
         // testing: Makes it into this statement
 
         // Make a vector with int entries of scint_a, scint_b, and those two +8. We will use this list to fill the above calculation vectors.
         std::vector<int> channelNumbers;
         channelNumbers.push_back(scint_a);
         channelNumbers.push_back(scint_a+8);
         channelNumbers.push_back(scint_b);
         channelNumbers.push_back(scint_b+8);
 
         // testing
         //printf("\n scint_a=%d \n", scint_a);
         //printf("Channel numbers 0, %d,num1=%d, num2=%d, number 3, %d", channelNumbers[0], channelNumbers[1],channelNumbers[2],channelNumbers[3]);

         for (int i=0; i<=3; ++i){
            for (int j=0; j<=3; ++j){
               // testing: Making it into this loop correctly
               int chan_i = channelNumbers[i];
               //int chan_i8 = channelNumbers[i] + 8;
               int chan_j = channelNumbers[j];
               int chan_j8 = channelNumbers[j] + 8;

               // Pair of the channels we are examining for this coincidence
               std::pair<int,int> chanPair_ij (channelNumbers[i],channelNumbers[j]);
               std::pair<int,int> chanPair_ji (channelNumbers[j],channelNumbers[i]);

               //if (channelNumbers[i] < 9){
               //   std::pair<int,int> chanPair_i8j (channelNumbers[i]+8,channelNumbers[j]);
               //}
               //else{
               //   std::pair<int,int> chanPair_i8j (channelNumbers[i],channelNumbers[j]);
               //}

               std::pair<int,int> chanPair_i8j (channelNumbers[i]+8,channelNumbers[j]);
               std::pair<int,int> chanPair_j8i (channelNumbers[j]+8,channelNumbers[i]);
               //std::cout<< "chanPair_ij:" << chanPair_ij.first << " " << chanPair_ij.second<<std::endl;
               //std::cout<< "chanPair_ji:" <<chanPair_ji.first  << " " << chanPair_ji.second<<std::endl;
               //std::cout<< "chanPair_i8j:" <<chanPair_i8j.first  << " " << chanPair_i8j.second<<std::endl;


                // Breaking it up into if statements here so we don't calculate the subtraction of a channel with itself or do the same calculation
                // twice but reversed (i.e. t9-t1 and t1-t9). Set so index i tracks the higher numbered channel
                // chan_i < 17 keeps us within the number of scintillator channels we actually have. No fake channels
                if ((chan_i==chan_j8) && (chan_i < 17)) {

                    //Calculations for time differences along a single scintillator
                    // testing: makes it here
                    //printf("Along a scintillator %d \n", chan_i);
                    //double timeDiff_pair = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].second).fLe, t.GetCh(scintPairs[channelNumbers[i]].first).fLe);
                    //subtract_alongPair.push_back({(channelNumbers[i], channelNumbers[j]), timeDiff_pair});
                    timeDiff_pair = subtract_ns(t.GetCh(scintPairs[channelNumbers[i]].second).fLe, t.GetCh(scintPairs[channelNumbers[i]].first).fLe);
                    
                    
                    double timeDiff_pair_twc = timeDiff_pair + (ww_twc / sqrt(ww_ns[channelNumbers[j]]) - ww_twc / sqrt(ww_ns[channelNumbers[i]]));
                    //subtract_alongPair_twc.push_back({(channelNumbers[i], channelNumbers[j]), timeDiff_pair_twc});
                     //printf("time diff calculations done \n");

                     //for (const auto & [key, value] : timeDiff_pa_papb){
                     //      std::cout<<key.first << ' ' <<key.second<<std::endl;
                     //      std::cout << value << std::endl;
                     //}

                    // Fill histograms
                    if (ww_ns[channelNumbers[i]] > width_cut && ww_ns[channelNumbers[i]+8] > width_cut && ww_ns[channelNumbers[j]] > width_cut && ww_ns[channelNumbers[j]+8] > width_cut) {
                        timeDiff_pa_papb_cut[chanPair_ij]->Fill(timeDiff_pair);
                        timeDiff_pa_papb_cut_twc[chanPair_ij]->Fill(timeDiff_pair_twc);
                        //printf("Filled histograms \n");
                    }

                    else { // No width cut
                        //tesing
                        //printf("Can enter no width cut statement \n");
                        //printf("time diff pair %f", timeDiff_pair);
                        //std::cout<<"Second part of chanPair "<< chanPair_ji.second<<std::endl;
                        //std::cout<< "First part of chanPair "<<chanPair_ji.first<<std::endl;
                        //std::cout<<timeDiff_pa_papb[chanPair_ji]<<std::endl;
                        timeDiff_pa_papb[chanPair_ji]->Fill(timeDiff_pair);
                        //printf("Filed no width cut histogram \n");
                    }
                }

                if ((i != j) && (chan_i < 9) && (chan_j < 9) && (chan_i > chan_j)){
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

                    double tof_a0a1b0b1 = -(timeDiff_a0_b0 + timeDiff_a1_b1);
                    double tof_a0a1b0b1_twc = tof_a0a1b0b1 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[i]]) + ww_twc / sqrt(ww_ns[channelNumbers[i]+8]) - ww_twc / sqrt(ww_ns[channelNumbers[j]]) - ww_twc / sqrt(ww_ns[channelNumbers[j]+8])); 
                    //printf("Calculates second block \n");

                    double absTime_a0 = 0.5*(timeDiff_a0_b0 + timeDiff_a0_b1);
                    double absTime_a1 = 0.5*(timeDiff_a1_b0 + timeDiff_a1_b1);
                    double absTime_b0 = -0.5*(timeDiff_a0_b0 + timeDiff_a1_b0);
                    double absTime_b1 = -0.5*(timeDiff_a0_b1 + timeDiff_a1_b1);
                    //printf("Calculates third block \n");

                    double absTime_a0_twc = absTime_a0 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[j]]) + ww_twc / sqrt(ww_ns[channelNumbers[j]+8]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[i]]));
                    double absTime_a1_twc = absTime_a1 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[j]]) + ww_twc / sqrt(ww_ns[channelNumbers[j]+8]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[i]+8]));
                    double absTime_b0_twc = absTime_b0 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[i]]) + ww_twc / sqrt(ww_ns[channelNumbers[i]+8]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[j]]));
                    double absTime_b1_twc = absTime_b1 + 0.5*(ww_twc / sqrt(ww_ns[channelNumbers[i]]) + ww_twc / sqrt(ww_ns[channelNumbers[i]+8]) - 2*ww_twc / sqrt(ww_ns[channelNumbers[j]+8]));
                    //printf("Calculates fourth block \n");


                    // Fill histograms with a width cut
                    if (ww_ns[channelNumbers[i]] > width_cut && ww_ns[channelNumbers[i]+8] > width_cut && ww_ns[channelNumbers[j]] > width_cut && ww_ns[channelNumbers[j]+8] > width_cut) {
                        //printf("Enters different scintillators width cut \n");
                        tof_papb_cut[chanPair_ij]->Fill(tof_a0a1b0b1);
                        tof_papb_cut_twc[chanPair_ij]->Fill(tof_a0a1b0b1_twc);
                        //printf("First block width calculations done \n");

                        absTime_chapb_cut_twc[chanPair_ij]->Fill(absTime_a0_twc);
                        absTime_chapb_cut_twc[chanPair_i8j]->Fill(absTime_a1_twc);
                        absTime_chapb_cut_twc[chanPair_ji]->Fill(absTime_b0_twc);
                        absTime_chapb_cut_twc[chanPair_j8i]->Fill(absTime_b1_twc);
                        //printf("Second block width calculations done \n");

                        timeDiff_pa_pb_papb_cut[chanPair_ij]->Fill(t_scinta_ns,t_scintb_ns);
                        //printf("Third block width calculations done \n");

                        tofwidth_cha_papb_cut[chanPair_ij]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]]);
                        tofwidth_cha_papb_cut[chanPair_i8j]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]+8]);
                        tofwidth_cha_papb_cut[chanPair_ji]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]]);
                        tofwidth_cha_papb_cut[chanPair_j8i]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]+8]);
                        //printf("Fourth block width calculations done \n");

                    }

                    else { // No width cut
                        //printf("Enters different scintillators NO width cut \n");
                        tof_papb[chanPair_ij]->Fill(tof_a0a1b0b1);
                        printf("First block no width cut calculations done \n");
                        std::cout << "chanPair_ij " << chanPair_ij.first << " " << chanPair_ij.second << std::endl;
                        std::cout << "absTime_a0 " << absTime_a0 << std::endl;
                        //std::cout << "absTime_chapb " << absTime_chapb << std::endl;

                        absTime_chapb[chanPair_ij]->Fill(absTime_a0);
                        //printf("Second block ij done \n");
                        //std::cout << "chanPair_i8j " << chanPair_i8j.first << " " << chanPair_i8j.second << std::endl;
                        absTime_chapb[chanPair_i8j]->Fill(absTime_a1);
                        //printf("Second block i8j done \n");
                        absTime_chapb[chanPair_ji]->Fill(absTime_b0);
                        //printf("Second block ji done \n");
                        absTime_chapb[chanPair_j8i]->Fill(absTime_b1);
                        //printf("Second block no width cut calculations done \n");

                        timeDiff_pa_pb_papb[chanPair_ij]->Fill(t_scinta_ns,t_scintb_ns);
                        //printf("Third block no width cut calculations done \n");

                        widths_pa_papb[chanPair_ij]->Fill(ww_ns[channelNumbers[i]],ww_ns[channelNumbers[j]]);
                        //printf("Fourth block no width cut calculations done \n");

                        tofwidth_cha_papb[chanPair_ij]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]]);
                        tofwidth_cha_papb[chanPair_i8j]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[i]+8]);
                        tofwidth_cha_papb[chanPair_ji]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]]);
                        tofwidth_cha_papb[chanPair_j8i]->Fill(tof_a0a1b0b1, ww_ns[channelNumbers[j]+8]);
                        //printf("Fifth block no width cut calculations done \n");
                    }
                    
                    // delete
                    /*

                    // 2D time difference
                    std::vector<(int,int),TH2D*> timeDiff_pa_pb_papb_cut_twc (24);

                    // TOF vs width
                    std::vector<(int,int),TH2D*> tofwidth_cha_papb_cut_twc (24);

                    // Absolute time vs width
                    std::vector<(int,int),TH2D*> absTimeWidth_cha_papb (24);
                    std::vector<(int,int),TH2D*> absTimeWidth_cha_papb_twc (24);
                    */

                }
            }
        }
    }
}


// Call function for each possible quad coincidence in this function
void AnalyzeTdcEvent(const DlTdcEvent& t)
   {
      if (fFlags->fDebug) {
         std::string s = "";
         for (size_t i=1; i<fMap8->fMap.size(); i++) {
            int tdc_ch = fMap8->fMap[i];
            if (t.HaveCh(tdc_ch))
               s += "H";
            else
               s += ".";
         }

         printf("EVENT %s, ABT %d%d%d\n", s.c_str(), t.HaveCh(fMap8->fChanA), t.HaveCh(fMap8->fChanB), t.HaveCh(fMap8->fChanT));
         printf("fMap", fMap8->fMap[1]);
      }

      ///////// check for triggered event ///////////
         
      if (fFlags->fTriggered && !t.HaveCh(fMap8->fChanT)) {
         return;
      }

///////// SET WIDTH CUTOFFS (NS)  /////////
      
      //double ww_cut_ns = 10;

      ///////// COMPUTE WIDTH AND PULSE HEIGHT ///////////

      std::vector<double> ww_ns(fMap8->fMap.size() + 1);

      for (size_t i=1; i<fMap8->fMap.size(); i++) {
         int tdc_ch = fMap8->fMap[i];

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

      // SET UP TIME WALK CORRECT VALUES //
      //double ww_twc = 9.0;

   // ALL MY COINCIDENCE FUNCTIONS
   do_quadcoinc(1, 5, t, ww_ns);

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
         AnalyzeTdcEvent(*f->fDlTdcEvent);
      }
      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      printf("MODULE - coinc, analyze special event \n");
      if (fTrace)
         printf("DlTdcModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
}; 

class DlTdc8CoincModuleFactory: public TAFactory
{
public:
   DlTdcFlags fFlags;

public:
   void Usage()
   {//ignore
      printf("DlTdc8CoincModuleFactory flags:\n");
      printf("--dltdc8 -- enable analysis of 4 paddle data\n");
      printf("--dltdc8-triggered -- analyze only events with hit in channel T (coincidence of A and B)\n");
      printf("--dltdc8-debug -- print detailed information\n");
      printf("--dltdc8-print -- print events\n");
      printf("--dltdc8-twc -- get twc W parameter");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("MODULE - coinc, init \n");
      printf("DlTdc8CoincModuleFactory::Init!\n");

      for (unsigned i=0; i<args.size(); i++) {
         if (args[i] == "--coinc") {
            fFlags.fEnabled = true;
         }
         if (args[i] == "--coinc-triggered") {
            fFlags.fTriggered = true;
         }
         if (args[i] == "--coinc-debug") {
            fFlags.fDebug = true;
         }
         if (args[i] == "--coinc-print") {
            fFlags.fPrint = true;
         }
         if (args[i] == "--coinc-twc") {
            fFlags.fTWC = true;
         }
      }
   }

   void Finish()
   {
      printf("DlTdc8CoincModuleFactory::Finish!\n");
      printf("GABBY GABBY GABBY GABBY GABBY GABBY GABBY GABBY GABBY GABBY");
   }
   
   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("DlTdc8CoincModuleFactory::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //return new DlTdc8Module(runinfo, &fFlags);
      return new DlTdc8CoincModule(runinfo, &fFlags);
   }
};

static TARegister tar(new DlTdc8CoincModuleFactory);
