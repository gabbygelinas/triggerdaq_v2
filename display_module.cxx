//
// display_module.cxx
//
// display analysis of TPC data
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
#include "TMath.h"
#include "TGraphPolar.h"

#include "AgFlow.h"
#include "ko_limits.h"

// histogram limit for number of hits in aw and pads
#define MAX_HITS 250

#define PLOT_MIN_TIME (900.0)
#define PLOT_MAX_TIME (5100.0)

#define NUM_PREAMP 32

#define DELETE(x) if (x) { delete (x); (x) = NULL; }

#define MEMZERO(p) memset((p), 0, sizeof(p))

class DisplayModule: public TARunObject
{
public:
   bool fPrint = false;

public:
   TCanvas* fC = NULL;

   bool fTrace = false;
   bool fDoEventDisplay = false;

   DisplayModule(TARunInfo* runinfo, bool do_print, bool do_event_display)
      : TARunObject(runinfo)
   {
      if (fTrace)
         printf("DisplayModule::ctor!\n");

      fPrint = do_print;
      fDoEventDisplay = do_event_display;

      if (fDoEventDisplay) {
         fC = new TCanvas("fC", "event display", 800, 800);
      }
   }

   ~DisplayModule()
   {
      if (fTrace)
         printf("DisplayModule::dtor!\n");
      DELETE(fC);
   }

   void BeginRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DisplayModule::BeginRun, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      //time_t run_start_time = runinfo->fOdb->odbReadUint32("/Runinfo/Start time binary", 0, 0);
      //printf("ODB Run start time: %d: %s", (int)run_start_time, ctime(&run_start_time));
   }

   void EndRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DisplayModule::EndRun, run %d\n", runinfo->fRunNo);
      //time_t run_stop_time = runinfo->fOdb->odbReadUint32("/Runinfo/Stop time binary", 0, 0);
      //printf("ODB Run stop time: %d: %s", (int)run_stop_time, ctime(&run_stop_time));
   }

   void PauseRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DisplayModule::PauseRun, run %d\n", runinfo->fRunNo);
   }

   void ResumeRun(TARunInfo* runinfo)
   {
      if (fTrace)
         printf("DisplayModule::ResumeRun, run %d\n", runinfo->fRunNo);
   }

   TAFlowEvent* AnalyzeFlowEvent(TARunInfo* runinfo, TAFlags* flags, TAFlowEvent* flow)
   {
      //printf("DisplayModule::Analyze, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);

      AgEventFlow *ef = flow->Find<AgEventFlow>();

      if (!ef || !ef->fEvent)
         return flow;

      AgAwHitsFlow* eawh = flow->Find<AgAwHitsFlow>();
      AgPadHitsFlow* eph = flow->Find<AgPadHitsFlow>();
      AgBscAdcHitsFlow* eba = flow->Find<AgBscAdcHitsFlow>();

      //int force_plot = false;

      AgEvent* age = ef->fEvent;

      //uint32_t adc16_coinc_dff = 0;
      uint32_t aw16_prompt = 0;
      uint32_t trig_bitmap  = 0;

      //if (age->trig && age->trig->udpData.size() > 7) {
      //   adc16_coinc_dff = (age->trig->udpData[6]>>8)&0xFFFF;
      //}

      if (age->trig && age->trig->udpData.size() > 9) {
         aw16_prompt = (age->trig->udpData[9])&0xFFFF;
      }

      if (age->trig && age->trig->udpData.size() > 7) {
         trig_bitmap = age->trig->udpData[6];
      }

      int trig_counter = -1;
      int adc_counter = -1;
      int pwb_counter = -1;
      int tdc_counter = -1;

      if (age->trig) {
         trig_counter = age->trig->counter;
      }

      if (age->a16) {
         adc_counter = age->a16->counter;
      }

      if (age->feam) {
         pwb_counter = age->feam->counter;
      }

      if (age->tdc) {
         tdc_counter = age->tdc->counter;
      }

      if (fPrint) {
         printf("Have AgEvent: %d %d %d %d\n", trig_counter, adc_counter, pwb_counter, tdc_counter);
      }

      bool do_plot = (runinfo->fRoot->fgApp != NULL);

      if (fDoEventDisplay && do_plot) {
         if (fC) {
            TVirtualPad* save_gpad = gPad;

            
            fC->Clear();
            fC->Divide(2,2);

            TVirtualPad *p = fC->cd(1);

            p->Divide(1, 2);

            TDirectory* dir = runinfo->fRoot->fgDir;
            dir->cd();

            p->cd(1);
            TH1D* hh = new TH1D("hh_aw_wire_time", "AW wire hit drift time; TPC wire number; drift time, ns", NUM_AW, -0.5, NUM_AW-0.5);
            hh->SetMinimum(0);
            hh->SetMaximum(MAX_TIME);
            hh->Draw();

            p->cd(2);
            TH1D* ha = new TH1D("hh_aw_wire_amp", "AW wire hit amplitude; TPC wire number; AW ADC counts", NUM_AW, -0.5, NUM_AW-0.5);
            ha->SetMinimum(0);
            ha->SetMaximum(66000);
            ha->Draw();

            TVirtualPad *p_pad = fC->cd(2);

            p_pad->Divide(1, 2);

            std::vector<Double_t> zpad_col[NUM_PC];
            std::vector<Double_t> zpad_row[NUM_PC];
            std::vector<Double_t> zpad_time[NUM_PC];
            std::vector<Double_t> zpad_amp[NUM_PC];
               
            std::vector<Double_t> pad_col;
            std::vector<Double_t> pad_row;
            std::vector<Double_t> pad_time;
            std::vector<Double_t> pad_amp;
               
            if (1) {
               std::vector<Double_t> theta;
               std::vector<Double_t> radius;
               std::vector<Double_t> etheta;
               std::vector<Double_t> eradius;
               
               std::vector<Double_t> aw_theta;
               std::vector<Double_t> aw_radius;
               std::vector<Double_t> aw_etheta;
               std::vector<Double_t> aw_eradius;
               
               std::vector<Double_t> pads_theta;
               std::vector<Double_t> pads_radius;
               std::vector<Double_t> pads_etheta;
               std::vector<Double_t> pads_eradius;
               
               bool preamp_hits[16];
               for (unsigned j=0; j<16; j++) {
                  preamp_hits[j] = 0;
               }

               std::vector<Double_t> preamp_theta;
               std::vector<Double_t> preamp_radius;
               std::vector<Double_t> preamp_etheta;
               std::vector<Double_t> preamp_eradius;

               std::vector<Double_t> bars_theta;
               std::vector<Double_t> bars_radius;
               std::vector<Double_t> bars_etheta;
               std::vector<Double_t> bars_eradius;
               
               double rmin = 0.6;
               double rmax = 1.0;
               double r_aw = 1.1;
               double r_preamp = 1.15;
               double r_pads = 1.2;
               double r_bars = 1.3;
               
               if (0) {
                  for (int i=0; i<8; i++) {
                     theta.push_back((i+1)*(TMath::Pi()/4.));
                     radius.push_back((i+1)*0.05);
                     etheta.push_back(TMath::Pi()/8.);
                     eradius.push_back(0.05);
                  }
               }
               
               if (eawh) {
                  for (unsigned j=0; j<eawh->fAwHits.size(); j++) {
                     int iwire = eawh->fAwHits[j].wire;
                     double time = eawh->fAwHits[j].time;
                     double amp = eawh->fAwHits[j].amp;

                     if (iwire < 0)
                        continue;

                     hh->SetBinContent(1+iwire, time);
                     ha->SetBinContent(1+iwire, amp);
                     
                     int num_wires = 256;
                     
                     int itbwire = iwire%num_wires;
                     int itb = iwire/num_wires;
                     
                     int ipreamp = itbwire/16;
                     
                     double dist = (time - 1000.0)/4000.0;
                     if (dist < 0)
                        dist = 0;
                     if (dist > 1)
                        dist = 1;
                     
                     double t = (itbwire-8.0)/(1.0*num_wires)*(2.0*TMath::Pi());
                     double r = rmax-dist*(rmax-rmin);
                     
                     printf("aw hit %3d, wire %3d, tb %d, preamp %2d, iwire %3d, time %f, amp %f, theta %f (%f), radius %f\n", j, iwire, itb, ipreamp, itbwire, time, amp, t, t/TMath::Pi(), r);

                     preamp_hits[ipreamp] = true;
                     
                     //theta.push_back(t+0.5*TMath::Pi());
                     //radius.push_back(r);
                     //etheta.push_back(2.0*TMath::Pi()/num_wires);
                     //eradius.push_back(0.05);

                     aw_theta.push_back(t+0.5*TMath::Pi());
                     aw_radius.push_back(r_aw);
                     aw_etheta.push_back(2.0*TMath::Pi()/num_wires);
                     aw_eradius.push_back(0.05);

                     aw_theta.push_back(t+0.5*TMath::Pi());
                     aw_radius.push_back(r);
                     aw_etheta.push_back(2.0*TMath::Pi()/num_wires);
                     aw_eradius.push_back(0.05);
                  }
               }

               if (eph) {
                  for (unsigned i=0; i<eph->fPadHits.size(); i++) {
                     int imodule = eph->fPadHits[i].imodule;
                     double time = eph->fPadHits[i].time_ns;
                     double amp = eph->fPadHits[i].amp;
                     int col = eph->fPadHits[i].tpc_col;
                     int row = eph->fPadHits[i].tpc_row;

                     if (eph->fPadHits.size() < 5000) {
                        printf("pad hit %d: pwb%02d, col %d, row %d, time %f, amp %f\n", i, imodule, col, row, time, amp);
                     }

                     if (col < 0 || row < 0) {
                        continue;
                     }
                     
                     pad_col.push_back(col);
                     pad_row.push_back(row);
                     pad_time.push_back(time);
                     pad_amp.push_back(amp);
                     
                     assert(col >= 0 && col < NUM_PC);
                     
                     zpad_col[col].push_back(col);
                     zpad_row[col].push_back(row);
                     zpad_time[col].push_back(time);
                     zpad_amp[col].push_back(amp);
                     
                     //double dist = -0.2;
                     
                     double t = ((col+0.5)/(1.0*NUM_PC))*(2.0*TMath::Pi());
                     //double r = rmax-dist*(rmax-rmin);
                     
                     //printf("hit %d, wire %d, tb %d, iwire %d, t %f (%f), r %f\n", j, eawh->fAwHits[j].wire, itb, iwire, t, t/TMath::Pi(), r);
                     
                     //theta.push_back(t+0.5*TMath::Pi());
                     //radius.push_back(r);
                     //etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                     //eradius.push_back(0.05);

                     pads_theta.push_back(t+0.5*TMath::Pi());
                     pads_radius.push_back(r_pads);
                     pads_etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                     pads_eradius.push_back(0.05);

                     if (1) {
                        double dist = (time - 1000.0)/4000.0;
                        if (dist < 0)
                           dist = 0;
                        if (dist > 1)
                           dist = 1;
                        double r = rmax-dist*(rmax-rmin);
                        
                        //theta.push_back(t+0.5*TMath::Pi());
                        //radius.push_back(r);
                        //etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                        //eradius.push_back(0.05);

                        pads_theta.push_back(t+0.5*TMath::Pi());
                        pads_radius.push_back(r);
                        pads_etheta.push_back(2.0*TMath::Pi()/NUM_PC/2.0);
                        pads_eradius.push_back(0.05);
                     }
                  }
               }
               
               if (eba) {
                  for (unsigned j=0; j<eba->fBscAdcHits.size(); j++) {
                     //int adc_module = eba->fBscAdcHits[j].adc_module;
                     //int adc_chan = eba->fBscAdcHits[j].adc_chan;
                     //int preamp = eawh->fAwHits[j].preamp_pos;
                     int bar = eba->fBscAdcHits[j].bar;
                     double time = eba->fBscAdcHits[j].time;
                     double amp = eba->fBscAdcHits[j].amp;

                     double t = ((bar%64)/(1.0*64))*(2.0*TMath::Pi());
                     double r = r_bars;

                     printf("bsc adc hit %3d, bar %3d, tb %d, bar %2d, time %f, amp %f, theta %f (%f), radius %f\n", j, bar, bar/64, bar%64, time, amp, t, t/TMath::Pi(), r);
                     
                     //if (bsc64_bus) {
                     //   for (uint64_t i=0; i<64; i++) {
                     //      if (bsc64_bus & (((uint64_t)1)<<i)) {
                     //         h_bsc64_vs_bsc_adc_hits->Fill(i, bar);
                     //      }
                     //   }
                     //}

                     bars_theta.push_back(t+(0.5-0.1-1+1.0/8.0)*TMath::Pi());
                     bars_radius.push_back(r);
                     bars_etheta.push_back(2.0*TMath::Pi()/64/2.0);
                     bars_eradius.push_back(0.05);
                  }
               }

               if (1) {
                  printf("preamp hits:");
                  for (unsigned ipreamp=0; ipreamp<16; ipreamp++) {
                     if (preamp_hits[ipreamp]) {
                        printf(" %d", ipreamp);
                        double t = ((ipreamp)/(1.0*16))*(2.0*TMath::Pi());
                        preamp_theta.push_back(t+0.5*TMath::Pi());
                        preamp_radius.push_back(r_preamp);
                        preamp_etheta.push_back(16*2.0*TMath::Pi()/256.0/2.0);
                        preamp_eradius.push_back(0.0);
                     }
                  }
                  printf("\n");
               }

               if (1) {
                  printf("aw16_prompt_bits: 0x%04x: link hits: ", aw16_prompt);
                  for (int i=0; i<16; i++) {
                     if (aw16_prompt & (1<<i)) {
                        printf(" %d", i);
                     }
                  }
                  printf("\n");
               }

               if (1) {
                  printf("(out of date!) trig_bitmap: 0x%08x: bits: ", trig_bitmap);
                  if (trig_bitmap & (1<<0)) printf("adc16_grand_or ");
                  if (trig_bitmap & (1<<1)) printf("adc32_grand_or ");
                  if (trig_bitmap & (1<<2)) printf("adc_grand_or ");
                  if (trig_bitmap & (1<<3)) printf("esata_nim_grand_or ");
                  if (trig_bitmap & (1<<7)) printf("MLU ");
                  printf("\n");
               }

               if (1) {
                  theta.push_back(0+0.5*TMath::Pi());
                  radius.push_back(0);
                  etheta.push_back(TMath::Pi()/8.);
                  eradius.push_back(0.10);
               }
               
               if (1) {
                  theta.push_back(0+0.5*TMath::Pi());
                  radius.push_back(rmin);
                  etheta.push_back(2.0*TMath::Pi());
                  eradius.push_back(0.10);
               }
               
               if (1) {
                  theta.push_back(2.5*TMath::Pi());
                  radius.push_back(rmax);
                  etheta.push_back(2.0*TMath::Pi());
                  eradius.push_back(0.10);
               }
               
               fC->cd(3);
               TGraphPolar * grP1 = new TGraphPolar(theta.size(), theta.data(), radius.data(), etheta.data(), eradius.data());
               grP1->SetTitle("TPC end view from T side, wire 0 is at pi/2");
               grP1->SetMarkerStyle(20);
               grP1->SetMarkerSize(0.75);
               grP1->SetMarkerColor(4);
               grP1->SetLineColor(2);
               grP1->SetLineWidth(3);
               grP1->SetMinPolar(0);
               grP1->SetMaxPolar(2.0*TMath::Pi());
               grP1->SetMinRadial(0);
               grP1->SetMaxRadial(1.0);
               //grP1->SetMinimum(0);
               //grP1->SetMaximum(1);
               //grP1->Draw("PRE");
               grP1->Draw("PE");
               // Update, otherwise GetPolargram returns 0
               gPad->Update();
               grP1->GetPolargram()->SetToRadian();

               if (bars_theta.size() > 0) {
                  TGraphPolar* gbars = new TGraphPolar(bars_theta.size(), bars_theta.data(), bars_radius.data(), bars_etheta.data(), bars_eradius.data());
                  gbars->SetMarkerStyle(20);
                  gbars->SetMarkerSize(0.75);
                  gbars->SetMarkerColor(7);
                  gbars->SetLineColor(4);
                  gbars->SetLineWidth(3);
                  gbars->Draw("PE");
               }

               if (pads_theta.size() > 0) {
                  TGraphPolar* gpads = new TGraphPolar(pads_theta.size(), pads_theta.data(), pads_radius.data(), pads_etheta.data(), pads_eradius.data());
                  gpads->SetMarkerStyle(20);
                  gpads->SetMarkerSize(0.75);
                  gpads->SetMarkerColor(7);
                  gpads->SetLineColor(3);
                  gpads->SetLineWidth(3);
                  gpads->Draw("PE");
               }

               if (aw_theta.size() > 0) {
                  TGraphPolar* gaw = new TGraphPolar(aw_theta.size(), aw_theta.data(), aw_radius.data(), aw_etheta.data(), aw_eradius.data());
                  gaw->SetMarkerStyle(20);
                  gaw->SetMarkerSize(0.75);
                  gaw->SetMarkerColor(4);
                  gaw->SetLineColor(2);
                  gaw->SetLineWidth(3);
                  gaw->Draw("PE");
               }

               if (preamp_theta.size() > 0) {
                  TGraphPolar* gpreamp = new TGraphPolar(preamp_theta.size(), preamp_theta.data(), preamp_radius.data(), preamp_etheta.data(), preamp_eradius.data());
                  gpreamp->SetMarkerStyle(20);
                  gpreamp->SetMarkerSize(0.75);
                  gpreamp->SetMarkerColor(4);
                  gpreamp->SetLineColor(4);
                  gpreamp->SetLineWidth(3);
                  gpreamp->Draw("PE");
               }
            }

            TVirtualPad *p_pad_row = fC->cd(4);
            p_pad_row->Divide(1, 4);

            int zpad_colour[NUM_PC];
            int zpad_side[NUM_PC];

            if (1) {
               int col = 1;
               for (int i=0; i<NUM_PC; i++) {
                  zpad_side[i] = -1;
                  if (zpad_col[i].size() == 0) {
                     zpad_colour[i] = 0;
                  } else {
                     zpad_colour[i] = col;
                     col++;
                     if (col == 5)
                        col++;
                     if (col > 8) {
                        col = 1;
                     }
                  }
               }
            }


            if (1) { // split into two groups based on a gap
               int ifirst = -1;
               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  ifirst = i;
                  break;
               }
               
               zpad_side[ifirst] = 1;

               int igap = 0;
               for (int j=0; j<NUM_PC; j++) {
                  int i = (ifirst+j)%NUM_PC;
                  if (zpad_col[i].size() == 0) {
                     igap++;
                     //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                     if (igap>1)
                        break;
                     continue;
                  }
                  igap=0;
                  //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                  zpad_side[i] = 1;
               }

               igap = 0;
               for (int j=0; j<NUM_PC; j++) {
                  int i = (ifirst-j)%NUM_PC;
                  if (i<0)
                     i+=NUM_PC;
                  if (zpad_col[i].size() == 0) {
                     igap++;
                     //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                     if (igap>1)
                        break;
                     continue;
                  }
                  igap=0;
                  //printf("ifirst %d, j %d, i %d, igap %d\n", ifirst, j, i, igap);
                  zpad_side[i] = 1;
               }
            }

            bool split_by_max_drift = false;

            if (1) {
               int count_plus = 0;
               int count_minus = 0;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  if (zpad_side[i]>0) count_plus++;
                  if (zpad_side[i]<0) count_minus++;
               }
               printf("counts: plus %d, minus %d\n", count_plus, count_minus);
               if (count_plus>0 && count_minus>0)
                  split_by_max_drift = false;
               else
                  split_by_max_drift = true;
            }

            if (split_by_max_drift) { // split into two groups separated by biggest drift time
               double max_time = -1;
               int max_pc = -1;
               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  for (unsigned j=0; j<zpad_col[i].size(); j++) {
                     if (zpad_time[i][j] > max_time) {
                        //printf("max_pc %d->%d, time %f->%f\n", max_pc, i, max_time, zpad_time[i][j]);
                        max_time = zpad_time[i][j];
                        max_pc = i;
                     }
                  }
               }
               if (max_pc >= 0) {
                  for (int i=0; i<NUM_PC; i++) {
                     if (i<max_pc) {
                        zpad_side[i] = 1;
                     } else {
                        zpad_side[i] = -1;
                     }
                  }
               }
            }

            if (pad_col.size() > 0) {
               p_pad->cd(1);
               TGraph* hpct = new TGraph(pad_col.size(), pad_col.data(), pad_time.data());
               hpct->SetTitle("Pads hit time per pad column; pad column; hit time, ns");
               hpct->SetMarkerStyle(20);
               hpct->SetMarkerSize(0.75);
               hpct->SetMarkerColor(4);
               hpct->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
               hpct->SetMinimum(0);
               hpct->SetMaximum(MAX_TIME);
               //hpct->Draw("AC*");
               hpct->Draw("A*");

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;

                  TGraph* g = new TGraph(zpad_col[i].size(), zpad_col[i].data(), zpad_time[i].data());
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_TIME);
                  //hpct->Draw("AC*");
                  g->Draw("*");
               }
            }

            if (pad_col.size() > 0) {
               p_pad->cd(2);
               TGraph* hpca = new TGraph(pad_col.size(), pad_col.data(), pad_amp.data());
               hpca->SetTitle("Pads hit amplitude per pad column; pad column; hit amplitude, ADC counts");
               hpca->SetMarkerStyle(20);
               hpca->SetMarkerSize(0.75);
               hpca->SetMarkerColor(4);
               hpca->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
               hpca->SetMinimum(0);
               hpca->SetMaximum(MAX_PAD_AMP);
               //hpct->Draw("AC*");
               hpca->Draw("A*");

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;

                  TGraph* g = new TGraph(zpad_col[i].size(), zpad_col[i].data(), zpad_amp[i].data());
                  g->SetTitle("Pads hit amplitude per pad row; pad row; hit amplitude, ADC counts");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PC-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_PAD_AMP);
                  //hpct->Draw("AC*");
                  g->Draw("*");
               }
            }

            if (1) {
               p_pad_row->cd(1);
               
               //TGraph* hprt = new TGraph();
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(0);
               //hprt->SetMaximum(MAX_PAD_AMP);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_col[i].size() == 0)
                     continue;
                  if (zpad_side[i] != 1)
                     continue;

                  //printf("111 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_amp[i].data());
                  g->SetTitle("Pads hit time per pad row; pad row; hit time, ns");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_PAD_AMP);
                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            if (1) {
               p_pad_row->cd(2);

               //TGraph* hprt = new TGraph(pad_row.size(), pad_row.data(), pad_time.data());
               //hprt->SetTitle("XXX; aaa; bbb");
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(-MAX_TIME);
               //hprt->SetMaximum(0);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_row[i].size() == 0)
                     continue;
                  if (zpad_side[i] != 1)
                     continue;

                  for (unsigned j=0; j<zpad_row[i].size(); j++) {
                     zpad_time[i][j] = -zpad_time[i][j];
                  }

                  //printf("222 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_time[i].data());
                  g->SetTitle("Pads hit time per pad row; pad row; hit time, ns");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(-PLOT_MAX_TIME);
                  g->SetMaximum(-PLOT_MIN_TIME);

                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            if (1) {
               p_pad_row->cd(3);
               
               //TGraph* hprt = new TGraph(pad_row.size(), pad_row.data(), pad_time.data());
               //hprt->SetTitle("XXX; aaa; bbb");
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(0);
               //hprt->SetMaximum(MAX_TIME);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_row[i].size() == 0)
                     continue;
                  if (zpad_side[i] != -1)
                     continue;

                  //printf("333 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_time[i].data());
                  g->SetTitle("Pads hit time per pad row; pad row; hit time, ns");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(PLOT_MIN_TIME);
                  g->SetMaximum(PLOT_MAX_TIME);
                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            if (1) {
               p_pad_row->cd(4);

               //TGraph* hprt = new TGraph(pad_row.size(), pad_row.data(), pad_amp.data());
               //hprt->SetMarkerStyle(20);
               //hprt->SetMarkerSize(0.75);
               //hprt->SetMarkerColor(4);
               //hprt->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
               //hprt->SetMinimum(0);
               //hprt->SetMaximum(MAX_PAD_AMP);
               //hprt->Draw("AC*");
               //hprt->Draw("A*");

               bool first = true;

               for (int i=0; i<NUM_PC; i++) {
                  if (zpad_row[i].size() == 0)
                     continue;
                  if (zpad_side[i] != -1)
                     continue;

                  //printf("444 col %d, size %d\n", i, zpad_row[i].size());
                  TGraph* g = new TGraph(zpad_row[i].size(), zpad_row[i].data(), zpad_amp[i].data());
                  g->SetTitle("Pads hit amplitude per pad row; pad row; hit amplitude, ADC counts");
                  g->SetMarkerStyle(20);
                  g->SetMarkerSize(0.75);
                  g->SetMarkerColor(zpad_colour[i]);
                  g->GetXaxis()->SetLimits(-0.5, NUM_PR-0.5);
                  g->SetMinimum(0);
                  g->SetMaximum(MAX_PAD_AMP);
                  if (first) {
                     g->Draw("A*");
                     first = false;
                  } else {
                     g->Draw("*");
                  }
               }
            }

            fC->Modified();
            fC->Draw();
            fC->Update();
            
            save_gpad->cd();
         }

         if (age->trig && age->trig->udpData.size() > 0) {
            for (unsigned i=0; i<age->trig->udpData.size(); i++) {
               printf("ATAT[%d]: 0x%08x (%d)\n", i, age->trig->udpData[i], age->trig->udpData[i]);
            }
         }

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

      //*flags |= TAFlag_DISPLAY;
#endif

      return flow;
   }

   void AnalyzeSpecialEvent(TARunInfo* runinfo, TMEvent* event)
   {
      if (fTrace)
         printf("DisplayModule::AnalyzeSpecialEvent, run %d, event serno %d, id 0x%04x, data size %d\n", runinfo->fRunNo, event->serial_number, (int)event->event_id, event->data_size);
   }
};

class DisplayModuleFactory: public TAFactory
{
public:
   bool fDoPrint = false;
   bool fDoEventDisplay = true;

public:
   void Usage()
   {
      printf("DisplayModuleFactory::Usage:\n");
      printf("--print ## be verbose\n");
      printf("--nodisplay ## disable the event display\n");
   }

   void Init(const std::vector<std::string> &args)
   {
      printf("DisplayModuleFactory::Init!\n");
      
      //fDoPads = true;
      //fPlotPad = -1;
      //fPlotPadCanvas = NULL;
      
      for (unsigned i=0; i<args.size(); i++) {
         //if (args[i] == "--nopads")
         //   fDoPads = false;
         //if (args[i] == "--plot1")
         //   fPlotPad = atoi(args[i+1].c_str());

         if (args[i] == "--print")
            fDoPrint = true;
         if (args[i] == "--nodisplay")
            fDoEventDisplay = false;
      }
   }

   void Finish()
   {
      printf("DisplayModuleFactory::Finish!\n");
   }

   TARunObject* NewRunObject(TARunInfo* runinfo)
   {
      printf("DisplayModule::NewRunObject, run %d, file %s\n", runinfo->fRunNo, runinfo->fFileName.c_str());
      return new DisplayModule(runinfo, fDoPrint, fDoEventDisplay);
   }
};

static TARegister tar(new DisplayModuleFactory);

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
