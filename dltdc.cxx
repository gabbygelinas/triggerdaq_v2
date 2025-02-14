
#undef NDEBUG // this program requires a working assert()

#include "dltdc.h"

#include <stdio.h>
#include <stdlib.h> // drand48()
#include <assert.h> // assert()
#include <math.h> // sqrt()

void DlTdcHit::Print() const
{
   printf("data 0x%08x 0x%08x, ch %2d lete %d%d, ts0 %d phase %3d, %5.1f ns, time: %.0f %.9f %.9f sec", data_hi, data_lo, ch, le, te, ts0, phase, fine_ns, coarse_epoch, coarse_sec, time_sec);
}

DlTdcFineCalib1::DlTdcFineCalib1() // ctor
{
   Resize(80);
   Reset();
}

DlTdcFineCalib1::~DlTdcFineCalib1() // dtor
{

}

void DlTdcFineCalib1::Resize(int nbins)
{
   fHistogram.resize(nbins);
   fBinWidthNs.resize(nbins);
   fBinTimeNs.resize(nbins);
}

void DlTdcFineCalib1::Reset()
{
   fMaxPhase = 0;
   fHits = 0;
   for (size_t i=0; i<fHistogram.size(); i++) {
      fHistogram[i] = 0;
      fBinWidthNs[i] = fTotalNs;
      fBinTimeNs[i] = 0;
   }
}

void DlTdcFineCalib1::AddHit(int phase, double fine_ns)
{
   if (phase == 0)
      return;

   assert(phase > 0);
   assert(phase < (int)fHistogram.size());
   
   fHits += 1;
   fHistogram[phase] += 1;

   if (phase > fMaxPhase) {
      //printf("max phase %d -> %d\n", fMaxPhase, phase);
      fMaxPhase = phase;
      Update();
      //Print();
   }
}

void DlTdcFineCalib1::Update()
{
   double sum = 0;
   for (size_t i=0; i<fHistogram.size(); i++) {
      sum += fHistogram[i];
   }

   //printf("sum %.0f\n", sum);

   if (sum < 10)
      return;

   double cumulative = 0;

   fBinMinNs = 99999;
   fBinMaxNs = 0;

   for (size_t i=0; i<fHistogram.size(); i++) {
      cumulative += fHistogram[i];
      fBinTimeNs[i+1] = cumulative/sum * fTotalNs;
      fBinWidthNs[i] = fHistogram[i]/sum * fTotalNs;

      if (fBinWidthNs[i] > fBinMaxNs)
         fBinMaxNs = fBinWidthNs[i];

      if (fHistogram[i] > 0) {
         if (fBinWidthNs[i] < fBinMinNs)
            fBinMinNs = fBinWidthNs[i];
      }
   }

   //Print();
}

void DlTdcFineCalib1::Print() const
{
   printf("hist: ");
   
   for (size_t i=0; i<fHistogram.size(); i++) {
      printf(" %3.0f", fHistogram[i]);
   }
   
   printf("\n");
   
   printf("binw: ");
   
   for (size_t i=0; i<fBinWidthNs.size(); i++) {
      if (fBinWidthNs[i] > 0) {
         printf(" %3.1f", fBinWidthNs[i]);
      } else {
         printf("   -");
      }
   }
   
   printf(", min %5.3f, max %3.1f ns\n", fBinMinNs, fBinMaxNs);
   
   printf("time: ");
   
   for (size_t i=0; i<fBinTimeNs.size(); i++) {
      printf(" %3.1f", fBinTimeNs[i]);
   }
   
   printf("\n");
}

std::string DlTdcFineCalib1::toJson() const
{
  char buf[256];

  std::string s;
  s += "{\n";
  sprintf(buf, "  \"TotalNs\":%.1f", fTotalNs);
  s += buf;
  s += ",\n";
  sprintf(buf, "  \"Hits\":%zu", fHits);
  s += buf;
  s += ",\n";
  sprintf(buf, "  \"MaxPhase\":%d", fMaxPhase);
  s += buf;
  s += ",\n";
  sprintf(buf, "  \"BinMinNs\":%.3f", fBinMinNs);
  s += buf;
  s += ",\n";
  sprintf(buf, "  \"BinMaxNs\":%.3f", fBinMaxNs);
  s += buf;
  s += ",\n";
  sprintf(buf, "  \"Bins\":%zu", fHistogram.size());
  s += buf;
  s += ",\n";
  s += "  \"Histogram\":[";
  for (size_t i=0; i<fHistogram.size(); i++) {
    if (i!=0)
      s += ",";
    sprintf(buf,"%.0f",fHistogram[i]);
    s += buf;
  }
  s += "],\n";
  s += "  \"BinWidthNs\":[";
  for (size_t i=0; i<fBinWidthNs.size(); i++) {
    if (i!=0)
      s += ",";
    sprintf(buf,"%.3f",fBinWidthNs[i]);
    s += buf;
  }
  s += "],\n";
  s += "  \"BinTimeNs\":[";
  for (size_t i=0; i<fBinTimeNs.size(); i++) {
    if (i!=0)
      s += ",";
    sprintf(buf,"%.3f",fBinTimeNs[i]);
    s += buf;
  }
  s += "]\n";
  s += "  }";
  return s;
}

double DlTdcFineCalib1::GetTime(int phase)
{
   if (phase == 0)
      return 0;
   //printf("phase %d\n", phase);
   assert(phase > 0);
   assert(phase < (int)fBinTimeNs.size());
   double time  = fBinTimeNs[phase];
   double width = fBinWidthNs[phase];
   double random = drand48();
   return time + width*random;
}

DlTdcFineCalib::DlTdcFineCalib() // ctor
{
   Reset();
}

DlTdcFineCalib::~DlTdcFineCalib() // dtor
{

}

void DlTdcFineCalib::Reset()
{
   lepos.Reset();
   leneg.Reset();
   tepos.Reset();
   teneg.Reset();
}

void DlTdcFineCalib::Print() const
{
   printf("min bin: le: %.3f %.3f, te: %.3f %.3f ns, max bin: le: %.3f %.3f, te: %.3f %.3f ns, max phase: le: %2d %2d, te: %2d %2d\n",
          lepos.fBinMinNs,
          leneg.fBinMinNs,
          tepos.fBinMinNs,
          teneg.fBinMinNs,
          lepos.fBinMaxNs,
          leneg.fBinMaxNs,
          tepos.fBinMaxNs,
          teneg.fBinMaxNs,
          lepos.fMaxPhase,
          leneg.fMaxPhase,
          tepos.fMaxPhase,
          teneg.fMaxPhase);
}

#if 0
static void SaveToFile1(FILE* fp, const DlTdcFineCalib1& c)
{
   fprintf(fp, "total %.3f, max_phase %d, bins %zu\n", c.fTotalNs, c.fMaxPhase, c.fHistogram.size());

   //fprintf(fp, "h");
   for (size_t i=0; i<c.fHistogram.size(); i++) {
      if ((i%10)==0)
         fprintf(fp, "\n");
      fprintf(fp, " %6.0f", c.fHistogram[i]);
   }
   fprintf(fp, "\n");

   //fprintf(fp, "w");
   for (size_t i=0; i<c.fBinWidthNs.size(); i++) {
      if ((i%10)==0)
         fprintf(fp, "\n");
      fprintf(fp, " %6.3f", c.fBinWidthNs[i]);
   }
   fprintf(fp, "\n");

   //fprintf(fp, "t");
   for (size_t i=0; i<c.fBinTimeNs.size(); i++) {
      if ((i%10)==0)
         fprintf(fp, "\n");
      fprintf(fp, " %6.3f", c.fBinTimeNs[i]);
   }
   fprintf(fp, "\n");

   fprintf(fp, "\n");
}
#endif

#include "mjson.h"

static void LoadFromJson1(const MJsonNode*j, DlTdcFineCalib1& c)
{
   if (!j)
      return;

   //printf("XXX: %s\n", j->Stringify().c_str());

   size_t bins = j->FindObjectNode("Bins")->GetInt();
   c.fTotalNs  = j->FindObjectNode("TotalNs")->GetDouble();
   c.fHits     = j->FindObjectNode("Hits")->GetInt();
   c.fMaxPhase = j->FindObjectNode("MaxPhase")->GetInt();
   c.fBinMinNs = j->FindObjectNode("BinMinNs")->GetDouble();
   c.fBinMaxNs = j->FindObjectNode("BinMaxNs")->GetDouble();

   c.fHistogram.clear();
   const MJsonNodeVector* v = j->FindObjectNode("Histogram")->GetArray();
   for (size_t i=0; i<v->size(); i++) {
      c.fHistogram.push_back((*v)[i]->GetDouble());
   }
   assert(c.fHistogram.size() == bins);

   c.fBinWidthNs.clear();
   v = j->FindObjectNode("BinWidthNs")->GetArray();
   for (size_t i=0; i<v->size(); i++) {
      c.fBinWidthNs.push_back((*v)[i]->GetDouble());
   }
   assert(c.fBinWidthNs.size() == bins);

   c.fBinTimeNs.clear();
   v = j->FindObjectNode("BinTimeNs")->GetArray();
   for (size_t i=0; i<v->size(); i++) {
      c.fBinTimeNs.push_back((*v)[i]->GetDouble());
   }
   assert(c.fBinTimeNs.size() == bins);

   //c.Print();
   //exit(123);
}

std::string DlTdcFineCalib::toJson() const
{
  std::string s;
  s += "{\n";
  s += " \"lepos\": ";
  s += lepos.toJson();
  s += ",\n";
  s += " \"leneg\": ";
  s += leneg.toJson();
  s += ",\n";
  s += " \"tepos\": ";
  s += tepos.toJson();
  s += ",\n";
  s += " \"teneg\": ";
  s += teneg.toJson();
  s += "\n";
  s += "}";
  return s;
}


void DlTdcFineCalib::SaveToFile(const char* filename) const
{
   Print();

   FILE *fp = fopen(filename, "w");
   assert(fp);
   fprintf(fp, "%s\n", toJson().c_str());
   fclose(fp);
}

bool DlTdcFineCalib::LoadFromJson(const MJsonNode* j)
{
   //printf("DlTdcFineCalib::LoadFromJson: %s\n", j->Stringify().c_str());

   LoadFromJson1(j->FindObjectNode("lepos"), lepos);
   LoadFromJson1(j->FindObjectNode("leneg"), leneg);
   LoadFromJson1(j->FindObjectNode("tepos"), tepos);
   LoadFromJson1(j->FindObjectNode("teneg"), teneg);

   return true;
}

bool DlTdcFineCalib::LoadFromFile(const char* filename)
{
   FILE *fp = fopen(filename, "r");
   if (!fp)
     return false;
   assert(fp);

   printf("Loading DL-TDC calibrations from \"%s\"\n", filename);

   std::string json;
   while (1) {
      char buf[1024];
      size_t rd = fread(buf, 1, sizeof(buf)-1, fp);
      //printf("rd %zu\n", rd);
      if (rd == 0)
         break;
      buf[rd] = 0;
      json += buf;
   }

   fclose(fp);

   //printf("json: %s\n", json.c_str());

   MJsonNode* j = MJsonNode::Parse(json.c_str());
   //printf("read: %s\n", j->Stringify().c_str());

   LoadFromJson(j);

   delete j;

   return true;
}

void DlTdcFineCalib::AddHit(const DlTdcHit& h)
{
   if (h.le) {
      if (h.phase > 0)
         lepos.AddHit(h.phase, h.fine_ns);
      else
         leneg.AddHit(-h.phase, h.fine_ns);
   }

   if (h.te) {
      if (h.phase > 0)
         tepos.AddHit(h.phase, h.fine_ns);
      else
         teneg.AddHit(-h.phase, h.fine_ns);
   }
}

void DlTdcFineCalib::Update()
{
   lepos.Update();
   leneg.Update();
   tepos.Update();
   teneg.Update();
}



DlTdcUnpack::DlTdcUnpack(int nchan) // ctor
{
   fLastCoarse.resize(nchan);
   fEpoch.resize(nchan);
   fEpochHits.resize(nchan);
   fCalib.resize(nchan);
};

DlTdcUnpack::~DlTdcUnpack() // dtor
{

};

void DlTdcUnpack::Reset()
{
   size_t n = fEpoch.size();
   for (size_t i=0; i<n; i++) {
      fLastCoarse[i] = 0;
      fEpoch[i] = 0;
      fEpochHits[i] = 0;
   }
   fFirstTimeSec = 0;

   for (auto& c: fCalib) {
      c.Reset();
   }
}

void DlTdcUnpack::PrintBits32(uint32_t v)
{
   for (int i=0; i<32; i++) {
      if (v & (1<<(31-i)))
         printf("1");
      else
         printf("0");
      
      if ((i%4)==3)
         printf("-");
   }
}

uint32_t DlTdcUnpack::FixHoles(uint32_t v)
{
   if ((v & 0xF) == 2)
      v |= 1;

   //if ((v & 0x3) == 1)
   //   v &= ~1;

   //if ((v & 0xF) == 0xE)
   //   v |= 1;

   for (int j=1; j<31; j++)
      if ((v & (1<<(j-1))) && ((v & (1<<(j+1)))))
         v |= (1<<j);

   for (int j=1; j<31; j++)
      if (((v & (1<<(j-1)))==0) && (((v & (1<<(j+1))))==0))
         v &= ~(1<<j);

   return v;
}

#if 0
int findLast1(uint32_t v)
{
   int i=0;
   for (; i<32; i++)
      if ((v & (1<<i)) != 0)
         break;
   for (; i<32; i++)
      if ((v & (1<<i)) == 0)
         return i;

   return i;
}

int findBits(uint32_t v)
{
   if ((v & 1) == 0)
      {
         int i=0;
         for (; i<32; i++)
            if ((v & (1<<i)) != 0)
               break;

         int start = i;

         for (; i<32; i++)
            if ((v & (1<<i)) == 0)
               return (start+i)/2;
      }
   else
      {
         int i=0;
         for (; i<32; i++)
            if ((v & (1<<i)) == 0)
               break;

         int start = i;

         for (; i<32; i++)
            if ((v & (1<<i)) != 0)
               return -(start+i)/2;
      }

   return 0;
}
#endif

int DlTdcUnpack::FindEdge(uint32_t v)
{
   if ((v & 1) == 0) {
      int i=0;
      for (; i<32; i++)
         if ((v & (1<<i)) != 0)
            return -i;
   } else {
      int i=0;
      for (; i<32; i++)
         if ((v & (1<<i)) == 0)
            return i;
   }

   return 0;
}

int DlTdcUnpack::FindEdge10(uint32_t v)
{
   if ((v & 1) == 0) {
      int i=0;
      for (; i<32; i++)
         if ((v & (1<<i)) != 0)
            break;
      for (; i<32; i++)
         if ((v & (1<<i)) == 0)
            return i;
   } else {
      int i=0;
      for (; i<32; i++)
         if ((v & (1<<i)) == 0)
            return i;
   }

   return 0;
}

bool DlTdcUnpack::Unpack(DlTdcHit*h, uint32_t lo, uint32_t hi)
{
   h->Clear();
   h->data_lo = lo;
   h->data_hi = hi;

   //printf("Unpack 0x%08x 0x%08x\n", hi, lo);

   h->le = hi & 0x80000000;
   h->te = hi & 0x40000000;

   h->coarse = (hi-1) & 0x3FFFFFFF;

   h->ts0 = h->coarse & 1;

   int ch = (lo>>24)&0xFF;
   //int ch = (lo>>28)&0xF;
   h->ch = ch;

   int ph = (lo>>16)&0xFF;

   //if (ph == 63)
   //   ph = 9999;
   //else if (ph & 0x20)
   //   ph = -(ph&~0x20);

   if (ph == 0xFF) // encoder error
      ph = 0;
   else if (ph == 0x7F) // all-1
      ph = 0;
   else if (ph == 0x00) // all-0
      ph = 0;
   else if (ph & 0x80)
      ph = -(ph&~0x80);

   if (h->coarse < fLastCoarse[ch]) {
      double max_epoch = 0;
      for (size_t cc=0; cc<fEpoch.size(); cc++) {
         if (fEpoch[cc] > max_epoch) {
            max_epoch = fEpoch[cc];
         }
      }

      //printf("ch %d epoch increment %6.0f with %4.0f hits, max epoch %6.0f: 0x%08x -> 0x%08x\n", ch, fEpoch[ch], fEpochHits[ch], max_epoch, fLastCoarse[ch], h->coarse);

      if (fEpoch[ch] == max_epoch) {
         // we are the first to see epoch increment
         fEpoch[ch] += 1;
      } else if (fEpoch[ch] + 1 == max_epoch) {
         // we will increment epoch and it will be same as others
         fEpoch[ch] += 1;
      } else {
         // we missed an epoch update!
         printf("missed epoch: ch %d forced update from %.0f to %.0f, have %.0f hits\n", ch, fEpoch[ch], max_epoch, fEpochHits[ch]);
         fEpoch[ch] = max_epoch;
      }
      fEpochHits[ch] = 0;
   }

   fLastCoarse[ch] = h->coarse;
   fEpochHits[ch] += 1;
   h->coarse_epoch = fEpoch[ch];

   uint32_t sr1 = lo & 0xFFFF;

   if (sr1 & 0x8000)
      sr1 |= 0xFFFF0000;

   //uint32_t sr1 = FixHoles(sr);

#if 1
   // version for delay chain fed by timestamp[0]


   h->phase = 0; // FindEdge10(sr1);

   h->fine_ns = 0; // h->phase * fNsPerBit[ch]; // this is wrong, should use min and max

   bool calib = true;
   double bits = 40.0;

   h->coarse_sec = fClkPeriodNs*h->coarse*1e-9 + h->coarse_epoch*0x40000000*fClkPeriodNs*1e-9;
   
   if (fFirstTimeSec == 0)
      fFirstTimeSec = h->coarse_sec;
   
   h->coarse_sec -= fFirstTimeSec;
   
   if (h->coarse_sec < -1.0 && h->coarse_epoch == 0) {
      printf("ch %d: correct epoch0\n", ch);
      fEpoch[ch] = 1;
      h->coarse_epoch = 1;
      h->coarse_sec = fClkPeriodNs*h->coarse*1e-9 + h->coarse_epoch*0x40000000*fClkPeriodNs*1e-9;
      h->coarse_sec -= fFirstTimeSec;
   }

   h->phase = ph; // FindEdge(sr1);

   if (h->phase > 0) {
      if (calib) {
         if (h->le) {
            h->fine_ns = fCalib[ch].lepos.GetTime(h->phase);
            h->offset_ns = fCalib[ch].lepos.fOffsetNs;
         }
         if (h->te) {
            h->fine_ns = fCalib[ch].tepos.GetTime(h->phase);
            h->offset_ns = fCalib[ch].tepos.fOffsetNs;
            //printf("ch %d: ", ch); fCalib[ch].tepos.Print();
         }
      } else {
         h->fine_ns = fClkPeriodNs/bits*h->phase;
         h->offset_ns = 0;
      }

      if (h->ts0==0) {
         if (h->fine_ns > fClkPeriodNs/2.0)
            h->offset_ns -= fClkPeriodNs;
         else
            h->offset_ns += fClkPeriodNs;
      }
   } else {
      if (calib) {
         if (h->le) {
            h->fine_ns = fCalib[ch].leneg.GetTime(-h->phase);
            h->offset_ns = fCalib[ch].leneg.fOffsetNs;
         }
         if (h->te) {
            h->fine_ns = fCalib[ch].teneg.GetTime(-h->phase);
            h->offset_ns = fCalib[ch].teneg.fOffsetNs;
         }
      } else {
         h->fine_ns = -fClkPeriodNs/bits*h->phase;
         h->offset_ns = 0;
      }

      if (h->ts0==1) {
         if (h->fine_ns > fClkPeriodNs/2.0)
            h->offset_ns -= fClkPeriodNs;
         else
            h->offset_ns += fClkPeriodNs;
      }
   }

#if 0
   if (h->ts0 == 0) {
      if (h->phase > 0) {
         if (calib) {
            if (h->le) {
               h->fine_ns = fCalib[ch].lepos0.GetTime(h->phase) - fClkPeriodNs;
               h->offset_ns = fCalib[ch].lepos0.fOffsetNs;
            }
            if (h->te) {
               h->fine_ns = fCalib[ch].tepos0.GetTime(h->phase) - fClkPeriodNs;
               h->offset_ns = fCalib[ch].tepos0.fOffsetNs;
               //printf("ch %d: ", ch); fCalib[ch].tepos0.Print();
            }
         } else {
            h->fine_ns = fClkPeriodNs/bits*h->phase - fClkPeriodNs;
            h->offset_ns = 0;
         }
      } else {
         if (calib) {
            if (h->le) {
               h->fine_ns = fCalib[ch].leneg0.GetTime(-h->phase);
               h->offset_ns = fCalib[ch].leneg0.fOffsetNs;
            }
            if (h->te) {
               h->fine_ns = fCalib[ch].teneg0.GetTime(-h->phase);
               h->offset_ns = fCalib[ch].teneg0.fOffsetNs;
            }
         } else {
            h->fine_ns = -fClkPeriodNs/bits*h->phase;
            h->offset_ns = 0;
         }
      }
   } else {
      if (h->phase > 0) {
         if (calib) {
            if (h->le) {
               h->fine_ns = fCalib[ch].lepos1.GetTime(h->phase);
               h->offset_ns = fCalib[ch].lepos1.fOffsetNs;
            }
            if (h->te) {
               h->fine_ns = fCalib[ch].tepos1.GetTime(h->phase);
               h->offset_ns = fCalib[ch].tepos1.fOffsetNs;
            }
         } else {
            h->fine_ns = fClkPeriodNs/bits*h->phase;
            h->offset_ns = 0;
         }
      } else {
         if (calib) {
            if (h->le) {
               h->fine_ns = fCalib[ch].leneg1.GetTime(-h->phase) - fClkPeriodNs;
               h->offset_ns = fCalib[ch].leneg1.fOffsetNs;
            }
            if (h->te) {
               h->fine_ns = fCalib[ch].teneg1.GetTime(-h->phase) - fClkPeriodNs;
               h->offset_ns = fCalib[ch].teneg1.fOffsetNs;
            }
         } else {
            h->fine_ns = -fClkPeriodNs/bits*h->phase - fClkPeriodNs;
            h->offset_ns = 0;
         }
      }
   }

   if (h->fine_ns < -2.0) {
      h->fine_ns += 2.0*fClkPeriodNs;
   }
#endif

   //if (h->fine_ns < -4.0) {
   //   h->Print();
   //   abort();
   //}
   
   //printf("phase %2d %2d\n", ph, h->phase);

   //if (ph == h->phase*2) {
   //   // good
   //} else if (ph == h->phase*2-1) {
   //   // good
   //} else {
   //   printf("BAD phase %2d %2d\n", ph, h->phase);
   //}

   //if (h->phase_ns < 0)
   //   h->phase_ns = -h->phase_ns + fClkPeriodNs/2.0;

   //if (h->le) {
   //   if (h->phase > 24) {
   //      h->fine_ns -= fClkPeriodNs;
   //   }
   //}

   //if (h->te) {
   //   if (h->phase > 25) {
   //      h->fine_ns -= fClkPeriodNs;
   //   }
   //}

   //int p1 = FindEdge(sr1);

   //double f1 = 0;

   //if (p1 > 0) {
   //   //f1 = (p1-1) * fNsPerBit[ch];
   //   if (h->le)
   //      f1 = fCalib[ch].lepos.GetTime(p1);
   //   if (h->te)
   //      f1 = fCalib[ch].tepos.GetTime(p1);
   //}

   //if (p1 < 0) {
   //   //f1 = - p1 * fNsPerBit[ch] + fClkPeriodNs/2.0;
   //   if (h->le)
   //      f1 = fCalib[ch].leneg.GetTime(-p1) + fClkPeriodNs/2.0;
   //   if (h->te)
   //      f1 = fCalib[ch].teneg.GetTime(-p1) + fClkPeriodNs/2.0;
   //}

   //printf("phase %2d fine %5.1f ns, phase %3d, fine %5.1f\n", h->phase, h->fine_ns, p1, f1);

   //h->phase = p1;
   //h->fine_ns = f1;

#endif

#if 0
   // version for delay chain fed by 100 MHz clock

   h->coarse_sec = fClkPeriodNs*h->coarse*1e-9 + h->coarse_epoch*0x40000000*fClkPeriodNs*1e-9;
   
   if (fFirstTimeSec == 0)
      fFirstTimeSec = h->coarse_sec;

   h->coarse_sec -= fFirstTimeSec;

   h->phase = 0; // FindEdge10(sr1);

   h->fine_ns = 0; // h->phase * fNsPerBit[ch]; // this is wrong, should use min and max

   //if (h->phase_ns < 0)
   //   h->phase_ns = -h->phase_ns + fClkPeriodNs/2.0;

   //if (h->le) {
   //   if (h->phase > 24) {
   //      h->fine_ns -= fClkPeriodNs;
   //   }
   //}

   //if (h->te) {
   //   if (h->phase > 25) {
   //      h->fine_ns -= fClkPeriodNs;
   //   }
   //}

   int p1 = FindEdge(sr1);

   double f1 = 0;

   if (p1 > 0) {
      //f1 = (p1-1) * fNsPerBit[ch];
      if (h->le)
         f1 = fCalib[ch].lepos.GetTime(p1);
      if (h->te)
         f1 = fCalib[ch].tepos.GetTime(p1);
   }

   if (p1 < 0) {
      //f1 = - p1 * fNsPerBit[ch] + fClkPeriodNs/2.0;
      if (h->le)
         f1 = fCalib[ch].leneg.GetTime(-p1) + fClkPeriodNs/2.0;
      if (h->te)
         f1 = fCalib[ch].teneg.GetTime(-p1) + fClkPeriodNs/2.0;
   }

   //printf("phase %2d fine %5.1f ns, phase %3d, fine %5.1f\n", h->phase, h->fine_ns, p1, f1);

   h->phase = p1;
   h->fine_ns = f1;

#endif

#if 0
   h->phase = FindEdge(sr1);

   if (h->coarse & 1) {
      if (h->phase > 0) {
         h->fine_ns = -h->phase * fNsPerBit + fClkPeriodNs;
         h->phase_ns = h->fine_ns;
         //h->coarse &= ~1;
      } else {
         h->fine_ns = h->phase * fNsPerBit;
         h->phase_ns = h->fine_ns;
         //h->coarse &= ~1;
      }
   } else {
      if (h->phase > 0) {
         h->fine_ns = h->phase * fNsPerBit;
         h->phase_ns = h->fine_ns;
      } else {
         h->fine_ns = -h->phase * fNsPerBit;
         h->phase_ns = h->fine_ns;
      }
   }
#endif

   //h->phase_ns = h->phase * fNsPerBit; // this is wrong, should use min and max

   //if (h->phase_ns < 0)
   //   h->phase_ns = -h->phase_ns + fClkPeriodNs/2.0;

   //if (h->le) {
   //   if (h->phase > 24) {
   //      h->phase_ns -= fClkPeriodNs;
   //   }
   //}

   //if (h->te) {
   //   if (h->phase > 25) {
   //      h->phase_ns -= fClkPeriodNs;
   //   }
   //}

#if 0
   // version from TDC4 block

   if ((h->coarse & 1) == 0) {
      if (h->phase_ns > fClkPeriodNs) {
         printf("subtract %f ns\n", 2*fClkPeriodNs);
         h->phase_ns -= 2*fClkPeriodNs;
      }
   }

   if (h->coarse & 1) {
      printf("subtrack %f ns coarse\n", fClkPeriodNs);
      h->phase_ns -= fClkPeriodNs;
   }
#endif
   
   //h->fine_ns = h->phase_ns;
   //h->time_sec = fClkPeriodNs*(h->coarse&(~1)) + h->phase_ns;
   h->time_sec = h->coarse_sec + (h->fine_ns + h->offset_ns)*1e-9;

   return true;
};

void DlTdcUnpack::UpdateCalib()
{
   for (size_t i = 0; i < fCalib.size(); i++) {
      fCalib[i].Update();
   }
}

std::string DlTdcUnpack::toJson() const
{
   std::string s;
   s += "{\n";
   char tmp[256];
   sprintf(tmp, "%d", (int)fCalib.size());
   s += "\"dltdc_num_channels\":";
   s += tmp;
   s += ",\n";
   s += "\"dltdc_fine_time\":[\n";

   for (size_t i = 0; i < fCalib.size(); i++) {
      if (i>0)
         s += ",\n";
      s += fCalib[i].toJson().c_str();
   }

   s += "]\n";
   s += "}\n";
   return s;
}

void DlTdcUnpack::SaveCalib(int runno) const
{
   char fname[256];
   sprintf(fname, "dlcfmdb/dltdc_finetime_%06d.json", runno);
   FILE *fp = fopen(fname, "w");
   fprintf(fp, "%s", toJson().c_str());
   fclose(fp);

   //for (size_t i = 0; i < fCalib.size(); i++) {
   //   if (fCalib[i].lepos.fMaxPhase > 0) {
   //      char fname[256];
   //      sprintf(fname, "dltdc_run%d_chan%02d.json", runno, (int)i);
   //      fCalib[i].SaveToFile(fname);
   //   }
   //   if (i>0)
   //      fprintf(fp,",\n");
   //   fprintf(fp, "%s\n", fCalib[i].toJson().c_str());
   //}
}

bool DlTdcUnpack::LoadCalib(const char* fname)
{
   FILE *fp = fopen(fname, "r");

   if (!fp)
      return false;

   std::string json;
   while (1) {
      char buf[1024];
      size_t rd = fread(buf, 1, sizeof(buf)-1, fp);
      //printf("rd %zu\n", rd);
      if (rd == 0)
         break;
      buf[rd] = 0;
      json += buf;
   }
   
   fclose(fp);
      
   //printf("json: %s\n", json.c_str());
   
   MJsonNode* j = MJsonNode::Parse(json.c_str());
   //printf("read: %s\n", j->Stringify().c_str());

   if (j == NULL)
      return false;
   
   int num_channels = j->FindObjectNode("dltdc_num_channels")->GetInt();
   
   printf("DlTdcUnpack::LoadCalib: Loading %d channels from %s\n", num_channels, fname);
   
   fCalib.resize(num_channels);
   
   const std::vector<MJsonNode*>* jj = j->FindObjectNode("dltdc_fine_time")->GetArray();
   
   if (jj->size() != fCalib.size()) {
      delete j;
      return false;
   }
   
   for (size_t i = 0; i < fCalib.size(); i++) {
      bool load_ok = fCalib[i].LoadFromJson((*jj)[i]);
      if (!load_ok) {
         delete j;
         return false;
      }
   }
   
   delete j;
   return true;
}

bool DlTdcUnpack::LoadOffsets(const char* fname)
{
   FILE *fp = fopen(fname, "r");

   if (!fp)
      return false;

   std::string json;
   while (1) {
      char buf[1024];
      size_t rd = fread(buf, 1, sizeof(buf)-1, fp);
      //printf("rd %zu\n", rd);
      if (rd == 0)
         break;
      buf[rd] = 0;
      json += buf;
   }
   
   fclose(fp);
      
   //printf("json: %s\n", json.c_str());
   
   MJsonNode* j = MJsonNode::Parse(json.c_str());
   //printf("read: %s\n", j->Stringify().c_str());

   if (j == NULL)
      return false;
   
   int num_channels = j->FindObjectNode("dltdc_num_channels")->GetInt();
   
   printf("DlTdcUnpack::LoadOffsets: Loading %d channels from %s\n", num_channels, fname);

   fCalib.resize(num_channels);
   
   const std::vector<MJsonNode*>* jj = j->FindObjectNode("offset_ns")->GetArray();
   
   if (jj->size() != fCalib.size()) {
      delete j;
      return false;
   }

   for (size_t i = 0; i < fCalib.size(); i++) {
      const std::vector<MJsonNode*>* jjj = (*jj)[i]->GetArray();
      double le_offset_ns = (*jjj)[0]->GetDouble();
      double te_offset_ns = (*jjj)[1]->GetDouble();
      fCalib[i].lepos.fOffsetNs = fCalib[i].leneg.fOffsetNs = -le_offset_ns;
      fCalib[i].tepos.fOffsetNs = fCalib[i].teneg.fOffsetNs = -te_offset_ns;
   }
   
   delete j;
   return true;
}

bool DlTdcUnpack::LoadCalib(int runno)
{
   for (int r=0; r<100; r++) {

      char fname[256];
      sprintf(fname, "dlcfmdb/dltdc_finetime_%06d.json", runno);

      bool load_ok = LoadCalib(fname);

      if (load_ok)
         return true;
      
      runno--;
   }

   return false;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
