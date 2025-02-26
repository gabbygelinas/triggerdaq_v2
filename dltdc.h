#ifndef DLTDC_H
#define DLTDC_H

#include <stdint.h> // uint32_t
#include <stdio.h> // printf()
#include <vector> // std::vector
#include <string> // std::string

struct MJsonNode;

struct DlTdcHit
{
   uint32_t data_lo = 0;
   uint32_t data_hi = 0;
   int ch = 0;
   bool le = false;
   bool te = false;
   bool ts0 = false;
   uint64_t coarse = 0;
   int phase = 0;

   double coarse_epoch = 0;
   double coarse_sec = 0;
   double fine_ns  = 0;
   double offset_ns = 0;
   double time_sec = 0;

   void Clear()
   {
      data_lo = 0;
      data_hi = 0;
      ch = 0;
      le = false;
      te = false;
      ts0 = false;
      coarse = 0;
      phase = 0;
      coarse_epoch = 0;
      coarse_sec = 0;
      fine_ns  = 0;
      offset_ns = 0;
      time_sec = 0;
   }

   void Unpack(uint32_t lo, uint32_t hi);
   void Print() const;

   static const uint64_t kEpoch = 0x40000000;
};

class DlTdcFineCalib1
{
public:
   double fTotalNs = 10.0;
   size_t fHits = 0;
   std::vector<double> fHistogram;
   std::vector<double> fBinWidthNs;
   std::vector<double> fBinTimeNs;
   int    fMaxPhase = 0;
   double fBinMinNs = 0;
   double fBinMaxNs = 0;
   double fOffsetNs = 0;

public:
   DlTdcFineCalib1(); // ctor
   ~DlTdcFineCalib1(); // dtor

public:
   void Resize(int nbins);
   void Reset();
   void AddHit(int phase, double fine_ns);
   void Update();
   void Print() const;
   double GetTime(int phase);
   std::string toJson() const;
};

class DlTdcFineCalib
{
public:
   DlTdcFineCalib1 lepos;
   DlTdcFineCalib1 leneg;
   DlTdcFineCalib1 tepos;
   DlTdcFineCalib1 teneg;

public:
   DlTdcFineCalib(); // ctor
   ~DlTdcFineCalib(); // dtor

public:
   void Reset();
   void AddHit(const DlTdcHit& h);
   void Update();
   void Print() const;
   void SaveToFile(const char* filename) const;
   bool LoadFromFile(const char* filename);
   bool LoadFromJson(const MJsonNode* j);
   std::string toJson() const;
};

class DlTdcUnpack
{
public: // configuration
   double fClkPeriodNs = 10.0;

public:
   DlTdcUnpack(int nchan); // ctor
   ~DlTdcUnpack(); // dtor

private:
   DlTdcUnpack(); // ctor

public:
   void Reset();
   //bool Unpack(DlTdcHit* h, uint32_t lo, uint32_t hi);
   bool ComputeTimes(DlTdcHit* h);
   void UpdateCalib();
   bool LoadCalib(const char* fname);
   bool LoadOffsets(const char* fname);
   bool LoadCalib(int runno);
   void SaveCalib(int runno) const;
   std::string toJson() const;

public: // internal state
   double fFirstTimeSec = 0;
   std::vector<uint64_t> fLastCoarse;
   std::vector<double>   fEpoch;
   std::vector<double>   fEpochHits;
   std::vector<DlTdcFineCalib> fCalib;

public:
   static void PrintBits32(uint32_t v);
   static uint32_t FixHoles(uint32_t sr);
   static int FindEdge(uint32_t sr);
   static int FindEdge10(uint32_t sr);
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */

