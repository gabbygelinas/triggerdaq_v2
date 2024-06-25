//
// DL TDC event. K.Olchanski
//

#ifndef DlTdcEvent_H
#define DlTdcEvent_H

#include "dltdc.h"

class DlTdcHit2
{
public:
   DlTdcHit fLe;
   DlTdcHit fTe;
   bool     fUp   = false;
   bool     fDown = false;
   int      fCount = 0;

public:
   double fTimeSec  = 0;
   double fWidthNs = 0;

public:
   void Clear();
   void AddHit(const DlTdcHit& h);
   void Print() const;
};

class DlTdcEvent
{
public: // matching with ADC data
   double time_sec = 0;
   double dt = 0;

public: // time of first and last TDC hits
   double first_time_sec = 0;
   double last_time_sec = 0;

   double min_time_sec = 0;
   double max_time_sec = 0;

public: // TDC hits
   std::vector<DlTdcHit2> fHits;

public:
   void Init(int num_chan);
   void Clear();
   void AddHit8(const DlTdcHit& h);
   bool HaveCh(int ch) const;
   const DlTdcHit2& GetCh(int ch) const;
};

double sec_to_ns(double t);
double subtract_ns(const DlTdcHit& h1, const DlTdcHit& h2);

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


