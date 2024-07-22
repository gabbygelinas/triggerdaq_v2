//
// DlTdcEvent.cxx
//
// K.Olchanski
//

#include <stdio.h>
#include <assert.h>
#include <unistd.h>

#include "DlTdcEvent.h"

double sec_to_ns(double t)
{
   return t*1e9;
}

double subtract_ns(const DlTdcHit& h1, const DlTdcHit& h2)
{
   return (h1.coarse_sec- h2.coarse_sec)*1e9 + ((h1.fine_ns + h1.offset_ns) - (h2.fine_ns + h2.offset_ns));
}

void DlTdcHit2::Clear()
{
   fLe.Clear();
   fTe.Clear();
   fUp   = false;
   fDown = false;
   fCount = 0;
   fTimeSec  = 0;
   fWidthNs = 0;
}

void DlTdcHit2::AddHit(const DlTdcHit& h)
{
   bool complain = true; // false;
   
   if (h.le) {
      if (!fUp) {
         fUp = true;
         fDown = false;
         if (fCount == 0) {
            fLe = h;
            fTimeSec = h.time_sec;
         } else {
            double te_to_le_ns = sec_to_ns(h.time_sec - fTe.time_sec);
            if (te_to_le_ns < 5.0) {
               //if (complain)
               //   printf("TTT: ch %d: LE-TE-LE, MERGE TE to LE %.3f ns\n", h.ch, te_to_le_ns);
               fCount = 0;
            } else {
               if (complain)
                  printf("TTT: ch %d: MULTIPLE HIT, time %.9f -> %.9f sec, dt %.3f ns, TE to LE %.3f ns, count %d\n", h.ch, fTimeSec, h.time_sec, sec_to_ns(h.time_sec - fTimeSec), sec_to_ns(h.time_sec - fTe.time_sec), fCount);
               //fLe.Print(); printf("\n");
               //fTe.Print(); printf("\n");
               //h.Print(); printf("\n");
            }
         }
      } else {
         double le_le_dt_ns = sec_to_ns(h.time_sec - fLe.time_sec);
         if (le_le_dt_ns < 80.0) {
            //if (complain)
            //   printf("TTT: ch %d: LE-LE %.3f ns\n", h.ch, le_le_dt_ns);
         } else {
            if (complain)
               printf("TTT: ch %d: LE-LE, LE FROM WRONG EVENT, LE-LE %.3f ns, count %d\n", h.ch, le_le_dt_ns, fCount);
         }
      }
   } else if (h.te) {
      if (fUp) {
         fUp = false;
         fDown = true;
         if (fCount == 0) {
            fTe = h;
            fWidthNs = subtract_ns(fTe, fLe);
         } else {
            if (complain)
               printf("TTT: ch %d: TE-TE %.3f ns, multiple hit, count %d\n", h.ch, sec_to_ns(h.time_sec - fTe.time_sec), fCount);
         }
         fCount++;
      } else {
         if (fCount == 0) {
            if (complain)
               printf("TTT: ch %d: TE without LE\n", h.ch);
         } else {
            double le_te_dt_ns = sec_to_ns(fTe.time_sec - fLe.time_sec);
            double te_te_dt_ns = sec_to_ns(h.time_sec - fTe.time_sec);
            if (te_te_dt_ns < 80.0) {
               if (complain)
                  printf("TTT: ch %d: LE-TE-TE, MISSING LE, LE-TE %.3f, TE-TE %.3f ns, count %d\n", h.ch, le_te_dt_ns, te_te_dt_ns, fCount);
            } else {
               if (complain)
                  printf("TTT: ch %d: LE-TE-TE, TE FROM WRONG EVENT, LE-TE %.3f, TE-TE %.3f ns, count %d\n", h.ch, le_te_dt_ns, te_te_dt_ns, fCount);
            }
            //fLe.Print(); printf("\n");
            //fTe.Print(); printf("\n");
            //h.Print(); printf("\n");
         }
      }
   }
};

void DlTdcHit2::Print() const
{
   printf("DlTdcHit2: u/d %d%d, time_sec %.9f, w_ns %.3f\n", fUp, fDown, fTimeSec, fWidthNs);
   fLe.Print();
   printf("\n");
   fTe.Print();
   printf("\n");
};

void DlTdcEvent::Init(int num_chan)
{
   fHits.resize(num_chan);
}

void DlTdcEvent::Clear()
{
   first_time_sec = 0;
   last_time_sec = 0;
   
   min_time_sec = 0;
   max_time_sec = 0;
   
   for (size_t i=0; i<fHits.size(); i++) {
      fHits[i].Clear();
   }
}

void DlTdcEvent::AddHit8(const DlTdcHit& h)
{
   if (first_time_sec == 0) {
      first_time_sec = h.time_sec;
      min_time_sec = h.time_sec;
      max_time_sec = h.time_sec;
   }
   
   last_time_sec = h.time_sec;
   
   if (h.time_sec < min_time_sec)
      min_time_sec = h.time_sec;
   
   if (h.time_sec > max_time_sec)
      max_time_sec = h.time_sec;
   
   assert(h.ch >= 0);
   assert(h.ch < (int)fHits.size());
   
   fHits[h.ch].AddHit(h);
}

bool DlTdcEvent::HaveCh(int ch) const
{
   return fHits[ch].fDown;
}

const DlTdcHit2& DlTdcEvent::GetCh(int ch) const
{
   return fHits[ch];
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
