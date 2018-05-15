//
// Unpacking FEAM data
// K.Olchanski
//

#include "FeamAsm.h"

#include <stdio.h> // NULL, printf()
#include <math.h> // fabs()
#include <assert.h> // assert()

FeamAsm::FeamAsm()
{
   //fAsm.push_back(new FeamModuleAsm);
}

FeamAsm::~FeamAsm()
{
   printf("FeamAsm::~FeamAsm: Total events: %d, complete: %d, incomplete: %d, with errors: %d, max timestamp difference: %.0f ns\n", fCounter, fCountComplete, fCountIncomplete, fCountError, fMaxDt*1e9);

   //Print();

   for (unsigned i=0; i<fAsm.size(); i++) {
      if (fAsm[i]) {
         delete fAsm[i];
         fAsm[i] = NULL;
      }
   }
}

void FeamAsm::AddPacket(int imodule, int icolumn, int iring, int format, const FeamPacket* p, const char* ptr, int size)
{
   assert(imodule >= 0);
   while (imodule >= fAsm.size()) {
      fAsm.push_back(new FeamModuleAsm());
   }
   assert(imodule < fAsm.size());

   fAsm[imodule]->AddPacket(p, imodule, imodule, icolumn, iring, format, ptr, size);

   delete p;
}

void FeamAsm::BuildEvent(FeamEvent* e)
{
   double ts_freq = 125.0e6; // 125MHz timestamp clock

   e->error = false;
   e->complete = true;

   bool first_ts = true;
   for (unsigned i=0; i<fAsm.size(); i++) {
      if (fAsm[i]) {
         fAsm[i]->Finalize();

         if (fAsm[i]->fBuffer.size() > 0) {
            assert(fAsm[i]->fBuffer.size() == 1);
            FeamModuleData* m = fAsm[i]->fBuffer.front();
            fAsm[i]->fBuffer.pop_front();

            FeamAdcData *a = new FeamAdcData;

            Unpack(a, m);

            e->modules.push_back(m);
            e->adcs.push_back(a);

            m->fTs = m->ts_trig;

            if (fCounter == 0) {
               fAsm[i]->fTsFirstEvent = m->fTs;
               fAsm[i]->fTsLastEvent = 0;
               fAsm[i]->fTsEpoch = 0;
               fAsm[i]->fTimeFirstEvent = m->fTs/ts_freq;
               fAsm[i]->fCntFirstEvent = m->cnt;
            }

            if (m->fTs < fAsm[i]->fTsLastEvent)
               fAsm[i]->fTsEpoch++;

            m->fTsEpoch = fAsm[i]->fTsEpoch;
            m->fTime = m->fTs/ts_freq - fAsm[i]->fTimeFirstEvent + fAsm[i]->fTsEpoch*2.0*0x80000000/ts_freq;
            m->fTimeIncr = m->fTime - fAsm[i]->fTimeLastEvent;

            m->fCntSeq = m->cnt - fAsm[i]->fCntFirstEvent;

            fAsm[i]->fTsLastEvent = m->fTs;
            fAsm[i]->fTimeLastEvent = m->fTime;

            if (first_ts) {
               first_ts = false;
               //printf("ts 0x%08x 0x%08x epoch %d time %f\n", m->fTs, m->ts_trig, m->fTsEpoch, m->fTime);
               e->time = m->fTime;
               e->timeIncr = m->fTimeIncr;
            }

            if (m->error)
               e->error = true;

            if (!m->complete)
               e->complete = false;
         }
      }
   }

   e->counter = ++fCounter;

   // check for consistency

   assert(e->modules.size() == e->adcs.size());

   for (unsigned i=0; i<e->modules.size(); i++) {
      double dt = e->modules[i]->fTime - e->time;
      double absdt = fabs(dt);
      //printf("YYY event time %f: module %d time %f diff %f\n", e->time, i, e->modules[i]->fTime, dt);
      if (absdt > fMaxDt)
         fMaxDt = absdt;
      if (absdt > fConfMaxDt) {
         printf("FeamAsm::BuildEvent: event %d timestamp mismatch module %d: time %f diff %f should be %f\n", e->counter, i, e->modules[i]->fTime, dt, e->time);
         e->error = true;
      }

      uint32_t cnt = e->modules[i]->cnt;
      int cntseq = e->modules[i]->fCntSeq;
      //printf("YYY event %d module %d cnt %d, seq %d\n", e->counter, i, cnt, cntseq);

      if (cntseq+1 != e->counter) {
         printf("FeamAsm::BuildEvent: event %d event counter mismatch module %d: cnt 0x%08x, seqcnt %d should be %d\n", e->counter, i, cnt, cntseq, e->counter);
         e->error = true;
      }
   }

   // increment counters

   if (e->error)
      fCountError++;
   if (e->complete)
      fCountComplete++;
   else
      fCountIncomplete++;
}

void FeamAsm::Print() const
{
   printf("FeamAsm::Print: FEAM event builder status:\n");

   printf("  event counter: %d\n", fCounter);

   printf("  complete events: %d\n", fCountComplete);
   printf("  incomplete events: %d\n", fCountIncomplete);
   printf("  events with error: %d\n", fCountError);
   
   printf("  timestamp check: max difference: %.0f ns, max permitted: %.0f ns\n", fMaxDt*1e9, fConfMaxDt*1e9);

   printf("  Assembler: %d entries\n", (int)fAsm.size());
   for (unsigned i=0; i<fAsm.size(); i++) {
      printf("    position %2d: ", i);
      fAsm[i]->Print();
      printf("\n");
   }

   printf("FeamAsm::Print: done.\n");
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


