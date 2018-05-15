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

            if (first_ts && fCounter == 0) {
               fAsm[i]->fTsFirstEvent = m->fTs;
               fAsm[i]->fTsLastEvent = 0;
               fAsm[i]->fTsEpoch = 0;
               fAsm[i]->fTimeFirstEvent = m->fTs/ts_freq;
            }

            if (m->fTs < fAsm[i]->fTsLastEvent)
               fAsm[i]->fTsEpoch++;

            m->fTsEpoch = fAsm[i]->fTsEpoch;
            m->fTime = m->fTs/ts_freq - fAsm[i]->fTimeFirstEvent + fAsm[i]->fTsEpoch*2.0*0x80000000/ts_freq;
            m->fTimeIncr = m->fTime - fAsm[i]->fTimeLastEvent;

            fAsm[i]->fTsLastEvent = m->fTs;
            fAsm[i]->fTimeLastEvent = m->fTime;

            if (first_ts) {
               first_ts = false;
               //printf("ts 0x%08x 0x%08x epoch %d time %f\n", m->fTs, m->ts_trig, m->fTsEpoch, m->fTime);
               e->time = m->fTime;
               e->timeIncr = m->fTimeIncr;
            }
         }
      }
   }

   assert(e->modules.size() == e->adcs.size());

   e->counter = ++fCounter;
   e->error = false;
   e->complete = true;

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


