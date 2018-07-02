//
// Unpacking trigger data
// K.Olchanski
//

#include "TrgAsm.h"

#include <stdio.h> // NULL, printf()
//#include <math.h> // fabs()
//#include <assert.h> // assert()


TrgAsm::TrgAsm() // ctor
{
   Reset();
}

TrgAsm::~TrgAsm() // dtor
{
   // empty
}

void TrgAsm::Reset()
{
   fCounter        = 0; // our packet counter
   fFirstUdpPacket = 0; // seqno of first udp packet
   fFirstTime = 0;
   fPrevTime  = 0;
   fLastTs    = 0;
   fEpoch     = 0;
}

TrigEvent* TrgAsm::UnpackBank(const char* ptr, int size)
{
   TrigEvent* e = new TrigEvent;

   //printf("TrgAsm::UnpackBank: ATAT: ptr %p, size %d\n", ptr, size);
   
   const uint32_t *p32 = (const uint32_t*)ptr;
   const unsigned n32 = size/4;
   for (unsigned i=0; i<n32; i++) {
      //printf("TrgAsm::UnpackBank: ATAT[%d]: 0x%08x (%d)\n", i, p32[i], p32[i]);
      e->udpData.push_back(p32[i]);
   }

   if (e->udpData.size() < 9) {
      e->complete = false;
      e->error = true;
      return e;
   }

   e->complete = true;
   e->error = false;

   if (fCounter == 0) {
      fFirstUdpPacket = e->udpData[0];
   }
   
   e->counter = e->udpData[0] + 1 - fFirstUdpPacket; // udp packet counter counts from 0, we want our counter to count from 1

   double ts_freq = 62.5*1e6; // timestamp is 62.5 MHz

   uint32_t ts = e->udpData[2];

   if (ts < fLastTs) {
      fEpoch++;
   }
   fLastTs = ts;

   if (fCounter == 0) {
      fFirstTime = ts/ts_freq;
      fPrevTime = 0;
   }
   
   e->time = ts/ts_freq - fFirstTime + fEpoch*2.0*0x80000000/ts_freq;
   e->timeIncr = e->time - fPrevTime;
   fPrevTime = e->time;

   fCounter++;

   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


