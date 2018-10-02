//
// ALPHA-g TPC
//
// TRB3 TDC functions
//
// Class functions for Tdc.h
//

#include "Tdc.h"

#include <stdio.h>
//#include <string.h>
//#include <assert.h>

static std::string toString(int v)
{
   char buf[256];
   sprintf(buf, "%d", v);
   return buf;
}

void TdcHit::Print() const
{
   printf("TdcHit: fpga %d, chan %2d, re %d, epoch 0x%06x, coarse_time %4d, fine_time %3d", fpga, chan, rising_edge, epoch, coarse_time, fine_time);
}

TdcEvent::TdcEvent() // ctor
{
}

TdcEvent::~TdcEvent() // dtor
{
   for (unsigned i=0; i<hits.size(); i++) {
      if (hits[i]) {
         delete hits[i];
         hits[i] = NULL;
      }
   }
}

void TdcEvent::Print(int level) const
{
   std::string e;
   e += toString(error);
   if (error) {
      e += " (";
      e += error_message;
      e += ")";
   }

   printf("TdcEvent %d, time %.6f, incr %.6f, complete %d, error %s, hits: %d", counter, time, timeIncr, complete, e.c_str(), (int)hits.size());
   if (level > 0) {
      printf("\n");
      for (unsigned i=0; i<hits.size(); i++) {
         printf("  hits[%d]: ", i);
         hits[i]->Print();
         printf("\n");
      }
   }
}

TdcAsm::TdcAsm() // ctor
{
}

TdcAsm::~TdcAsm() // dtor
{
   printf("TdcAsm: %d events\n", fEventCount);
}

void TdcAsm::Print() const
{
   printf("TdcAsm::Print!\n");
}

static uint32_t getUint32(const void* ptr, int offset)
{
  uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
  return (ptr8[0]<<24) | (ptr8[1]<<16) | (ptr8[2]<<8) | ptr8[3];
}

TdcEvent* TdcAsm::UnpackBank(const void* bkptr, int bklen)
{
   TdcEvent* e = new TdcEvent();

   e->counter = ++fEventCount;

   unsigned n32 = bklen/4;

   bool print = false;

   if (print)
      printf("TRBA: length %d bytes, %d words\n", bklen, bklen/4);
   
   int state = 0;
   unsigned cts_count = 0;
   unsigned fpga_count = 0;

   int xfpga = 0;
   uint32_t xepoch = 0;

   uint32_t event_epoch = 0;
   uint32_t event_coarse = 0;
   
   for (unsigned i=0; i<n32; i++) {
      uint32_t v = getUint32(bkptr, i*4);
      //printf("TRBA[%d]: 0x%08x (%d)\n", i, p32[i], p32[i]);
      unsigned threebits = (v>>29)&0x7;
      if (print)
         printf("TRBA[%2d]: 0x%08x: ", i, v);
      
      if (state == 0) {
         if ((i>=6)&&((v&0xFFFF) == 0xC001)) {
            cts_count = (v>>16)&0xFF;
            if (print)
               printf(" CTS data, count %d", cts_count);
            state = 2;
         } else if ((i>=6) && ((v&0xFFFC) == 0x0100)) {
            unsigned fpga = v&0x3;
            xfpga = fpga;
            fpga_count = (v>>16)&0xFF;
            if (print)
               printf(" fpga %d tdc data, count %d", fpga, fpga_count);
            state = 1;
         } else if ((i>=6) && ((v&0xFFFF) == 0x5555)) {
            if (print)
               printf(" end of subevent data");
            state = 3;
         }
      } else if (state == 1) {
         if (print)
            printf(" 3bits: %d", threebits);
         if (threebits == 0) {
            unsigned trigger_type = (v>>(8+16))&0xF; // 4 bits
            unsigned random_code = (v>>16)&0xFF; // 8 bits
            unsigned error_bits = v&0xFFFF; // 16 bits
            if (print)
               printf(" tdc trailer: tt %d, random 0x%02x, errors 0x%04x", trigger_type, random_code, error_bits);
         } else if (threebits == 1) {
            if (print)
               printf(" tdc header");
         } else if (threebits == 3) {
            uint32_t epoch = v & 0x0FFFFFFF;
            xepoch = epoch;
            if (print)
               printf(" epoch counter: %d", epoch);
         } else if (threebits == 4) {
            unsigned chan = (v>>22)&0x7F; // 7 bits
            unsigned fine_time = (v>>12)&0x3FF; // 10 bits
            unsigned rising_edge = (v>>11)&1; // 1 bit
            unsigned coarse_time = (v>>0)&0x7FF; // 11 bits
            if (print)
               printf(" time data: chan %2d, fine %3d, re %1d, coarse %4d", chan, fine_time, rising_edge, coarse_time);
            TdcHit* h = new TdcHit();
            h->fpga = xfpga;
            h->epoch = xepoch;
            h->chan = chan;
            h->fine_time = fine_time;
            h->rising_edge = rising_edge;
            h->coarse_time = coarse_time;
            e->hits.push_back(h);

            if (chan == 0) {
               if (event_epoch == 0) {
                  event_epoch = xepoch;
                  event_coarse = coarse_time;
               } else if (xepoch != event_epoch) {
                  printf("TdcAsm::UnpackBank: Error: channel 0 epoch time mismatch fpga %d has 0x%06x should be 0x%06x\n", xfpga, xepoch, event_epoch);
                  e->error = true;
               } else if ((coarse_time != event_coarse) && (coarse_time+1 != event_coarse) && (coarse_time != event_coarse+1)) {
                  printf("TdcAsm::UnpackBank: Error: channel 0 coarse time mismatch fpga %d has 0x%06x+%d should be 0x%06x+%d\n", xfpga, xepoch, coarse_time, event_epoch, event_coarse);
                  e->error = true;
               }
            }
         }
         if (print)
            printf(" fpga word %d", fpga_count);
         if (fpga_count == 1) {
            state = 0;
         } else {
            fpga_count--;
         }
      } else if (state == 2) {
         if (print)
            printf(" cts word %d", cts_count);
         if (cts_count == 1) {
            state = 0;
         } else {
            cts_count--;
         }
      }
      if (print)
         printf(" state %d\n", state);
   }

   double epoch_freq = 97656.25; // 200MHz/(2<<11);
   double coarse_freq = 200.0e6; // 200MHz

   double event_time = event_epoch/epoch_freq + event_coarse/coarse_freq;

   if (fFirstEventTime == 0) {
      fFirstEventTime = event_time;
      fLastEventTs = 0;
      fLastEventTime = 0;
   }

   e->time = event_time - fFirstEventTime;
   e->timeIncr = e->time - fLastEventTime;

   fLastEventTs = event_epoch;
   fLastEventTime = e->time;

   e->complete = true;

   return e;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
