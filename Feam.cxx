//
// Unpacking FEAM data
// K.Olchanski
//

#include "Feam.h"

#include <stdio.h> // NULL
#include <stdlib.h> // realloc()
#include <string.h> // memcpy()
#include <assert.h> // assert()
#include <utility>  // std::pair
#include <iostream>

#if 0
static uint8_t getUint8(const void* ptr, int offset)
{
   return *(uint8_t*)(((char*)ptr)+offset);
}

static uint16_t getUint16be(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return ((ptr8[0]<<8) | ptr8[1]);
}

static uint32_t getUint32be(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return (ptr8[0]<<24) | (ptr8[1]<<16) | (ptr8[2]<<8) | ptr8[3];
}
#endif

static uint16_t getUint16le(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return ((ptr8[1]<<8) | ptr8[0]);
}

static uint32_t getUint32le(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return (ptr8[3]<<24) | (ptr8[2]<<16) | (ptr8[1]<<8) | ptr8[0];
}

padMap::padMap(){
   const int first = 0;  // change if index is numbered from 1
   for(int i = 0; i < MAX_FEAM_READOUT; i++){
      if(i < 4) channel[i] = -99;
      else {
         if(i == 16) channel[i] = -1;
         else if(i == 29) channel[i] = -2;
         else if(i == 54) channel[i] = -3;
         else if(i == 67) channel[i] = -4;
         else {
            int ch = i - int(i > 16) - int(i > 29) - int(i > 54) - int(i > 67) - 3;
            channel[i] = ch;
            readout[ch] = i;
         }
      }
   }
   for(int isca = 0; isca < MAX_FEAM_SCA; isca++){
      padcol[isca][0] = -99;
      padrow[isca][0] = -99;
   }
   for(int isca = 0; isca < MAX_FEAM_SCA; isca++){
      int offset = (isca%2)*MAX_FEAM_PAD_ROWS/2;
      for(int ch = 1; ch <= MAX_FEAM_CHAN; ch++){
         int col, row;
         if(ch <= 36){
            if(ch > 18){
               col = 0;
               row = 36-ch+offset;
            } else {
               col = 1;
               row = ch-1+offset;
            }
         } else {
            if(ch < 55){
               col = 0;
               row = 72-ch+offset;
            } else {
               col = 1;
               row = ch-37+offset;
            }
         }
         if(isca > 1){
            col = 3-col;
            row = 71-row;
         }
         padcol[isca][ch] = col;
         padrow[isca][ch] = row;
         sca[col][row] = isca;
         sca_chan[col][row] = ch;
      }
   }
};

//static int x1count = 0;

FeamPacket::FeamPacket()
{
   //printf("FeamPacket: ctor!\n");
   //printf("FeamPacket: count %d!\n", x1count++);
   error = true;
}

FeamPacket::~FeamPacket()
{
   //printf("FeamPacket: dtor!\n");
   //x1count--;
}

void FeamPacket::Unpack(const char* data, int size)
{
   error = true;

   off = 0;
   cnt = getUint32le(data, off); off += 4;
   n   = getUint16le(data, off); off += 2;
   x511 = getUint16le(data, off); off += 2;
   buf_len = getUint16le(data, off); off += 2;
   if (n == 0) {
      ts_start = getUint32le(data, off); off += 8;
      ts_trig  = getUint32le(data, off); off += 8;
   } else {
      ts_start = 0;
      ts_trig  = 0;
   }

   error = false;
}

void FeamPacket::Print() const
{
   printf("decoded %2d bytes, ", off);
   printf("cnt %6d, n %3d, x511 %3d, buf_len %4d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          n,
          x511,
          buf_len,
          ts_start,
          ts_trig);
   printf("error %d", error);
}

//static int x2count = 0;

void FeamModuleData::Finalize()
{
   if (error) {
      return;
   }

   if (next_n != 256) {
      error = true;
      return;
   }

   if (fSize != 310688) {
      error = true;
      return;
   }
}

/*
   ZZZ Run 421
   ZZZ Processing FEAM event: module  4, cnt     84, ts_start 0xc8f30f8a, ts_trig 0xc8f311fa, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  5, cnt     97, ts_start 0xc8f6d5d4, ts_trig 0xc8f6d844, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  6, cnt     93, ts_start 0xc8f91bb4, ts_trig 0xc8f91e24, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  7, cnt     83, ts_start 0xc8f5c6dc, ts_trig 0xc8f5c94c, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  0, cnt     80, ts_start 0xc8f46ecc, ts_trig 0xc8f4713c, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  1, cnt   1465, ts_start 0xc8f2ce84, ts_trig 0xc8f2d0f4, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  2, cnt     63, ts_start 0xc8f2892c, ts_trig 0xc8f28b9c, next_n 256, size 310688, error 0
   ZZZ Processing FEAM event: module  3, cnt     62, ts_start 0xc8ef4a56, ts_trig 0xc8ef4cc6, next_n 256, size 310688, error 0
*/

void FeamModuleData::Print(int level) const
{
   printf("bank %s, module %2d, ", fBank.c_str(), fPosition);
   printf("cnt %6d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          ts_start,
          ts_trig);
   printf("next_n %d, ", next_n);
   printf("size %d, ", fSize);
   printf("error %d", error);

   if (level > 0) {
      printf("\n");
      printf("ADC data:\n");
      const uint16_t* aptr16 = (uint16_t*)fPtr;
      int asize16 = fSize/2;
      for (int i=0; i<asize16; i++) {
         if (i%4 == 0)
            printf("%d: ", i);
         printf(" 0x%04x", aptr16[i]);
         if (i%4 == 3)
            printf("\n");
      }
      printf("ADC data done\n");
   }
}

void FeamModuleData::AddData(const FeamPacket*p, int position, const char* ptr, int size)
{
   assert(position == fPosition);

   //printf("add %d size %d\n", p->n, size);
   if (p->n != next_n) {
      printf("FeamModuleData::AddData: position %2d, cnt %6d, wrong packet sequence: expected %d, got %d!\n", position, cnt, next_n, p->n);
      next_n = p->n + 1;
      error = true;
      return;
   }

   assert(size >= 0);
   assert(size < 12000); // UDP packet size is 1500 bytes, jumbo frame up to 9000 bytes

   int new_size = fSize + size;
   char* new_ptr = (char*)realloc(fPtr, new_size);

   if (!new_ptr) {
      printf("FeamModuleData::AddData: cannot reallocate ADC buffer from %d to %d bytes!\n", fSize, new_size);
      error = true;
      return;
   }

   memcpy(new_ptr + fSize, ptr, size);

   fPtr = new_ptr;
   fSize = new_size;

   next_n = p->n + 1;
}

FeamModuleData::FeamModuleData(const FeamPacket* p, const char* bank, int position)
{
   //printf("FeamModuleData: ctor! %d\n", x2count++);
   assert(p->n == 0);

   fBank = bank;
   fPosition = position;

   cnt = p->cnt;
   ts_start = p->ts_start;
   ts_trig = p->ts_trig;

   fSize = 0;
   fPtr = NULL;

   next_n = 0;

   error = false;

   fTs = 0;
   fTsEpoch = 0;
   fTime = 0;
   fTimeIncr = 0;
}

FeamModuleData::~FeamModuleData() // dtor
{
   //printf("FeamModuleData: dtor!\n"); x2count--;
   if (fPtr)
      free(fPtr);
   fPtr = NULL;
   fSize = 0;
}

FeamEvent::FeamEvent() // ctor
{
   complete = false;
   error = false;
   counter = 0;
   time = 0;
   timeIncr = 0;
}

FeamEvent::~FeamEvent() // dtor
{
   for (unsigned i=0; i<modules.size(); i++)
      if (modules[i]) {
         delete modules[i];
         modules[i] = NULL;
      }

   for (unsigned i=0; i<adcs.size(); i++)
      if (adcs[i]) {
         delete adcs[i];
         adcs[i] = NULL;
      }

   for (unsigned i=0; i<hits.size(); i++)
      if (hits[i]) {
         delete hits[i];
         hits[i] = NULL;
      }
}

void FeamEvent::Print(int level) const
{
   printf("FeamEvent %d, time %f, incr %f, complete %d, error %d, modules: ", counter, time, timeIncr, complete, error);
   for (unsigned i=0; i<modules.size(); i++) {
      if (modules[i] == NULL) {
            printf(" null");
      } else {
         printf(" %d", modules[i]->cnt);
      }
   }
}

void Unpack(FeamAdcData* a, FeamModuleData* m)
{
   a->nsca  = 4; // 0,1,2,3
   a->nchan = 79; // 1..79 readout channels per table 2 in AFTER SCA manual
   a->nbins = 511; // 0..511

   memset(a->adc, 0, sizeof(a->adc));

   const unsigned char* ptr = (const unsigned char*)m->fPtr;
   int count = 0;

   for (int ibin = 0; ibin < 511; ibin++) {
      for (int ichan = 0; ichan <= 3; ichan++) {
         for (int isca = 0; isca < 4; isca++) {
            a->adc[isca][ichan][ibin] = 0;
         }
      }
   }

   for (int ibin = 0; ibin < 511; ibin++) {
      for (int ichan = 4; ichan <= 79; ichan++) {
         for (int isca = 0; isca < 4; isca++) {
            unsigned v = ptr[0] | ((ptr[1])<<8);
            // manual sign extension
            if (v & 0x8000)
               v |= 0xffff0000;
            //if (isca == 0) {
            //   adc[ichan][ibin] = v;
            //}
            a->adc[isca][ichan][ibin] = v;
            ptr += 2;
            count += 2;
         }
      }
   }

   //printf("count %d\n", count);
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
