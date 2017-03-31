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

void FeamModuleData::Print() const
{
   printf("module %2d, ", module);
   printf("cnt %6d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          ts_start,
          ts_trig);
   printf("next_n %d, ", next_n);
   printf("size %d, ", fSize);
   printf("error %d", error);
}

void FeamModuleData::AddData(const FeamPacket*p, int xmodule, const char* ptr, int size)
{
   assert(xmodule == module);

   //printf("add %d size %d\n", p->n, size);
   if (p->n != next_n) {
      printf("module %2d, cnt %6d, wrong packet sequence: expected %d, got %d!\n", module, cnt, next_n, p->n);
      next_n = p->n + 1;
      error = true;
      return;
   }

   assert(size >= 0);
   assert(size < 12000); // UDP packet size is 1500 bytes, jumbo frame up to 9000 bytes

   int new_size = fSize + size;
   char* new_ptr = (char*)realloc(fPtr, new_size);

   if (!new_ptr) {
      printf("cannot reallocate ADC buffer from %d to %d bytes!\n", fSize, new_size);
      error = true;
      return;
   }

   memcpy(new_ptr + fSize, ptr, size);

   fPtr = new_ptr;
   fSize = new_size;

   next_n = p->n + 1;
}

FeamModuleData::FeamModuleData(const FeamPacket* p, int xmodule)
{
   //printf("FeamModuleData: ctor! %d\n", x2count++);
   assert(p->n == 0);

   module = xmodule;

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
}

void FeamEvent::Print() const
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
   a->nsca  = 4;
   a->nchan = 76;
   a->nbins = 511;

   const unsigned char* ptr = (const unsigned char*)m->fPtr;
   int count = 0;

   for (int ibin = 0; ibin < 511; ibin++) {
      for (int ichan = 0; ichan < 76; ichan++) {
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

std::pair<int,int> getPad(short fAFTER, int index){
    // int channel = index; // FIXME: change to next line once basic testing is done
    const int first = 0;  // change if index is numbered from 1
    index -= first;

    // Remove "special" channels to get 72 consecutive indices
    int channel = index;
    if(index == 12) channel = -1;
    else if(index == 25) channel = -2;
    else if(index == 50) channel = -3;
    else if(index == 63) channel = -4;
    else if(index > 76){
       std::cerr << "Index " << index << " too large!" << std::endl;
       return std::pair<int,int>(-999,-999);
    } else
        channel -= int(index > 11) + int(index > 24) + int(index > 49) + int(index > 62);

    if(channel < 0){
       return std::pair<int,int>(-(fAFTER+1), channel);
    }
    int col = -1;
    int pad = -1;
    int offset = 0;
    if(fAFTER < 2){
        if(fAFTER == 0) offset = 36;
        if(channel < 36){
            if(channel >= 18){
                col = 0;
                pad = channel+offset;
            } else {
                col = 1;
                pad = 35-channel+offset;
            }
        } else {
            if(channel < 54){
                col = 0;
                pad = channel-36+offset;
            } else {
                col = 1;
                pad = 71-channel+offset;
            }
        }
    } else {
        if(fAFTER == 3) offset = 36;
        if(channel < 36){
            if(channel < 18){
                col = 2;
                pad = channel+offset;
            } else {
                col = 3;
                pad = 35-channel+offset;
            }
        } else {
            if(channel >= 54){
                col = 2;
                pad = channel-36+offset;
            } else {
                col = 3;
                pad = 71-channel+offset;
            }
        }
    }


    return std::pair<int,int>(col,pad);
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
