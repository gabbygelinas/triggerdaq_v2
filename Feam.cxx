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
   printf("fmt %d, ", fDataFormat);
   printf("cnt %6d, ts_start 0x%08x, ts_trig 0x%08x, ",
          cnt,
          ts_start,
          ts_trig);
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

int FeamModuleData::fgMaxAlloc = 0;

void FeamModuleData::AddData(const FeamPacket*p, const char* ptr, int size)
{
   //printf("add %d size %d\n", p->n, size);

   assert(size >= 0);
   assert(size < 12000); // UDP packet size is 1500 bytes, jumbo frame up to 9000 bytes

   int new_size = fSize + size;

   if (new_size > fAlloc) {
      if (new_size < fgMaxAlloc)
         new_size = fgMaxAlloc;
      //printf("realloc %d -> %d, max %d\n", fAlloc, new_size, fgMaxAlloc);
      char* new_ptr = (char*)realloc(fPtr, new_size);

      if (!new_ptr) {
         printf("FeamModuleData::AddData: cannot reallocate ADC buffer from %d to %d bytes!\n", fSize, new_size);
         error = true;
         return;
      }

      fAlloc = new_size;
      fPtr = new_ptr;

      if (fAlloc > fgMaxAlloc)
         fgMaxAlloc = fAlloc;
   }

   memcpy(fPtr + fSize, ptr, size);
   fSize += size;
}

FeamModuleData::FeamModuleData(const FeamPacket* p, const char* bank, int position, int format)
{
   //printf("FeamModuleData: ctor! %d\n", x2count++);
   assert(p->n == 0);

   fBank = bank;
   fPosition = position;
   fDataFormat = format;

   cnt = p->cnt;
   ts_start = p->ts_start;
   ts_trig = p->ts_trig;
}

FeamModuleData::~FeamModuleData() // dtor
{
   //printf("FeamModuleData: dtor!\n"); x2count--;
   if (fPacket) {
      delete fPacket;
      fPacket = NULL;
   }
   if (fPtr)
      free(fPtr);
   fPtr = NULL;
   fSize = 0;
   fAlloc = 0;
}

#define ST_INIT  0
#define ST_DATA  1
#define ST_WAIT  2
#define ST_DONE  3

FeamAsm::~FeamAsm()
{
   fState = -1;
   fCnt = 0;
   fNextN = -1;
   if (fCurrent) {
      delete fCurrent;
      fCurrent = NULL;
   }
   for (unsigned i=0; i<fBuffer.size(); i++) {
      if (fBuffer[i]) {
         delete fBuffer[i];
         fBuffer[i] = NULL;
      }
   }
}

void FeamAsm::Print() const
{
   int countComplete = 0;
   int countError = 0;
   for (unsigned i=0; i<fBuffer.size(); i++) {
      if (fBuffer[i]->complete)
         countComplete++;
      if (fBuffer[i]->error)
         countError++;
   }
   printf("pos %d, bank %s, state %d, cnt %d, nextn %d, ig %d, fi %d, do %d (sy %d, tr %d, sk %d, wcnt %d), cur %p, buf %d (com %d, err %d)", fPosition, fBank.c_str(), fState, fCnt, fNextN, fCountIgnoredBeforeFirst, fCountFirst, fCountDone, fCountLostSync, fCountTruncated, fCountSkip, fCountWrongCnt, fCurrent, (int)fBuffer.size(), countComplete, countError);
}

void FeamAsm::StFirstPacket(const FeamPacket* p, const char* bank, int position, int format, const char* ptr, int size)
{
   fState = ST_DATA;
   fCnt = p->cnt;
   fNextN = 1;
   fCountFirst++;

   fPosition = position;
   fDataFormat = format;
   fBank = bank;

   assert(fCurrent == NULL);

   fCurrent = new FeamModuleData(p, bank, position, format);
   fCurrent->fPacket = new FeamPacket(*p); // copy constructor

   fCurrent->AddData(p, ptr, size);
}

void FeamAsm::StLastPacket()
{
   fState = ST_DONE;
   fCountDone++;

   assert(fCurrent != NULL);

   fCurrent->complete = true;
   fBuffer.push_back(fCurrent);
   fCurrent = NULL;
}

void FeamAsm::AddData(const FeamPacket* p, const char* ptr, int size)
{
   assert(fCurrent != NULL);
   fCurrent->AddData(p, ptr, size);
}

void FeamAsm::FlushIncomplete()
{
   assert(fCurrent != NULL);

   fCurrent->complete = false;
   fCurrent->error = true;
   fBuffer.push_back(fCurrent);
   fCurrent = NULL;
}

void FeamAsm::AddPacket(const FeamPacket* p, const char* bank, int position, int format, const char* ptr, int size)
{
   bool trace = false;
   bool traceNormal = false;
   
   switch (fState) {
   default: {
      assert(!"invalid state!");
      break;
   }
   case ST_INIT: { // initial state, before received first first packet
      if (p->n==0) { // first packet
         if (trace)
            printf("ST_INIT: bank %s, packet cnt %d, n %d ---> first packet\n", bank, p->cnt, p->n);
         StFirstPacket(p, bank, position, format, ptr, size);
      } else {
         if (traceNormal)
            printf("ST_INIT: bank %s, packet cnt %d, n %d\n", bank, p->cnt, p->n);
         fCountIgnoredBeforeFirst++;
      }
      break;
   }
   case ST_DATA: { // receiving data
      //printf("ST_FIRST: bank %s, packet cnt %d, n %d\n", bank, p->cnt, p->n);
      if (p->n == 0) { // unexpected first packet
         if (trace)
            printf("ST_DATA: bank %s, packet cnt %d, n %d ---> unexpected first packet\n", bank, p->cnt, p->n);
         fCountTruncated++;
         FlushIncomplete();
         StFirstPacket(p, bank, position, format, ptr, size);
      } else if (p->cnt != fCnt) { // packet from wrong event
         if (trace)
            printf("ST_DATA: bank %s, packet cnt %d, n %d ---> wrong cnt expected cnt %d\n", bank, p->cnt, p->n, fCnt);
         fState = ST_WAIT;
         FlushIncomplete();
         fCountWrongCnt++;
      } else if (p->n != fNextN) { // out of sequence packet
         if (trace)
            printf("ST_DATA: bank %s, packet cnt %d, n %d ---> out of sequence expected n %d\n", bank, p->cnt, p->n, fNextN);
         fState = ST_WAIT;
         FlushIncomplete();
         fCountLostSync++;
      } else if (p->n == 255) { // last packet
         if (trace)
            printf("ST_DATA: bank %s, packet cnt %d, n %d ---> last packet\n", bank, p->cnt, p->n);
         AddData(p, ptr, size);
         fNextN++;
         StLastPacket();
      } else {
         if (traceNormal)
            printf("ST_DATA: bank %s, packet cnt %d, n %d\n", bank, p->cnt, p->n);
         AddData(p, ptr, size);
         fNextN++;
      }
      break;
   }
   case ST_WAIT: { // skipping bad data
      if (p->n == 0) { // first packet
         if (trace)
            printf("ST_WAIT: bank %s, packet cnt %d, n %d ---> first packet\n", bank, p->cnt, p->n);
         StFirstPacket(p, bank, position, format, ptr, size);
      } else {
         if (traceNormal)
            printf("ST_WAIT: bank %s, packet cnt %d, n %d\n", bank, p->cnt, p->n);
         fCountSkip++;
      }
      break;
   }
   case ST_DONE: { // received last packet
      if (p->n == 0) { // first packet
         if (trace)
            printf("ST_DONE: bank %s, packet cnt %d, n %d ---> first packet\n", bank, p->cnt, p->n);
         StFirstPacket(p, bank, position, format, ptr, size);
      } else {
         if (trace)
            printf("ST_DONE: bank %s, packet cnt %d, n %d ---> lost first packet\n", bank, p->cnt, p->n);
         fState = ST_WAIT;
         fCountLostFirst++;
      }
      break;
   }
   } // switch (fState)
}

void FeamAsm::Finalize()
{
   switch (fState) {
   default: {
      assert(!"invalid state!");
      break;
   }
   case ST_INIT: { // initial state, before received first first packet
      break;
   }
   case ST_DATA: { // receiving data
      FlushIncomplete();
      fState = ST_DONE;
      break;
   }
   case ST_WAIT: { // skipping bad data
      break;
   }
   case ST_DONE: { // received last packet
      break;
   }
   } // switch (fState)
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
   printf("PwbEvent %d, time %f, incr %f, complete %d, error %d, modules: ", counter, time, timeIncr, complete, error);
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
   int f = m->fDataFormat;
   assert(f==1 || f==2);

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
            if (f==1)
               a->adc[isca][ichan][ibin] = ((int)v)/16;
            else if (f==2)
               a->adc[isca][ichan][ibin] = v;
            else
               a->adc[isca][ichan][ibin] = 0xdead;
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
