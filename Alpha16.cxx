//
// ALPHA-g TPC
//
// GRIF-16/ALPHA-16 ADC functions
//
// Class functions for Alpha16.h
//

#include "Alpha16.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

#define MEMZERO(array) memset((array),0,sizeof(array))

static uint8_t getUint8(const void* ptr, int offset)
{
   return *(uint8_t*)(((char*)ptr)+offset);
}

static uint16_t getUint16(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return ((ptr8[0]<<8) | ptr8[1]);
}

static uint32_t getUint32(const void* ptr, int offset)
{
   uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
   return (ptr8[0]<<24) | (ptr8[1]<<16) | (ptr8[2]<<8) | ptr8[3];
}

static std::string toString(int v)
{
   char buf[256];
   sprintf(buf, "%d", v);
   return buf;
}

// CRC16 from http://stackoverflow.com/questions/10564491/function-to-calculate-a-crc16-checksum
static unsigned short crc16(const unsigned char* data_p, unsigned char length){
   unsigned char x;
   unsigned short crc = 0xFFFF;
   
   while (length--){
      x = crc >> 8 ^ *data_p++;
      x ^= x>>4;
      crc = (crc << 8) ^ ((unsigned short)(x << 12)) ^ ((unsigned short)(x <<5)) ^ ((unsigned short)x);
   }
   return crc;
}

Alpha16Packet* Alpha16Packet::Unpack(const void*ptr, int bklen8)
{
   Alpha16Packet* p = new Alpha16Packet();
   /*
     ALPHA16 UDP packet data format from Bryerton: this is packet version 1:
     
     Date: Mon, 20 Jun 2016 15:07:08 -0700
     From bryerton@triumf.ca  Mon Jun 20 15:07:05 2016
     From: Bryerton Shaw <bryerton@triumf.ca>
     To: Konstantin Olchanski <olchansk@triumf.ca>
     Subject: Re: eta on udp data?
     
     Hi Konstantin,
     
     Just trying to iron out one last bug, it was (is?) locking up if I
     saturated the link, but I think I just resolved that! So we ve got all
     16 channels outputting pretty steadily up to 40kHz or so, if I reduce
     it to 3 channels, we can get 150+ kHz event rates per channel.
     
     I am going to add a checksum onto the packet structure but it looks as
     follows, broken down by BYTE offset
      
      0 Packet Type - Currently fixed at0x01
      1 Packet Version - Currently fixed at0x01
      2 Accepted Trigger MSB - Inside the firmware logic Accepted Trigger is unsigned 32bits, providing the lower 16bits here for useful as a dropped UDP packet check
      3 Accepted Trigger LSB
      4 MSB Hardware ID - Currently the lower 48 bits of the ArriaV ChipID. It will be the MAC address shortly however
      5 "" ""
      6 "" ""
      7 "" ""
      8 "" ""
      9 LSB Hardware ID
      10 Build Timestamp (UNIX timestamp, aka seconds since Jan 01, 1980 UTC)
      11 "" ""
      12 "" ""
      13 "" ""
      14 0x00
      15 0x00
      16 MSB Event Timestamp
      17 "" ""
      18 "" ""
      19 "" ""
      20 "" ""
      21 LSB Event Timestamp
      22 MSB Trigger Offset - Trigger Point in relation to the start of the waveform packet. Signed 32bit integer
      23
      24
      25 LSB Trigger Offset
      26 ModuleID - Logical Identifier for the Alpha16 board. unsigned byte
      27 [7] Channel Type - Either 0 or 1, for BScint or Anode. The MSB of this byte
      27 [6:0] Channel ID - Unsigned 7bits, identifies the ADC channel (0-15) used
      28 MSB Sample Count - Unsigned 16 bit value indicating the number of samples (1-N)
      29 LSB Sample Count
      30 MSB First Waveform Sample - Signed 16 bit value
      31 LSB First Waveform Sample
      ....
      30 + (SampleCount*2) MSB Checksum
      31 + (SampleCount*2) LSB Checksum
      
      I will give you the checksum details in a moment, I am just adding it in
      now. Most likely will be a crc16 based on 1+x^2+x^15+x^16 .
      The byte positions may not be ideal, but we can sort that out.
      
      Cheers,
      
      Bryerton
   */
   
   p->bankLength = bklen8;
   p->packetType = getUint8(ptr, 0);
   p->packetVersion = getUint8(ptr, 1);
   p->acceptedTrigger = getUint16(ptr, 2);
   p->hardwareId = getUint32(ptr, 4);
   p->buildTimestamp = getUint32(ptr, 10);
   //int zero = getUint16(ptr, 14);
   p->eventTimestamp = getUint32(ptr, 18);
   p->triggerOffset = getUint32(ptr, 22);
   p->moduleId = getUint8(ptr, 26);
   int chanX = getUint8(ptr, 27);
   p->channelType = chanX & 0x80;
   p->channelId = chanX & 0x7F;
   p->nsamples = getUint16(ptr, 28);
   p->checksum = getUint16(ptr, 30 + p->nsamples*2);
   p->length = 32 + p->nsamples*2;
   
   p->xcrc16 = crc16((const unsigned char*)ptr, 32 + p->nsamples*2);

   return p;
}

int Alpha16Packet::PacketType(const void*ptr, int bklen8)
{
   return getUint8(ptr, 0);
}

int Alpha16Packet::PacketVersion(const void*ptr, int bklen8)
{
   return getUint8(ptr, 1);
}

uint32_t Alpha16Packet::PacketTimestamp(const void*ptr, int bklen8)
{
   return getUint32(ptr, 18);
}

int Alpha16Packet::PacketChannel(const void*ptr, int bklen8)
{
   int chanX = getUint8(ptr, 27);
   //int channelType = chanX & 0x80;
   int channelId = chanX & 0x7F;
   return channelId;
}

void Alpha16Packet::Print() const
{
   printf("ALPHA16 data packet:\n");
   printf("  packet type:    0x%02x (%d)\n", packetType, packetType);
   printf("  packet version: 0x%02x (%d)\n", packetVersion, packetVersion);
   printf("  hwid:     0x%08x\n", hardwareId);
   printf("  buildts:  0x%08x\n", buildTimestamp);
   printf("  mod id:   0x%02x\n", moduleId);
   printf("  trig no:  0x%04x (%d)\n", acceptedTrigger, acceptedTrigger);
   printf("  event ts: 0x%08x (%d)\n", eventTimestamp, eventTimestamp);
   printf("  trig offset:   %d\n", triggerOffset);
   printf("  channel: type: %d, id: %d\n", channelType, channelId);
   printf("  nsamples: %d\n", nsamples);
   printf("  checksum: 0x%04x, computed checksum 0x%04x\n", checksum, xcrc16);
   printf("length: %d, bank length %d\n", length, bankLength);
};

void Alpha16Channel::Print() const
{
   printf("Alpha16Channel: bank %s, adc module %2d chan %2d, preamp pos %2d, wire %2d, tpc_wire %3d, first_bin %d, samples %d", bank.c_str(), adc_module, adc_chan, preamp_pos, preamp_wire, tpc_wire, first_bin, (int)adc_samples.size());
}

Alpha16Channel* Unpack(const char* bankname, int module, const Alpha16Packet* p, const void* bkptr, int bklen8)
{
   Alpha16Channel* c = new Alpha16Channel;

   c->bank = bankname;
   c->adc_module = module;
   if (p->channelType == 0) {
      c->adc_chan = p->channelId;
   } else {
      c->adc_chan = 16 + p->channelId;
   }
   c->preamp_pos  = -1;
   c->preamp_wire = -1;
   c->tpc_wire    = -1;
   c->first_bin   = 0;

   int nsamples = getUint16(bkptr, 28);

   c->adc_samples.reserve(nsamples);
   c->adc_samples.clear();
   
   for (int i=0; i<nsamples; i++) {
      unsigned v = getUint16(bkptr, 30 + i*2);
      // manual sign extension
      if (v & 0x8000)
         v |= 0xffff0000;
      c->adc_samples.push_back(v);
   }

   return c;
};

Alpha16Event::Alpha16Event() // ctor
{
}

Alpha16Event::~Alpha16Event() // dtor
{
}

void Alpha16Event::Print() const
{
   std::string e;
   e += toString(error);
   if (error) {
      e += " (";
      e += error_message;
      e += ")";
   }

   printf("AdcEvent %d, time %.6f, incr %.6f, complete %d, error %s, hits: %d, udp: %d", counter, time, timeIncr, complete, e.c_str(), (int)hits.size(), (int)udp.size());
}

static std::vector<std::string> split(const std::string& s, char seperator)
{
   std::vector<std::string> output;
   std::string::size_type prev_pos = 0, pos = 0;
   while((pos = s.find(seperator, pos)) != std::string::npos) {
      std::string substring( s.substr(prev_pos, pos-prev_pos) );
      output.push_back(substring);
      prev_pos = ++pos;
   }
   output.push_back(s.substr(prev_pos, pos-prev_pos));
   return output;
}

static int xatoi(const char* s)
{
   while (*s) {
      if (*s >= '0' && *s <= '9') {
         return atoi(s);
      }
      s++;
   }
   return 0;
}

void Alpha16Map::Init(const std::vector<std::string>& map)
{
   fNumChan = 0;

   for (unsigned i=0; i<map.size(); i++) {
      std::vector<std::string> s = split(map[i], ' ');
      if (s.size() != 4) {
         printf("invalid adc map entry %d: [%s]\n", i, map[i].c_str());
         abort();
         continue;
      }

      int imodule = xatoi(s[0].c_str());

      if (imodule < 1 || imodule > 20) {
         printf("invalid adc map entry %d: [%s], bad imodule %d\n", i, map[i].c_str(), imodule);
         abort();
      }

      int iconn = atoi(s[1].c_str());
      if (iconn < 0 || iconn > 2) {
         printf("invalid adc map entry %d: [%s], bad iconn %d\n", i, map[i].c_str(), iconn);
         abort();
      }

      int ipreamp = xatoi(s[3].c_str());
      if (s[3][0] == 'T') {
         ipreamp += 16;
         if (ipreamp < 16 || ipreamp > 31) {
            printf("invalid adc map entry %d: [%s], bad ipreamp %d\n", i, map[i].c_str(), ipreamp);
            abort();
         }
      } else if (s[3][0] == 'B') {
         ipreamp += 0;
         if (ipreamp < 0 || ipreamp > 15) {
            printf("invalid adc map entry %d: [%s], bad ipreamp %d\n", i, map[i].c_str(), ipreamp);
            abort();
         }
      } else if (s[3][0] < '0' || s[3][0] > '9') {
         printf("invalid adc map entry %d: [%s], bad preamp \"%s\"\n", i, map[i].c_str(), s[3].c_str());
         abort();
      }

      if (ipreamp < 0 || ipreamp > 31) {
         printf("invalid adc map entry %d: [%s], bad ipreamp %d\n", i, map[i].c_str(), ipreamp);
         abort();
      }

      printf("map %d: [%s] split [%s] [%s] [%s] [%s], module %2d, connector %1d, preamp %2d\n", i, map[i].c_str(), s[0].c_str(), s[1].c_str(), s[2].c_str(), s[3].c_str(), imodule, iconn, ipreamp);

      while ((int)fMap.size() <= imodule) {
         Alpha16MapEntry e;
         fMap.push_back(e);
      }

      fMap[imodule].module = imodule;
      if (iconn == 0) {
         fMap[imodule].preamp_0 = ipreamp;
         fNumChan += 16;
      } else if (iconn == 1) {
         fMap[imodule].preamp_1 = ipreamp;
         fNumChan += 16;
      } else if (iconn == 2) {
         fMap[imodule].preamp_2 = ipreamp;
         fNumChan += 16;
      } else {
         abort();
      }

      if (fFirstModule <= 0) {
         fFirstModule = imodule;
      }
   }
}

void Alpha16Map::Print() const
{
   printf("ADC map, %d modules, %d channels, first module %d:\n", (int)fMap.size(), fNumChan, fFirstModule);
   for (unsigned i=0; i<fMap.size(); i++) {
      printf(" slot %2d: module %2d preamps %3d and %3d %3d\n", i, fMap[i].module, fMap[i].preamp_0, fMap[i].preamp_1, fMap[i].preamp_2);
   }
}

Alpha16Asm::Alpha16Asm() // ctor
{
   Init();
}

// wire number -> adc16 channel
//static const int chanmap_top[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
static const int chanmap_top[] = { 7, 15, 6, 14, 5, 13, 4, 12, 3, 11, 2, 10, 1, 9, 0, 8 };
static const int chanmap_bot[] = { 8, 0, 9, 1, 10, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7 };

// adc16 channel -> wire number
static const int inv_chanmap_top[] = { 14, 12, 10, 8, 6, 4, 2, 0, 15, 13, 11, 9, 7, 5, 3, 1 };
static const int inv_chanmap_bot[] = { 1, 3, 5, 7, 9, 11, 13, 15, 0, 2, 4, 6, 8, 10, 12, 14 };

#if 0
// wire number -> adc32 channel, map by K.O.
static const int adc32_chanmap_ko[32] = {
   0,
   10,//1
   6,//2,
   8,//3,
   5,//9,//4,//5,//4,
   9,//1,//4,//9,//4,//5,
   2,//1,//9,//5,//4,//1,//2,//6,
   1,//2,//14,//7,
   14,//2,//3,//8,
   3,//4,//5,//4,//9,
   13,//4,//3,//2,//1,//10,
   15,//11,
   11,//20,//11,//15,//12,
   12,//4,//13,
   20,//18,//20,//24,//18,//25,//28,//25,//20,//7,//14,
   28,//20,//28,//16,//24,//20,//11,//20,//4,//25,//28,//18,//4,//12,//15,
   4,//28,//4,//28,//24,//16,//25,//18,//4,//20,//16,
   18,//4,//7,//4,//24,//28,//25,//7,//24,//17,
   7,//18,//7,//24,//4,//7,//25,//4,//18,
   25,//18,//24,//7,//4,//16,//20,//25,//19,
   16,//25,//18,//24,//7,//20,
   24,//25,//16,//4,//18,//28,//21,
   17,//30,//22,
   30,//17,//23,
   23,//7,//24,
   19,//25,
   21,//28, //26,
   26,//29,//22,//30, //27,
   27,//29,//26, //28,
   29,//27, //30, //29,
   22,//29, //30,
   31
};
#endif

// adc32 channel -> wire number, map from Keith Ong
// verified to be correct by K.O. by pulsing preamp inputs, one at a time
static const int inv_adc32_chanmap_top[32] = {
   1, 6, 3, 8,
   16+1, 5, 2, 16+2,
   4, 10, 0, 12,
   13, 9, 7, 11,
   16+3, 16+5, 16+6, 16+9,
   15, 16+11, 16+14, 16+8,
   16+4, 14, 16+10, 16+12,
   16+0, 16+13, 16+7, 16+15
};

#if 0
static int inv_adc32_chanmap_ko[32];
#endif

void Alpha16Asm::Init()
{
   // construct or check the inverted adc16 map
   
   for (int xchan=0; xchan<16; xchan++) {
      int ychan = -1;
      for (int i=0; i<16; i++)
         if (chanmap_bot[i] == xchan) {
            ychan = i;
            break;
         }
      assert(inv_chanmap_bot[xchan] == ychan);

      ychan = -1;
      for (int i=0; i<16; i++)
         if (chanmap_top[i] == xchan) {
            ychan = i;
            break;
         }
      assert(inv_chanmap_top[xchan] == ychan);
   }

#if 0
   // construct the inverted adc32 map

   for (int xchan=0; xchan<32; xchan++) {
      inv_adc32_chanmap[xchan] = -1;
   }

   for (int xchan=0; xchan<32; xchan++) {
      int ychan = -1;
      for (int i=0; i<32; i++)
         if (adc32_chanmap[i] == xchan) {
            ychan = i;
            break;
         }
      assert(inv_adc32_chanmap[xchan] == -1);
      inv_adc32_chanmap[xchan] = ychan;
   }
#endif

#if 0
   printf("inv_adc32_chanmap: ");
   for (int xchan=0; xchan<32; xchan++) {
      printf(" %d,", inv_adc32_chanmap[xchan]);
   }
   printf("\n");
#endif

#if 0
   printf("inv_adc32_chanmap: ");
   for (int xchan=0; xchan<32; xchan++) {
      int wire = inv_adc32_chanmap[xchan];
      int conn = wire/16;
      int wire16 = wire%16;
      int chan16 = chanmap_top[wire16];
      printf(" %d/%d/%d,", wire, conn, chan16);
   }
   printf("\n");
#endif

#if 0
   printf("inv_chanmap_bot: ");
   for (int xchan=0; xchan<16; xchan++) {
      printf(" %d,", inv_chanmap_bot[xchan]);
   }
   printf("\n");

   printf("inv_chanmap_top: ");
   for (int xchan=0; xchan<16; xchan++) {
      printf(" %d,", inv_chanmap_top[xchan]);
   }
   printf("\n");
#endif
}

void Alpha16Asm::AddChannel(Alpha16Event* e, Alpha16Packet* p, Alpha16Channel* c)
{
   int imodule = c->adc_module;
   if (imodule > 0 && imodule < (int)fMap.fMap.size() && fMap.fMap[imodule].module == imodule) {
      int pos0 = fMap.fMap[imodule].preamp_0;
      int pos1 = fMap.fMap[imodule].preamp_1;
      int pos2 = fMap.fMap[imodule].preamp_2;
      int xchan = -1;
      if (c->adc_chan < 16) {
         c->preamp_pos = pos0;
         xchan = c->adc_chan;
      } else if (c->adc_chan < 48) {
         int ichan = c->adc_chan-16; // 32ch ADC channel number
         //int zchan = inv_adc32_chanmap[ichan];
         //int zchan = ichan;
         int zchan = inv_adc32_chanmap_top[ichan];
         if (zchan < 16) {
            c->preamp_pos = pos1;
            xchan = zchan;
         } else if (zchan < 32) {
            c->preamp_pos = pos2;
            xchan = zchan-16;
         } else {
            abort(); // cannot happen
         }
         c->preamp_wire = xchan;
         c->tpc_wire = c->preamp_pos*16 + c->preamp_wire;
         if (0 && c->tpc_wire == 288) {
            printf("remap: module %d, adc_chan %d, ichan %d, zchan %d, xchan %d, preamp pos %d, wire %d, tpc_wire %d\n",
                   imodule,
                   c->adc_chan,
                   ichan,
                   zchan,
                   xchan,
                   c->preamp_pos,
                   c->preamp_wire,
                   c->tpc_wire
                   );
         }
         //if (ichan == 0) {
         //   c->tpc_wire = 1;
         //}
         xchan = -1;
      }

      if (xchan >= 0 && xchan < 16) {
         if (c->preamp_pos < 16) {
            // bot
            c->preamp_wire = inv_chanmap_bot[xchan];
         } else {
            // top
            c->preamp_wire = inv_chanmap_top[xchan];
         }
         c->tpc_wire = c->preamp_pos*16 + c->preamp_wire;
      }
   }
   
   //printf("AddBank: channel: "); c->Print(); printf("\n");

   e->udp.push_back(p);
   e->hits.push_back(c);
}

#if 0
void Alpha16EVB::AddBank(Alpha16Event* e, Alpha16Packet* p, Alpha16Channel* c)
{
   int ymodule = -1;
   int top_bot = 0;
   
   for (unsigned x=0; x<fConfModMap.size(); x++)
      if (fConfModMap[x] == xmodule) {
         ymodule = x;
         top_bot = 1;
         break;
      } else if (fConfModMap[x] == -xmodule) {
         ymodule = x;
         top_bot = -1;
         break;
      }
   
   if (ymodule < 0)
      return;

   int xchan = Alpha16Packet::PacketChannel(bkptr, bklen);

   if (xchan<0 || xchan >= NUM_CHAN_ALPHA16) {
      // invalid channel number
      e->error_message += "| invalid channel number ";
      e->error_message += toString(xchan);
      e->error = true;
      return;
   }

   int ychan = -1; 

   if (top_bot < 0)
      for (int i=0; i<16; i++)
         if (chanmap_bot[i] == xchan) {
            ychan = i;
            break;
         }

   if (top_bot > 0)
      for (int i=0; i<16; i++)
         if (chanmap_top[i] == xchan) {
            ychan = i;
            break;
         }

   assert(ychan>=0);

   int chan = ymodule*NUM_CHAN_ALPHA16 + ychan;

   if (e->udpPresent[chan]) {
      // duplicate udp packet
      e->error_message += "| duplicate udp channel";
      e->error = true;
      return;
   }
   
   e->udpPacket[chan].Unpack(bkptr, bklen);
   e->waveform[chan].Unpack(bkptr, bklen);
   e->udpPresent[chan] = true;
   e->udpEventTs[chan] = e->udpPacket[chan].eventTimestamp;

   e->udpEventTsIncr[chan] = e->udpEventTs[chan] - fLastUdpEventTs[chan];
   fLastUdpEventTs[chan] = e->udpEventTs[chan];

   e->numChan++;
}
#endif

Alpha16Event* Alpha16Asm::NewEvent()
{
   Alpha16Event *e = new Alpha16Event();
   fEventCount++;
   e->counter = fEventCount;
   return e;
}

//const double TSCLK = 0.1; // GHz

void Alpha16Asm::CheckEvent(Alpha16Event* e)
{
   //e->Print(); printf("\n");

   if ((int)e->udp.size() != fMap.fNumChan) {
      e->error = true;
      e->error_message = "incomplete";
      e->complete = false;
      return;
   }

   e->complete = true;
   
   assert(e->udp.size() == e->hits.size());

   if (e->udp.size() > 0) {
      int i=-1;
      for (unsigned j=0; j<e->udp.size(); j++) {
         if (e->hits[j]->adc_module == fMap.fFirstModule) {
            i = j;
            break;
         }
      }
      assert(i>=0);
      //int module = e->hits[i]->adc_module;
      uint32_t ets = e->udp[i]->eventTimestamp;

      bool wrap = false;
      if (ets < fLastEventTs) {
         wrap = true;
         //printf("wrap!\n");
      }

      if (wrap) {
         fTsEpoch++;
      }

      double eventTime = ets/fTsFreq + fTsEpoch*(2.0*0x80000000/fTsFreq);

      if (e->counter <= 1) {
         fFirstEventTime = eventTime;
      }

      eventTime -= fFirstEventTime;

      e->time = eventTime;
      e->timeIncr = eventTime - fLastEventTime;

      fLastEventTs = ets;
      fLastEventTime = eventTime;
   }
   
#if 0
   int count = 0;
   unsigned num_samples = 0;

   for (int i=0; i<MAX_ALPHA16 * NUM_CHAN_ALPHA16; i++) {
      if (e->udpPresent[i]) {
         count++;
         if (num_samples == 0)
            num_samples = e->waveform[i].size();
         if (e->waveform[i].size() != num_samples) {
            // wrong number of ADC samples
            e->error_message += "| wrong number of ADC samples ";
            e->error_message += toString(e->waveform[i].size());
            e->error_message += " should be ";
            e->error_message += toString(num_samples);
            e->error = true;
         }
      }
   }

   if (fConfNumSamples == 0)
      fConfNumSamples = num_samples;

   if ((int)num_samples != fConfNumSamples) {
      // wrong number of ADC samples
      e->error_message += "| wrong number of ADC samples ";
      e->error_message += toString(num_samples);
      e->error_message += " should be fConfNumSamples ";
      e->error_message += toString(fConfNumSamples);
      e->error = true;
   }

   if ((count == fConfNumChan) && (e->numChan == fConfNumChan)) {
      e->complete = true;
   }

   if (!fHaveEventTs && e->complete) {
      fHaveEventTs = true;
      // set timestamp offsets from the first complete event
      for (int i=0; i<MAX_ALPHA16 * NUM_CHAN_ALPHA16; i++) {
         if (e->udpPresent[i]) {
            fFirstEventTs[i] = e->udpEventTs[i];
            //printf("XXX %d -> 0x%08x\n", i, fFirstEventTs[i]);
         }
      }
   }

   if (fHaveEventTs && e->complete) {
      uint32_t ets = 0;
      // check timestamps
      for (int i=0; i<MAX_ALPHA16 * NUM_CHAN_ALPHA16; i++) {
         if (e->udpPresent[i]) {
            uint32_t ts = e->udpEventTs[i] - fFirstEventTs[i];
            if (ets == 0)
               ets = ts;
            if (ts != ets && ts != ets+1 && ts+1 != ets) {
               double ts125 = ts*100.0/125.0;
               ts = ts125;
               if (ts != ets && ts != ets+1 && ts+1 != ets) {
                  printf("ts mismatch %d: 0x%08x vs 0x%08x diff %d\n", i, ts, ets, ts-ets);
                  e->error_message += "| timestamp mismatch ";
                  e->error = true;
               }
            }
         }
      }

      bool wrap = false;
      if (ets < fLastEventTs) {
         wrap = true;
         //printf("wrap!\n");
      }

      if (wrap) {
         fTsEpoch++;
      }

      double eventTime = ets/TSCLK + fTsEpoch*(2.0*0x80000000/TSCLK);

      e->eventTime = eventTime;
      e->prevEventTime = fLastEventTime;

      e->time = e->eventTime * 1e-9; // ns to sec
      e->timeIncr = (e->eventTime-e->prevEventTime) * 1e-9; // ns to sec

      fLastEventTs = ets;
      fLastEventTime = eventTime;
   } else if (e->numChan <= 16) { // unsynchronized A16 event
      e->eventTime = -1;
      e->prevEventTime = -1;
   } else { // incomplete event
      e->eventTime = -2;
      e->prevEventTime = -2;
   }
#endif
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
