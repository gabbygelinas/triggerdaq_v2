#ifndef ALPHA16H
#define ALPHA16H

#include "Waveform.h"

uint8_t getUint8(const void* ptr, int offset)
{
  return *(uint8_t*)(((char*)ptr)+offset);
}

uint16_t getUint16(const void* ptr, int offset)
{
  uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
  return ((ptr8[0]<<8) | ptr8[1]);
}

uint32_t getUint32(const void* ptr, int offset)
{
  uint8_t *ptr8 = (uint8_t*)(((char*)ptr)+offset);
  return (ptr8[0]<<24) | (ptr8[1]<<16) | (ptr8[2]<<8) | ptr8[3];
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

struct Alpha16Packet
{
  int bankLength;
  int packetType;
  int packetVersion;
  int acceptedTrigger;
  uint32_t hardwareId;
  uint32_t buildTimestamp;
  uint32_t eventTimestamp;
  uint32_t triggerOffset;
  int moduleId;
  int channelType;
  int channelId;
  int nsamples;
  int checksum;
  int length;
  int xcrc16;

  Waveform* w;

  Alpha16Packet(const void*ptr, int bklen8) // ctor
  {
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

    bankLength = bklen8;
    packetType = getUint8(ptr, 0);
    packetVersion = getUint8(ptr, 1);
    acceptedTrigger = getUint16(ptr, 2);
    hardwareId = getUint32(ptr, 4);
    buildTimestamp = getUint32(ptr, 10);
    //int zero = getUint16(ptr, 14);
    eventTimestamp = getUint32(ptr, 18);
    triggerOffset = getUint32(ptr, 22);
    moduleId = getUint8(ptr, 26);
    int chanX = getUint8(ptr, 27);
    channelType = chanX & 0x80;
    channelId = chanX & 0x7F;
    nsamples = getUint16(ptr, 28);
    checksum = getUint16(ptr, 30 + nsamples*2);
    length = 32 + nsamples*2;

    xcrc16 = crc16((const unsigned char*)ptr, 32 + nsamples*2);

    w = new Waveform(nsamples);

    //printf("Unpacking: "); Print();
    for (int i=0; i<nsamples; i++) {
      int16_t v = getUint16(ptr, 30 + i*2);
      //printf("sample %d: 0x%02x (%d)\n", i, v, v);
      w->samples[i] = v;
    }

    //Print();
  };

  ~Alpha16Packet() // dtor
  {
    if (w)
      delete w;
    w = NULL;
  }

  void Print() const
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
};

struct Alpha16Event
{
   uint16_t eventNo;
   uint32_t eventTs; // syncronized timestamp
   uint32_t prevEventTs; // same from previous event
   uint32_t udpEventTs;  // timestamp from udp packet
   int numChan;
   bool complete;
   Alpha16Packet* chan[16];
   
   Alpha16Event() // ctor
   {
      eventNo = 0;
      eventTs = 0;
      prevEventTs = 0;
      udpEventTs = 0;
      numChan = 0;
      complete = false;
      for (int i=0; i<16; i++)
         chan[i] = NULL;
   }
   
   ~Alpha16Event() // dtor
   {
      for (int i=0; i<16; i++) {
         if (chan[i])
            delete chan[i];
         chan[i] = NULL;
      }
   }

   void Print() const
   {
      printf("ALPHA16 event: %d, ts 0x%08x (incr 0x%08x), channels: %d, complete %d\n", eventNo, eventTs, eventTs - prevEventTs, numChan, complete);
   }
};

#include <deque>

struct Alpha16EVB
{
  int fEventCount; // event counter
  uint32_t fFirstEventTs; // udp timestamp of first event
  uint32_t fLastEventTs;  // udp timestamp of last seen event
  std::deque<Alpha16Event*> fEvents;

  Alpha16EVB() // ctor
  {
    Reset();
  }

  void Reset();
  void AddPacket(Alpha16Packet* p);
  Alpha16Event* GetNextEvent();

  static bool Match(const Alpha16Event* e, const Alpha16Packet* p);
  static bool AddToEvent(Alpha16Event* e, Alpha16Packet* p);
};

void Alpha16EVB::Reset()
{
  fEventCount = 0;
  fFirstEventTs = 0;
  fLastEventTs = 0;

  while (fEvents.size() > 0) {
    Alpha16Event* e = fEvents.back();
    fEvents.pop_back();
    if (e)
      delete e;
  }
}

bool Alpha16EVB::Match(const Alpha16Event* e, const Alpha16Packet* p)
{
  return (e->udpEventTs == p->eventTimestamp);
}

bool Alpha16EVB::AddToEvent(Alpha16Event* e, Alpha16Packet* p)
{
  int chan = p->channelId;
  if (chan < 0)
    return false;
  if (chan >= 16)
    return false;
  if (e->chan[chan] != NULL)
    return false; // check for duplicate channel
  e->chan[chan] = p;
  e->numChan++;
  if (e->numChan == 16)
    e->complete = true;
  //printf("Event %d ts 0x%08x added packet chan %2d, no %d, ts 0x%08x, have %2d, complete %d\n", e->eventNo, e->eventTs, chan, p->acceptedTrigger, p->eventTimestamp, e->numChan, e->complete);
  //if (p->triggerOffset != 15)
  //p->Print();
  return true;
}

void Alpha16EVB::AddPacket(Alpha16Packet*p )
{
  // loop over buffered events, look for match

  for (unsigned i=0; i<fEvents.size(); i++) {
    if (Match(fEvents[i], p)) {
      bool ok = AddToEvent(fEvents[i], p);
      if (!ok)
	delete p;
      return;
    }
  }

  Alpha16Event *e = new Alpha16Event();
  fEventCount++;
  e->eventNo = fEventCount;
  e->udpEventTs = p->eventTimestamp;
  if (fFirstEventTs == 0) {
    fFirstEventTs = p->eventTimestamp;
    fLastEventTs = fFirstEventTs;
  }
  e->eventTs = p->eventTimestamp - fFirstEventTs;
  e->prevEventTs = fLastEventTs - fFirstEventTs;

  fLastEventTs = p->eventTimestamp;

  //printf("Event %d ts 0x%08x: new!\n", e->eventNo, e->eventTs);

  bool ok = AddToEvent(e, p);
  if (!ok)
    delete p;
  p = NULL;

  fEvents.push_front(e);
}

Alpha16Event* Alpha16EVB::GetNextEvent()
{
  if (fEvents.size() < 1)
    return NULL;

  Alpha16Event* e = fEvents.back();

  if (e->complete || fEvents.size() > 10) {
    fEvents.pop_back();
    return e;
  }
  
  return NULL;
}

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
