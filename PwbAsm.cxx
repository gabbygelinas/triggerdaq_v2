//
// Unpacking PWB data
// K.Olchanski
//

#include "PwbAsm.h"

#include <stdio.h> // NULL, printf()
//#include <math.h> // fabs()
#include <assert.h> // assert()

PwbUdpPacket::PwbUdpPacket(const char* ptr, int size) // ctor
{
   const uint32_t* p32 = (const uint32_t*)ptr;
   int n32 = size/4;
   
   fError = false;
   fPacketSize = size;

   if (n32 < 6) {
      fError = true;
      return;
   }

   MYSTERY     = p32[0];
   PKT_SEQ     = p32[1];
   CHANNEL_SEQ = (p32[2] >>  0) & 0xFFFF;
   CHANNEL_ID  = (p32[2] >> 16) & 0xFF;
   FLAGS       = (p32[2] >> 24) & 0xFF;
   CHUNK_ID    = (p32[3] >>  0) & 0xFFFF;
   CHUNK_LEN   = (p32[3] >> 16) & 0xFFFF;
   HEADER_CRC  = p32[4];
   start_of_payload = 5*4;
   end_of_payload = start_of_payload + CHUNK_LEN;

   if (end_of_payload + 4 > fPacketSize) {
      fError = true;
      return;
   }
   
   payload_crc = p32[end_of_payload/4];
}

void PwbUdpPacket::Print() const
{
   printf("PwbUdpPacket: M 0x%08x, PKT_SEQ 0x%08x, CHAN SEQ 0x%04x, ID 0x%02x, FLAGS 0x%02x, CHUNK ID 0x%04x, LEN 0x%04x, CRC 0x%08x, bank bytes %d, end of payload %d, CRC 0x%08x\n",
          MYSTERY,
          PKT_SEQ,
          CHANNEL_SEQ,
          CHANNEL_ID,
          FLAGS,
          CHUNK_ID,
          CHUNK_LEN,
          HEADER_CRC,
          fPacketSize,
          end_of_payload,
          payload_crc);
};

PwbEventHeader::PwbEventHeader(const char* ptr, int size)
{
   const uint32_t* p32 = (const uint32_t*)ptr;
   int nw = size/4;

   fError = false;
   
   FormatRevision  = (p32[5]>> 0) & 0xFF;
   ScaId           = (p32[5]>> 8) & 0xFF;
   CompressionType = (p32[5]>>16) & 0xFF;
   TriggerSource   = (p32[5]>>24) & 0xFF;
   
   HardwareId1 = p32[6];
   
   HardwareId2 = (p32[7]>> 0) & 0xFFFF;
   TriggerDelay     = (p32[7]>>16) & 0xFFFF;
   
   // NB timestamp clock is 125 MHz
   
   TriggerTimestamp1 = p32[8];
   
   TriggerTimestamp2 = (p32[9]>> 0) & 0xFFFF;
   Reserved1         = (p32[9]>>16) & 0xFFFF;
   
   ScaLastCell = (p32[10]>> 0) & 0xFFFF;
   ScaSamples  = (p32[10]>>16) & 0xFFFF;
   
   ScaChannelsSent1 = p32[11];
   ScaChannelsSent2 = p32[12];
   
   ScaChannelsSent3 = (p32[13]>> 0) & 0xFF;
   ScaChannelsThreshold1 = (p32[13]>> 8) & 0xFFFFFF;
   
   ScaChannelsThreshold2 = p32[14];
   
   ScaChannelsThreshold3 = p32[15] & 0xFFFF;
   Reserved2             = (p32[15]>>16) & 0xFFFF;

   start_of_data = 16*4;
}

void PwbEventHeader::Print() const
{
   if (1) {
      printf("PwbEventHeader: F 0x%02x, Sca 0x%02x, C 0x%02x, T 0x%02x, H 0x%08x, 0x%04x, Delay 0x%04x, TS 0x%08x, 0x%04x, R1 0x%04x, SCA LastCell 0x%04x, Samples 0x%04x, Sent 0x%08x 0x%08x 0x%08x, Thr 0x%08x 0x%08x 0x%08x, R2 0x%04x\n",
             FormatRevision,
             ScaId,
             CompressionType,
             TriggerSource,
             HardwareId1, HardwareId2,
             TriggerDelay,
             TriggerTimestamp1, TriggerTimestamp2,
             Reserved1,
             ScaLastCell,
             ScaSamples,
             ScaChannelsSent1,
             ScaChannelsSent2,
             ScaChannelsSent3,
             ScaChannelsThreshold1,
             ScaChannelsThreshold2,
             ScaChannelsThreshold3,
             Reserved2);
   }

   if (0) {
      printf("S Sca 0x%02x, TS 0x%08x, 0x%04x, SCA LastCell 0x%04x, Samples 0x%04x\n",
             ScaId,
             TriggerTimestamp1, TriggerTimestamp2,
             ScaLastCell,
             ScaSamples);
   }
}

#define PWB_CA_ST_INIT 0
#define PWB_CA_ST_DATA 1
#define PWB_CA_ST_LAST 2

PwbChannelAsm::PwbChannelAsm(int module, int sca)
{
   fModule = module;
   fSca = sca;
   Reset();
}

void PwbChannelAsm::Reset()
{
   fLast_CHANNEL_SEQ = 0;
   fState = 0;
   fSaveChannel = 0;
   fSaveSamples = 0;
   fSaveNw = 0;
   fSavePos = 0;
}

void PwbChannelAsm::AddSamples(int channel, const uint16_t* samples, int count)
{
   printf("pwb module %d, sca %d, channel %d, add %d samples\n", fModule, fSca, channel, count);
}

void PwbChannelAsm::CopyData(const uint16_t* s, const uint16_t* e)
{
   const uint16_t* p = s;

   while (1) {
      int r = e-p;
      if (r < 2) {
         printf("need to en-buffer header!\n");
         break;
      }
      if (p[0] == 0xCCCC && p[1] == 0xCCCC) {
         printf("CopyData: unexpected 0xCCCC words, r %d!\n", r);
         break;
      }
      int channel = p[0];
      int samples = p[1];
      int nw = samples;
      if (samples&1)
         nw+=1;
      if (nw <= 0) {
         printf("invalid word counter!\n");
         break;
      }
      printf("adc samples ptr %d, end %d, r %d, channel %d, samples %d, nw %d\n", (int)(p-s), (int)(e-s), r, channel, samples, nw);
      p += 2;
      r = e-p;
      int h = nw;
      bool truncated = false;
      if (nw > r) {
         h = r;
         printf("split data nw %d, with %d here, %d to follow\n", nw, h, nw-h);
         fSaveChannel = channel;
         fSaveSamples = samples;
         fSaveNw = nw;
         fSavePos = h;
         truncated = true;
      }
      AddSamples(channel, p, h);
      p += h;
      if (truncated) {
         assert(p == e);
         break;
      }
      if (p == e) {
         break;
      }
   }
}

void PwbChannelAsm::BeginData(const char* ptr, int size, int start_of_data, int end_of_data)
{
   const uint16_t* s = (const uint16_t*)(ptr+start_of_data);
   const uint16_t* e = (const uint16_t*)(ptr+end_of_data);
   CopyData(s, e);
}

void PwbChannelAsm::AddData(const char* ptr, int size, int start_of_data, int end_of_data)
{
   const uint16_t* s = (const uint16_t*)(ptr+start_of_data);
   const uint16_t* e = (const uint16_t*)(ptr+end_of_data);
   const uint16_t* p = s;

   if (fSaveNw > 0) {
      int r = (e-p);
      
      int h = fSaveNw - fSavePos;
      bool truncated = false;
      
      if (h > r) {
         h = r;
         truncated = true;
      }
      
      printf("AddData: save channel %d, samples %d, nw %d, pos %d, remaining %d, have %d, truncated %d\n", fSaveChannel, fSaveSamples, fSaveNw, fSavePos, r, h, truncated);
      
      AddSamples(fSaveChannel, p, h);
      
      if (!truncated) {
         fSaveChannel = 0;
         fSaveSamples = 0;
         fSaveNw = 0;
         fSavePos = 0;
      }
      
      p += h;
   }

   if (p < e) {
      CopyData(p, e);
   } else if (p == e) {
      // good, no more data in this packet
   } else {
      printf("AddData: error!\n");
      assert(!"this cannot happen");
   }
}

void PwbChannelAsm::EndData()
{
   if (fSaveNw > 0) {
      printf("EndData: missing some data at the end\n");
   } else {
      printf("EndData: ok!\n");
   }

   fSaveChannel = 0;
   fSaveSamples = 0;
   fSaveNw = 0;
   fSavePos = 0;
}

void PwbChannelAsm::AddPacket(PwbUdpPacket* udp, const char* ptr, int size)
{
   if (fLast_CHANNEL_SEQ == 0) {
      fLast_CHANNEL_SEQ = udp->CHANNEL_SEQ;
   } else if (udp->CHANNEL_SEQ != fLast_CHANNEL_SEQ + 1) {
      printf("PwbChannelAsm::AddPacket(): misordered or lost UDP packet: CHANNEL_SEQ jump 0x%08x to 0x%08x\n", fLast_CHANNEL_SEQ, udp->CHANNEL_SEQ);
      fLast_CHANNEL_SEQ = udp->CHANNEL_SEQ;
   } else {
      fLast_CHANNEL_SEQ++;
   }

   if (fState == PWB_CA_ST_INIT || fState == PWB_CA_ST_LAST) {
      PwbEventHeader* eh = new PwbEventHeader(ptr, size);
      eh->Print();
      delete eh;
      fState = PWB_CA_ST_DATA;
      BeginData(ptr, size, eh->start_of_data, udp->end_of_payload);
   } else if (fState == PWB_CA_ST_DATA) {
      AddData(ptr, size, udp->start_of_payload, udp->end_of_payload);
      if (udp->FLAGS & 1) {
         fState = PWB_CA_ST_LAST;
         EndData();
      }
   } else {
      printf("PwbChannelAsm::AddPacket(): invalid state %d\n", fState);
   }
}

PwbModuleAsm::PwbModuleAsm(int module)
{
   fModule = module;
   Reset();
}

void PwbModuleAsm::Reset()
{
   fLast_PKT_SEQ = 0;
   for (unsigned i=0; i<fChannels.size(); i++) {
      if (fChannels[i])
         fChannels[i]->Reset();
   }
}

void PwbModuleAsm::AddPacket(const char* ptr, int size)
{
   PwbUdpPacket* udp = new PwbUdpPacket(ptr, size);
   udp->Print();

   if (fLast_PKT_SEQ == 0) {
      fLast_PKT_SEQ = udp->PKT_SEQ;
   } else if (udp->PKT_SEQ != fLast_PKT_SEQ + 1) {
      printf("PwbModuleAsm::AddPacket(): misordered or lost UDP packet: PKT_SEQ jump 0x%08x to 0x%08x\n", fLast_PKT_SEQ, udp->PKT_SEQ);
      fLast_PKT_SEQ = udp->PKT_SEQ;
   } else {
      fLast_PKT_SEQ++;
   }

   int s = udp->CHANNEL_ID;

   if (s < 0 || s > 4) {
      printf("PwbModuleAsm::AddPacket(): invalid CHANNEL_ID 0x%08x\n", udp->CHANNEL_ID);
   } else {
      while (s >= fChannels.size()) {
         fChannels.push_back(NULL);
      }

      if (!fChannels[s]) {
         fChannels[s] = new PwbChannelAsm(fModule, s);
      }

      fChannels[s]->AddPacket(udp, ptr, size);
   }
   
   delete udp;
}

PwbAsm::PwbAsm() // ctor
{
   // empty
}

PwbAsm::~PwbAsm() // dtor
{
   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i]) {
         delete fModules[i];
         fModules[i] = NULL;
      }
   }
}

void PwbAsm::AddPacket(int module, const char* ptr, int size)
{
   while (module >= fModules.size()) {
      fModules.push_back(NULL);
   }

   if (!fModules[module]) {
      fModules[module] = new PwbModuleAsm(module);
   }

   fModules[module]->AddPacket(ptr, size);
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


