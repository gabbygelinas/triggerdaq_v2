//
// Unpacking PWB data
// K.Olchanski
//

#ifndef PwbAsm_H
#define PwbAsm_H

#include <stdint.h> // uint32_t & co
#include <vector>

#include "Feam.h"

struct PwbUdpPacket
{
public:
   bool fError = false;
   int fPacketSize = 0;

public: // packet data
   uint32_t DEVICE_ID;
   uint32_t PKT_SEQ;
   uint16_t CHANNEL_SEQ;
   uint32_t CHANNEL_ID;
   uint32_t FLAGS;
   uint32_t CHUNK_ID;
   uint32_t CHUNK_LEN;
   uint32_t HEADER_CRC;
   uint32_t start_of_payload;
   uint32_t end_of_payload;
   uint32_t payload_crc;

public:
   PwbUdpPacket(const char* ptr, int size); // ctor
   void Print() const;
};
   
struct PwbEventHeader
{
public:
   bool fError = false;

public:
   int FormatRevision;
   int ScaId;
   int CompressionType;
   int TriggerSource;
   uint32_t HardwareId1;
   uint32_t HardwareId2;
   int TriggerDelay;
   uint32_t TriggerTimestamp1;
   uint32_t TriggerTimestamp2;
   uint32_t Reserved1;
   int ScaLastCell;
   int ScaSamples;
   uint32_t ScaChannelsSent1;
   uint32_t ScaChannelsSent2;
   uint32_t ScaChannelsSent3;
   uint32_t ScaChannelsThreshold1;
   uint32_t ScaChannelsThreshold2;
   uint32_t ScaChannelsThreshold3;
   uint32_t Reserved2;
   int start_of_data;

public:
   PwbEventHeader(const char* ptr, int size);
   void Print() const;
};

class PwbChannelAsm
{
public:
   int fModule = 0;
   int fColumn = 0;
   int fRing = 0;
   int fSca = 0;
   
public: // state
   uint16_t fLast_CHANNEL_SEQ = 0;
   int fState = 0;
   int fCountErrors = 0;
   uint32_t fTs = 0;

public: // configuration
   bool fTrace = false;

public:
   int fSaveChannel;
   int fSaveSamples;
   int fSaveNw;
   int fSavePos;

public: // channel currently filled
   FeamChannel* fCurrent = NULL;

public: // output
   std::vector<FeamChannel*> fOutput;

public:
   PwbChannelAsm(int module, int column, int ring, int sca); // ctor
   ~PwbChannelAsm(); // dtor
   
   void Reset();
   void AddPacket(PwbUdpPacket* udp, const char* ptr, int size);

   void BeginData(const char* ptr, int size, int start_of_data, int end_of_data, uint32_t ts);
   void AddData(const char* ptr, int size, int start_of_data, int end_of_data);
   void EndData();

   void CopyData(const uint16_t* s, const uint16_t* e);

   void AddSamples(int channel, const uint16_t* samples, int count);

   bool CheckComplete() const;

   void BuildEvent(FeamEvent* e);
};

class PwbModuleAsm
{
public:
   int fModule = 0;
   int fColumn = 0;
   int fRing = 0;
   std::vector<PwbChannelAsm*> fChannels;

public: // state
   int fCountErrors = 0;
   uint32_t fLast_PKT_SEQ = 0;
   uint32_t fTs = 0;
   double   fTime = 0;
   double   fTimeIncr = 0;

public: // timestamps
   uint32_t fTsFirstEvent = 0;
   uint32_t fTsLastEvent = 0;
   int      fTsEpoch = 0;
   double   fTimeFirstEvent = 0;
   double   fTimeLastEvent = 0;

public: // configuration
   bool fTrace = false;

public:
   PwbModuleAsm(int module, int column, int ring); // ctor
   ~PwbModuleAsm(); // dtor
   void Reset();
   void AddPacket(const char* ptr, int size);
   bool CheckComplete() const;
   void BuildEvent(FeamEvent* e);
};

class PwbAsm
{
public:
   std::vector<PwbModuleAsm*> fModules;

public: // state
   int fCounter = 0;
   double fLastTime = 0;

public:
   PwbAsm();
   ~PwbAsm();

public:
   void Reset();
   void AddPacket(int module, int column, int ring, const char* ptr, int size);
   bool CheckComplete() const;
   void BuildEvent(FeamEvent* e);
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


