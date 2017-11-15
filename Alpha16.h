//
// ALPHA-g TPC
//
// GRIF-16/ALPHA-16 ADC functions
//
// Class definitions
//

#ifndef ALPHA16H
#define ALPHA16H

#include <stdint.h>
#include <vector>
#include <string>

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

   static int PacketType(const void*ptr, int bklen8);
   static int PacketVersion(const void*ptr, int bklen8);
   static uint32_t PacketTimestamp(const void*ptr, int bklen8);
   static int PacketChannel(const void*ptr, int bklen8);
   static Alpha16Packet* Unpack(const void* bkptr, int bklen8);
   void Print() const;
};

struct Alpha16Channel
{
   std::string bank;
   int adc_module = -1;  // ADC module mod1..mod20
   int adc_chan = -1;    // ADC channel: 0..15: GRIF16 onboard 100MHz, 16..47: FMC-ADC32 62.5MHz
   int preamp_pos = -1;  // preamp position: 0..15: B0..B15, 16..31: T0..T15
   int preamp_wire = -1; // preamp wire number 0..15
   int tpc_wire = -1;    // TPC anode wire 0..127 bottom, 128..255 top

   int first_bin = 0; /* usually 0 */
   std::vector<int> adc_samples;

   void Print() const;
};

#if 0
struct Alpha16ModuleConfig
{
   double adc16_ts_freq = 0;
   double adc32_ts_freq = 0;
   
   int adc16_num_samples = 0;
   int adc32_num_samples = 0;

   int adc16_preamp_pos = 0;
   int adc32a_preamp_pos = 0;
   int adc32b_preamp_pos = 0;
};
#endif

struct Alpha16MapEntry
{
   int module = -1;   // ADC module number 1..20
   int preamp_0 = -1; // preamp connector 0 ch 0..15, 100MHz ADC
   int preamp_1 = -1; // preamp connector 1 ch 16..31, 62.5MHz ADC
   int preamp_2 = -1; // preamp connector 2 ch 32..47, 62.5MHz ADC
};

class Alpha16Map
{
 public:
   std::vector<Alpha16MapEntry> fMap;
 public:
   void Init(const std::vector<std::string>& map);
   void Print() const;
};

Alpha16Channel* Unpack(const char* bankname, int module, const Alpha16Packet* p, const void* bkptr, int bklen8);

#define NUM_CHAN_ALPHA16 16

struct Alpha16Event
{
   bool complete = false; // event is complete
   bool error = false;    // event has an error
   std::string error_message; // error message
   int  counter = 0;      // event sequential counter
   double time = 0;       // event time, sec
   double timeIncr = 0;   // time from previous event, sec

   std::vector<Alpha16Packet*>  udp;
   std::vector<Alpha16Channel*> hits;

   Alpha16Event(); // ctor
   ~Alpha16Event(); // dtor

   void Print() const;
};

struct Alpha16EVB
{
   int fEventCount = 0; // event counter

   bool     fHaveEventTs = false;
#if 0
   uint32_t fFirstEventTs[MAX_ALPHA16*NUM_CHAN_ALPHA16]; // udp timestamp of first event
   uint32_t fLastUdpEventTs[MAX_ALPHA16*NUM_CHAN_ALPHA16];  // udp timestamp of last seen event
#endif
   uint32_t fLastEventTs = 0;
   double   fLastEventTime = 0;
   int      fTsEpoch = 0;

#if 0
   int fConfNumChan;
   int fConfNumSamples;
   std::vector<int> fConfModMap;
#endif
   
   Alpha16EVB(); // ctor
   
   void Reset();
   //void Print() const;

   void Configure(int runno);

   Alpha16Event* NewEvent();
   void AddBank(Alpha16Event* e, Alpha16Packet* p, Alpha16Channel* c);
   void CheckEvent(Alpha16Event* e);

   //Alpha16Event* FindEvent(int imodule, uint32_t udpTs);
   //Alpha16Event* GetNextEvent();
   //static bool Match(const Alpha16Event* e, int imodule, uint32_t udpTs);

   Alpha16Map fMap;
};

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
