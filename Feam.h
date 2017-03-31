//
// Unpacking FEAM data
// K.Olchanski
//

#ifndef Feam_H
#define Feam_H

#include <stdint.h>
#include <vector>
#include <utility>  // std::pair

//
// AFTER ASIC basic information:
//
// 511 SCA time bins
// 76 SCA channels:
// 72 data channels
// 4 fixed pattern noise (fpn) channels: 13, 26, 51 and 64.
//
// Readout indices:
//
// 1, 2, 3: read
// 4..15 chan 1..12 pins 36..25
// 16 fpn1
// 17..28 chan 13..24 pins 24..13
// 29 fpn2
// 30..53 chan 25..48 pins 21..1, 120..109
// 54 fpn3
// 55..66 chan 49..60 pins 108..97
// 67 fpn4
// 68..79 chan 61..72 pins 96..85
//
// followed by:
// 8 bits of 0
// 9 bits of last read cell index (MSB first)
// 61 bits of 0
//

class FeamPacket
{
public:
   uint32_t cnt;
   uint16_t n;
   uint16_t x511;
   uint16_t buf_len;
   uint32_t ts_start;
   uint32_t ts_trig;
   int off;
   bool error;

public:
   FeamPacket(); // ctor
   ~FeamPacket(); // dtor
   void Unpack(const char* data, int size);
   void Print() const;
};

class FeamModuleData
{
public:
   int module;

   uint32_t cnt;
   uint32_t ts_start;
   uint32_t ts_trig;

   bool error;

   uint32_t next_n;

   int fSize;
   char* fPtr;

   uint32_t fTs;
   int      fTsEpoch;
   double   fTime;
   double   fTimeIncr;

public:
   FeamModuleData(const FeamPacket* p, int xmodule);
   ~FeamModuleData(); // dtor
   void AddData(const FeamPacket*p, int module, const char* ptr, int size);
   void Finalize();
   void Print() const;
};

#define MAX_FEAM   8
#define MAX_FEAM_SCA    4
#define MAX_FEAM_CHAN  80
#define MAX_FEAM_BINS 900
#define MAX_FEAM_PAD_COL    4
#define MAX_FEAM_PADS 72

struct FeamAdcData
{
   int nsca;
   int nchan;
   int nbins;

   int adc[MAX_FEAM_SCA][MAX_FEAM_CHAN][MAX_FEAM_BINS];
};

struct FeamEvent
{
   bool   complete; // event is complete
   bool   error;    // event has an error
   int    counter;  // event sequential counter
   double time;     // event time, sec
   double timeIncr; // time from previous event, sec

   std::vector<FeamModuleData*> modules;
   std::vector<FeamAdcData*> adcs;

   FeamEvent(); // ctor
   ~FeamEvent(); // dtor
   void Print() const;
};

extern void Unpack(FeamAdcData* a, FeamModuleData* m);
extern std::pair<int,int> getPad(short fAFTER, int index);

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
