//
// Unpacking FEAM data
// K.Olchanski
//

#ifndef Feam_H
#define Feam_H

#include <stdint.h>
#include <string>
#include <vector>
#include <utility>  // std::pair

//
// AFTER ASIC basic information:
//
// 511 SCA time bins
// 79 readout index entries
// 76 SCA channels
// 72 data channels
// 4 fixed pattern noise (fpn) channels, readout index 16, 29, 54, 67)
//
// Readout index, channel index, pin number:
//
// 1, 2, 3: reset1, reset2, reset3
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
   std::string fBank; // MIDAS bank, FEAM serial number
   int fPosition; // geographical position 0..63

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
   FeamModuleData(const FeamPacket* p, const char* bank, int position);
   ~FeamModuleData(); // dtor
   void AddData(const FeamPacket*p, int position, const char* ptr, int size);
   void Finalize();
   void Print() const;
};

#define MAX_FEAM           8 /* 64: 0..63 */
#define MAX_FEAM_SCA       4 /* 0,1,2,3 is A,B,C,D */
#define MAX_FEAM_READOUT  80 /* 1..79 */
#define MAX_FEAM_BINS    511 /* SCA time bins     */
#define MAX_FEAM_FPN       4 /* 4 FPN channels    */
#define MAX_FEAM_CHAN     72 /* 72 data channels  */
#define MAX_FEAM_PAD_COL   4 /* 4 TPC pad columns */
#define MAX_FEAM_PAD_ROWS (4*18) /* 4*18 TPC pad rows across 2 ASICs on one wing */

struct FeamAdcData
{
   int nsca;
   int nchan;
   int nbins;

   int adc[MAX_FEAM_SCA][MAX_FEAM_READOUT][MAX_FEAM_BINS];
};

struct FeamChannel
{
   std::string bank;
   int position; /* 0..63, div8 is vertical position, mod8 is azimuthal position */
   int sca; /* 0..3 */
   int sca_readout; /* 1..79 */
   int sca_chan; /* 1..72 */
   int tpc_col; /* 0..3 */
   int tpc_row; /* 0..72 */
   int first_bin; /* usually 0 */
   std::vector<int> adc_samples;
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
   std::vector<FeamChannel*> hits;

   FeamEvent(); // ctor
   ~FeamEvent(); // dtor
   void Print(int level=0) const;
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
