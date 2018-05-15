//
// Unpacking FEAM data
// K.Olchanski
//

#ifndef Feam_H
#define Feam_H

#include <stdint.h>
#include <string>
#include <vector>
#include <deque>
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
   int fPosition = 0; // array index into the PWB event builder and event assembler
   int fModule = 0;   // PWB serial number NN from pwbNN
   int fColumn = 0;   // geographical position column 0..7 (as labeled on the detector)
   int fRing   = 0;   // geographical position ring 0..7 (bottom to top)
   int fDataFormat = 0; // FEAM/PWB data format, FEAM rev0 is 1, PWB rev1 is 2
   FeamPacket* fPacket = NULL; // event counter, timestamps, etc

   uint32_t cnt = 0;
   uint32_t ts_start = 0;
   uint32_t ts_trig  = 0;

   bool complete = false;
   bool error = false;

 public: // ADC data stream data

   static int fgMaxAlloc;
   int fAlloc = 0;
   int fSize = 0;
   char* fPtr = NULL;

 public: // FeamEVB and FeamAsm data

   uint32_t fTs = 0;
   int      fTsEpoch = 0;
   double   fTime = 0;
   double   fTimeIncr = 0;

   uint32_t fCntSeq = 0; // cnt - cnt_of_first_event

public:
   FeamModuleData(const FeamPacket* p, int position, int imodule, int icolumn, int iring, int format);
   ~FeamModuleData(); // dtor
   void AddData(const FeamPacket*p, const char* ptr, int size);
   void Finalize();
   void Print(int level=0) const;
};

class FeamModuleAsm
{
 public: // state
   int fState = 0;
   uint32_t fCnt = 0;
   int fNextN = 0;
   FeamModuleData* fCurrent = NULL;

 public: // timestamps
   uint32_t fTsFirstEvent = 0;
   uint32_t fTsLastEvent = 0;
   int      fTsEpoch = 0;
   double   fTimeFirstEvent = 0;
   double   fTimeLastEvent = 0;

   uint32_t fCntFirstEvent = 0;

 public: // context
   int fModule = 0;
   int fDataFormat = 0;

 public: // output buffer
   std::deque<FeamModuleData*> fBuffer;

 public: // counters
   int fCountIgnoredBeforeFirst = 0;
   int fCountFirst = 0;
   int fCountDone  = 0;
   int fCountLostFirst = 0;
   int fCountLostSync = 0;
   int fCountTruncated = 0;
   int fCountSkip = 0;
   int fCountWrongCnt = 0;

 public: // API
   ~FeamModuleAsm();
   void Print() const;
   void AddPacket(const FeamPacket* p, int position, int imodule, int icolumn, int iring, int format, const char* ptr, int size);
   void Finalize();

 public: // internal functions
   void StFirstPacket(const FeamPacket* p, int position, int imodule, int icolumn, int iring, int format, const char* ptr, int size);
   void StLastPacket();
   void AddData(const FeamPacket* p, const char* ptr, int size);
   void FlushIncomplete();
};

#define MAX_FEAM_SCA       4 /* 0,1,2,3 is A,B,C,D */
#define MAX_FEAM_READOUT  80 /* 1..79 */
#define MAX_FEAM_BINS    511 /* SCA time bins     */
#define MAX_FEAM_FPN       4 /* 4 FPN channels    */
#define MAX_FEAM_CHAN     72 /* 72 data channels  */
#define MAX_FEAM_PAD_COL   4 /* 4 TPC pad columns */
#define MAX_FEAM_PAD_ROWS (4*18) /* 4*18 TPC pad rows across 2 ASICs on one wing */

#define PWB_CHAN_INVALID (0) // invalid data from the device
#define PWB_CHAN_FPN1 (-1)   // SCA FPN channel
#define PWB_CHAN_FPN2 (-2)   // SCA FPN channel
#define PWB_CHAN_FPN3 (-3)   // SCA FPN channel
#define PWB_CHAN_FPN4 (-4)   // SCA FPN channel
#define PWB_CHAN_RESET1 (-5) // SCA Reset channel
#define PWB_CHAN_RESET2 (-6) // SCA Reset channel
#define PWB_CHAN_RESET3 (-7) // SCA Reset channel

#define PWB_MODULE_LAST 79 // serial number of last PWB module (last used pwbNN name)

class PwbPadMap
{
 public:
   PwbPadMap();
   bool CheckMap() const;

 public: // get the map singleton
   static const PwbPadMap* Map();
   static bool chan_is_pad(int chan);
   static bool chan_is_fpn(int chan);
   static bool chan_is_reset(int chan);

 public:
   int channel[MAX_FEAM_READOUT];
   int readout[MAX_FEAM_CHAN+1];
   int padcol[MAX_FEAM_SCA][MAX_FEAM_CHAN+1];
   int padrow[MAX_FEAM_SCA][MAX_FEAM_CHAN+1];
   int sca[MAX_FEAM_PAD_COL][MAX_FEAM_PAD_ROWS];
   int sca_chan[MAX_FEAM_PAD_COL][MAX_FEAM_PAD_ROWS];

 public:
   //std::pair<int,int> getPad(short sca, int readout_index) const;
};

class PwbModuleMapEntry
{
 public:
   int fModule; // PWB module pwbNN
   int fColumn; // TPC column, 0..7
   int fRing;   // TPC ring, 0..7
};

class PwbModuleMap
{
 public:
   std::vector<PwbModuleMapEntry*> fMap;
   int fNumModules = 0;
 public:
   PwbModuleMap(); // ctor
   ~PwbModuleMap(); // dtor
   void LoadFeamBanks(const std::vector<std::string> banks);
   const PwbModuleMapEntry* FindPwb(int module);
   void Print() const;
};

struct FeamAdcData
{
   int nsca;
   int nchan;
   int nbins;

   int adc[MAX_FEAM_SCA][MAX_FEAM_READOUT][MAX_FEAM_BINS];
};

struct FeamChannel
{
   int imodule;     /* pwb serial number NN from pwbNN, 0..80 or so */
   int pwb_column;  /* pwb column position, 0..7 */
   int pwb_ring;    /* pwb ring position, 0..7 */
   int sca;         /* 0..3 */
   int sca_readout; /* 1..79, includes 72 pad channels, 3 reset channels and 4 FPN channels */
   int sca_chan;    /* 1..72 */
   int pad_col;     /* local pad column position, 0..3 */
   int pad_row;     /* local pad row position, 0..71 */
   int first_bin;   /* usually 0 */
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

void PrintFeamChannels(const std::vector<FeamChannel*>& v);

extern void Unpack(FeamAdcData* a, FeamModuleData* m);

#endif

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
