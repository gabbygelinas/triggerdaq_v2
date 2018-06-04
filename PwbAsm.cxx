//
// Unpacking PWB data
// K.Olchanski
//

#include "PwbAsm.h"

#include <stdio.h> // NULL, printf()
#include <math.h> // fabs()
#include <assert.h> // assert()

/* crc32c.c -- compute CRC-32C using the Intel crc32 instruction
 * Copyright (C) 2013 Mark Adler
 * Version 1.1  1 Aug 2013  Mark Adler
 */

#ifndef CRC32C_INCLUDE
#define CRC32C_INCLUDE

#include <stdint.h>
#include <string.h>
uint32_t crc32c_hw(uint32_t crc, const void *buf, size_t len);
uint32_t crc32c_sw(uint32_t crc, const void *buf, size_t len);
uint32_t crc32c(uint32_t crc, const void *buf, size_t len);

#endif

// CRC32C code taken from
// http://stackoverflow.com/questions/17645167/implementing-sse-4-2s-crc32c-in-software/17646775#17646775
// on 29 July 2015
//
// See also:
// https://tools.ietf.org/html/rfc3309
// https://tools.ietf.org/html/rfc3720#appendix-B.4
// and
// http://stackoverflow.com/questions/20963944/test-vectors-for-crc32c
// 

//#include "crc32c.h"

#ifdef __SSE__
#define HAVE_HWCRC32C 1
#endif

/* crc32c.c -- compute CRC-32C using the Intel crc32 instruction
 * Copyright (C) 2013 Mark Adler
 * Version 1.1  1 Aug 2013  Mark Adler
 */

/*
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the author be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Mark Adler
  madler@alumni.caltech.edu
 */

/* Use hardware CRC instruction on Intel SSE 4.2 processors.  This computes a
   CRC-32C, *not* the CRC-32 used by Ethernet and zip, gzip, etc.  A software
   version is provided as a fall-back, as well as for speed comparisons. */

/* Version history:
   1.0  10 Feb 2013  First version
   1.1   1 Aug 2013  Correct comments on why three crc instructions in parallel
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <pthread.h>

/* CRC-32C (iSCSI) polynomial in reversed bit order. */
#define POLY 0x82f63b78

/* Table for a quadword-at-a-time software crc. */
static pthread_once_t crc32c_once_sw = PTHREAD_ONCE_INIT;
static uint32_t crc32c_table[8][256];

/* Construct table for software CRC-32C calculation. */
static void crc32c_init_sw(void)
{
    uint32_t n, crc, k;

    for (n = 0; n < 256; n++) {
        crc = n;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc = crc & 1 ? (crc >> 1) ^ POLY : crc >> 1;
        crc32c_table[0][n] = crc;
    }
    for (n = 0; n < 256; n++) {
        crc = crc32c_table[0][n];
        for (k = 1; k < 8; k++) {
            crc = crc32c_table[0][crc & 0xff] ^ (crc >> 8);
            crc32c_table[k][n] = crc;
        }
    }
}

/* Table-driven software version as a fall-back.  This is about 15 times slower
   than using the hardware instructions.  This assumes little-endian integers,
   as is the case on Intel processors that the assembler code here is for. */
uint32_t crc32c_sw(uint32_t crci, const void *buf, size_t len)
{
   const unsigned char *next = (const unsigned char*)buf;
    uint64_t crc;

    pthread_once(&crc32c_once_sw, crc32c_init_sw);
    crc = crci ^ 0xffffffff;
    while (len && ((uintptr_t)next & 7) != 0) {
        crc = crc32c_table[0][(crc ^ *next++) & 0xff] ^ (crc >> 8);
        len--;
    }
    while (len >= 8) {
        crc ^= *(uint64_t *)next;
        crc = crc32c_table[7][crc & 0xff] ^
              crc32c_table[6][(crc >> 8) & 0xff] ^
              crc32c_table[5][(crc >> 16) & 0xff] ^
              crc32c_table[4][(crc >> 24) & 0xff] ^
              crc32c_table[3][(crc >> 32) & 0xff] ^
              crc32c_table[2][(crc >> 40) & 0xff] ^
              crc32c_table[1][(crc >> 48) & 0xff] ^
              crc32c_table[0][crc >> 56];
        next += 8;
        len -= 8;
    }
    while (len) {
        crc = crc32c_table[0][(crc ^ *next++) & 0xff] ^ (crc >> 8);
        len--;
    }
    return (uint32_t)crc ^ 0xffffffff;
}

/* Multiply a matrix times a vector over the Galois field of two elements,
   GF(2).  Each element is a bit in an unsigned integer.  mat must have at
   least as many entries as the power of two for most significant one bit in
   vec. */
static inline uint32_t gf2_matrix_times(uint32_t *mat, uint32_t vec)
{
    uint32_t sum;

    sum = 0;
    while (vec) {
        if (vec & 1)
            sum ^= *mat;
        vec >>= 1;
        mat++;
    }
    return sum;
}

/* Multiply a matrix by itself over GF(2).  Both mat and square must have 32
   rows. */
static inline void gf2_matrix_square(uint32_t *square, uint32_t *mat)
{
    int n;

    for (n = 0; n < 32; n++)
        square[n] = gf2_matrix_times(mat, mat[n]);
}

#ifdef HAVE_HWCRC32C

/* Construct an operator to apply len zeros to a crc.  len must be a power of
   two.  If len is not a power of two, then the result is the same as for the
   largest power of two less than len.  The result for len == 0 is the same as
   for len == 1.  A version of this routine could be easily written for any
   len, but that is not needed for this application. */
static void crc32c_zeros_op(uint32_t *even, size_t len)
{
    int n;
    uint32_t row;
    uint32_t odd[32];       /* odd-power-of-two zeros operator */

    /* put operator for one zero bit in odd */
    odd[0] = POLY;              /* CRC-32C polynomial */
    row = 1;
    for (n = 1; n < 32; n++) {
        odd[n] = row;
        row <<= 1;
    }

    /* put operator for two zero bits in even */
    gf2_matrix_square(even, odd);

    /* put operator for four zero bits in odd */
    gf2_matrix_square(odd, even);

    /* first square will put the operator for one zero byte (eight zero bits),
       in even -- next square puts operator for two zero bytes in odd, and so
       on, until len has been rotated down to zero */
    do {
        gf2_matrix_square(even, odd);
        len >>= 1;
        if (len == 0)
            return;
        gf2_matrix_square(odd, even);
        len >>= 1;
    } while (len);

    /* answer ended up in odd -- copy to even */
    for (n = 0; n < 32; n++)
        even[n] = odd[n];
}

/* Take a length and build four lookup tables for applying the zeros operator
   for that length, byte-by-byte on the operand. */
static void crc32c_zeros(uint32_t zeros[][256], size_t len)
{
    uint32_t n;
    uint32_t op[32];

    crc32c_zeros_op(op, len);
    for (n = 0; n < 256; n++) {
        zeros[0][n] = gf2_matrix_times(op, n);
        zeros[1][n] = gf2_matrix_times(op, n << 8);
        zeros[2][n] = gf2_matrix_times(op, n << 16);
        zeros[3][n] = gf2_matrix_times(op, n << 24);
    }
}

#endif

/* Apply the zeros operator table to crc. */
static inline uint32_t crc32c_shift(uint32_t zeros[][256], uint32_t crc)
{
    return zeros[0][crc & 0xff] ^ zeros[1][(crc >> 8) & 0xff] ^
           zeros[2][(crc >> 16) & 0xff] ^ zeros[3][crc >> 24];
}

/* Block sizes for three-way parallel crc computation.  LONG and SHORT must
   both be powers of two.  The associated string constants must be set
   accordingly, for use in constructing the assembler instructions. */
#define LONG 8192
#define LONGx1 "8192"
#define LONGx2 "16384"
#define SHORT 256
#define SHORTx1 "256"
#define SHORTx2 "512"

#ifdef HAVE_HWCRC32C

/* Tables for hardware crc that shift a crc by LONG and SHORT zeros. */
static pthread_once_t crc32c_once_hw = PTHREAD_ONCE_INIT;
static uint32_t crc32c_long[4][256];
static uint32_t crc32c_short[4][256];

/* Initialize tables for shifting crcs. */
static void crc32c_init_hw(void)
{
    crc32c_zeros(crc32c_long, LONG);
    crc32c_zeros(crc32c_short, SHORT);
}

/* Compute CRC-32C using the Intel hardware instruction. */
uint32_t crc32c_hw(uint32_t crc, const void *buf, size_t len)
{
   const unsigned char *next = (const unsigned char*)buf;
    const unsigned char *end;
    uint64_t crc0, crc1, crc2;      /* need to be 64 bits for crc32q */

    /* populate shift tables the first time through */
    pthread_once(&crc32c_once_hw, crc32c_init_hw);

    /* pre-process the crc */
    crc0 = crc ^ 0xffffffff;

    /* compute the crc for up to seven leading bytes to bring the data pointer
       to an eight-byte boundary */
    while (len && ((uintptr_t)next & 7) != 0) {
        __asm__("crc32b\t" "(%1), %0"
                : "=r"(crc0)
                : "r"(next), "0"(crc0));
        next++;
        len--;
    }

    /* compute the crc on sets of LONG*3 bytes, executing three independent crc
       instructions, each on LONG bytes -- this is optimized for the Nehalem,
       Westmere, Sandy Bridge, and Ivy Bridge architectures, which have a
       throughput of one crc per cycle, but a latency of three cycles */
    while (len >= LONG*3) {
        crc1 = 0;
        crc2 = 0;
        end = next + LONG;
        do {
            __asm__("crc32q\t" "(%3), %0\n\t"
                    "crc32q\t" LONGx1 "(%3), %1\n\t"
                    "crc32q\t" LONGx2 "(%3), %2"
                    : "=r"(crc0), "=r"(crc1), "=r"(crc2)
                    : "r"(next), "0"(crc0), "1"(crc1), "2"(crc2));
            next += 8;
        } while (next < end);
        crc0 = crc32c_shift(crc32c_long, crc0) ^ crc1;
        crc0 = crc32c_shift(crc32c_long, crc0) ^ crc2;
        next += LONG*2;
        len -= LONG*3;
    }

    /* do the same thing, but now on SHORT*3 blocks for the remaining data less
       than a LONG*3 block */
    while (len >= SHORT*3) {
        crc1 = 0;
        crc2 = 0;
        end = next + SHORT;
        do {
            __asm__("crc32q\t" "(%3), %0\n\t"
                    "crc32q\t" SHORTx1 "(%3), %1\n\t"
                    "crc32q\t" SHORTx2 "(%3), %2"
                    : "=r"(crc0), "=r"(crc1), "=r"(crc2)
                    : "r"(next), "0"(crc0), "1"(crc1), "2"(crc2));
            next += 8;
        } while (next < end);
        crc0 = crc32c_shift(crc32c_short, crc0) ^ crc1;
        crc0 = crc32c_shift(crc32c_short, crc0) ^ crc2;
        next += SHORT*2;
        len -= SHORT*3;
    }

    /* compute the crc on the remaining eight-byte units less than a SHORT*3
       block */
    end = next + (len - (len & 7));
    while (next < end) {
        __asm__("crc32q\t" "(%1), %0"
                : "=r"(crc0)
                : "r"(next), "0"(crc0));
        next += 8;
    }
    len &= 7;

    /* compute the crc for up to seven trailing bytes */
    while (len) {
        __asm__("crc32b\t" "(%1), %0"
                : "=r"(crc0)
                : "r"(next), "0"(crc0));
        next++;
        len--;
    }

    /* return a post-processed crc */
    return (uint32_t)crc0 ^ 0xffffffff;
}

/* Check for SSE 4.2.  SSE 4.2 was first supported in Nehalem processors
   introduced in November, 2008.  This does not check for the existence of the
   cpuid instruction itself, which was introduced on the 486SL in 1992, so this
   will fail on earlier x86 processors.  cpuid works on all Pentium and later
   processors. */
#define SSE42(have) \
    do { \
        uint32_t eax, ecx; \
        eax = 1; \
        __asm__("cpuid" \
                : "=c"(ecx) \
                : "a"(eax) \
                : "%ebx", "%edx"); \
        (have) = (ecx >> 20) & 1; \
    } while (0)

#endif // HAVE_HWCRC32C


/* Compute a CRC-32C.  If the crc32 instruction is available, use the hardware
   version.  Otherwise, use the software version. */
uint32_t crc32c(uint32_t crc, const void *buf, size_t len)
{
#ifdef HAVE_HWCRC32C
    int sse42;

    SSE42(sse42);
    return sse42 ? crc32c_hw(crc, buf, len) : crc32c_sw(crc, buf, len);
#else
#warning Hardware accelerated CRC32C is not available.
    return crc32c_sw(crc, buf, len);
#endif // HAVE_HWCRC32C
}

#ifdef TEST

#define SIZE (262144*3)
#define CHUNK SIZE

int main(int argc, char **argv)
{
    char *buf;
    ssize_t got;
    size_t off, n;
    uint32_t crc;

    (void)argv;
    crc = 0;
    buf = (char*)malloc(SIZE);
    if (buf == NULL) {
        fputs("out of memory", stderr);
        return 1;
    }
    while ((got = read(0, buf, SIZE)) > 0) {
        off = 0;
        do {
            n = (size_t)got - off;
            if (n > CHUNK)
                n = CHUNK;
            crc = argc > 1 ? crc32c_sw(crc, buf + off, n) :
                             crc32c(crc, buf + off, n);
            off += n;
        } while (off < (size_t)got);
    }
    free(buf);
    if (got == -1) {
        fputs("read error\n", stderr);
        return 1;
    }
    printf("%08x\n", crc);
    return 0;
}

#endif /* TEST */


PwbUdpPacket::PwbUdpPacket(const char* ptr, int size) // ctor
{
   const uint32_t* p32 = (const uint32_t*)ptr;
   int n32 = size/4;
   
   fError = false;
   fPacketSize = size;

   if (n32 < 6) {
      fError = true;
      printf("PwbUdpPacket: Error: invalid value of n32: %d\n", n32);
      return;
   }

   uint32_t crc = crc32c(0, ptr+0, 4*4);

   DEVICE_ID   = p32[0];
   PKT_SEQ     = p32[1];
   CHANNEL_SEQ = (p32[2] >>  0) & 0xFFFF;
   CHANNEL_ID  = (p32[2] >> 16) & 0xFF;
   FLAGS       = (p32[2] >> 24) & 0xFF;
   CHUNK_ID    = (p32[3] >>  0) & 0xFFFF;
   CHUNK_LEN   = (p32[3] >> 16) & 0xFFFF;
   HEADER_CRC  = p32[4];
   start_of_payload = 5*4;
   end_of_payload = start_of_payload + CHUNK_LEN;

   if (HEADER_CRC != ~crc) {
      printf("PwbUdpPacket: Error: header CRC mismatch: HEADER_CRC 0x%08x, computed 0x%08x\n", HEADER_CRC, ~crc);
      fError = true;
      return;
   }

   if (end_of_payload + 4 > fPacketSize) {
      printf("PwbUdpPacket: Error: invalid byte counter, end_of_payload %d, packet size %d\n", end_of_payload, fPacketSize);
      fError = true;
      return;
   }
   
   payload_crc = p32[end_of_payload/4];

   uint32_t payload_crc_ours = crc32c(0, ptr+start_of_payload, CHUNK_LEN);

   if (payload_crc != ~payload_crc_ours) {
      printf("PwbUdpPacket: Error: payload CRC mismatch: CRC 0x%08x, computed 0x%08x\n", payload_crc, ~payload_crc_ours);
      fError = true;
      return;
   }
}

void PwbUdpPacket::Print() const
{
   printf("PwbUdpPacket: DEVICE_ID 0x%08x, PKT_SEQ 0x%08x, CHAN SEQ 0x%02x, ID 0x%02x, FLAGS 0x%02x, CHUNK ID 0x%04x, LEN 0x%04x, CRC 0x%08x, bank bytes %d, end of payload %d, CRC 0x%08x\n",
          DEVICE_ID,
          PKT_SEQ,
          ((int)CHANNEL_SEQ)&0xFFFF,
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
   int nw32 = size/4;
   
   if (nw32 < 16) {
      printf("PwbEventHeader::ctor: Error: event header size %d words is too small\n", nw32);
      fError = true;
      return;
   }

   fError = false;
   
   FormatRevision  = (p32[5]>> 0) & 0xFF;
   ScaId           = (p32[5]>> 8) & 0xFF;
   CompressionType = (p32[5]>>16) & 0xFF;
   TriggerSource   = (p32[5]>>24) & 0xFF;

   if (FormatRevision == 0) {
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
   } else if (FormatRevision == 1) {
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
      
      ScaChannelsSent3 = (p32[13]>> 0) & 0xFFFF;
      ScaChannelsThreshold1 = (p32[13]>>16) & 0xFFFF;
      
      ScaChannelsThreshold2 = p32[14];
      
      ScaChannelsThreshold3 = p32[15];
      Reserved2             = 0;
      
      start_of_data = 16*4;
   } else {
      printf("PwbEventHeader::ctor: Error: invalid FormatRevision %d, expected 0 or 1\n", FormatRevision);
      fError = true;
   }
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

#define PWB_CA_ST_ERROR -1
#define PWB_CA_ST_INIT 0
#define PWB_CA_ST_DATA 1
#define PWB_CA_ST_LAST 2

PwbChannelAsm::PwbChannelAsm(int module, int column, int ring, int sca)
{
   fModule = module;
   fColumn = column;
   fRing = ring;
   fSca = sca;
   Reset();
}

PwbChannelAsm::~PwbChannelAsm()
{
   if (fCountErrors) {
      printf("PwbChannelAsm: module %d sca %d: %d errors\n", fModule, fSca, fCountErrors);
   }
}

void PwbChannelAsm::Reset()
{
   fLast_CHANNEL_SEQ = 0;
   fState = 0;
   fSaveChannel = 0;
   fSaveSamples = 0;
   fSaveNw = 0;
   fSavePos = 0;
   if (fCurrent) {
      delete fCurrent;
      fCurrent = NULL;
   }
   for (unsigned i=0; i<fOutput.size(); i++) {
      if (fOutput[i]) {
         delete fOutput[i];
         fOutput[i] = NULL;
      }
   }
}

void PwbChannelAsm::AddSamples(int channel, const uint16_t* samples, int count)
{
   const int16_t* signed_samples = (const int16_t*)samples;
   int ri = -1;
   if (fFormatRevision == 0) {
      ri = channel + 1;
   } else if (fFormatRevision == 1) {
      ri = channel;
   } else {
      printf("PwbChannelAsm::AddSamples: Error: module %d sca %d state %d: invalid FormatRevision %d\n", fModule, fSca, fState, fFormatRevision);
      fCountErrors++;
      fState = PWB_CA_ST_ERROR;
      fError = true;
      return;
   }

   if (fTrace) {
      printf("pwb module %d, sca %d, channel %d, ri %d, add %d samples\n", fModule, fSca, channel, ri, count);
   }

   if (fCurrent) {
      if (ri != fCurrent->sca_readout) {
         fOutput.push_back(fCurrent);
         fCurrent = NULL;
      }
   }
   
   if (!fCurrent) {
      fCurrent = new FeamChannel;
      fCurrent->imodule = fModule;
      fCurrent->pwb_column = fColumn;
      fCurrent->pwb_ring = fRing;
      fCurrent->sca = fSca;
      fCurrent->sca_readout = ri;
      fCurrent->sca_chan = 0; // filled by Build()
      fCurrent->pad_col = 0; // filled by Build()
      fCurrent->pad_row = 0; // filled by Build()
      fCurrent->first_bin = 0;
   }

   for (int i=0; i<count; i++) {
      fCurrent->adc_samples.push_back(signed_samples[i]);
   }

   if (0) {
      if (fSca == 0 && ri==1) {
         printf("pwb module %d, sca %d, channel %d, ri %d, add %d samples: \n", fModule, fSca, channel, ri, count);
         for (int i=0; i<count; i++) {
            printf(" %4d", signed_samples[i]);
         }
         printf("\n");
      }
   }
}

void PwbChannelAsm::CopyData(const uint16_t* s, const uint16_t* e)
{
   const uint16_t* p = s;

   while (1) {
      int r = e-p;
      if (r < 2) {
         printf("PwbChannelAsm::CopyData: Error: need to en-buffer a partial header, r: %d!\n", r);
         fCountErrors++;
         break;
      }
      if (p[0] == 0xCCCC && p[1] == 0xCCCC && r == 2) {
         printf("PwbChannelAsm::CopyData: module %d sca %d: ignoring unexpected 0xCCCC words at the end of a packet\n", fModule, fSca);
         //fCountErrors++;
         break;
      }
      if (p[0] == 0xCCCC && p[1] == 0xCCCC) {
         printf("PwbChannelAsm::CopyData: Error: module %d, unexpected 0xCCCC words, r %d!\n", fModule, r);
         fCountErrors++;
         break;
      }
      int channel = p[0];
      int samples = p[1];

      if (channel < 0 || channel >= 80) {
         printf("PwbChannelAsm::CopyData: Error: module %d, invalid channel %d\n", fModule, channel);
         fCountErrors++;
      }

      if (samples != 511) {
         printf("PwbChannelAsm::CopyData: Error: invalid samples %d\n", samples);
         fCountErrors++;
      }

      int nw = samples;
      if (samples&1)
         nw+=1;
      if (nw <= 0) {
         printf("PwbChannelAsm::CopyData: Error: invalid word counter nw: %d\n", nw);
         fCountErrors++;
         break;
      }
      if (fTrace) {
         printf("adc samples ptr %d, end %d, r %d, channel %d, samples %d, nw %d\n", (int)(p-s), (int)(e-s), r, channel, samples, nw);
      }
      p += 2;
      r = e-p;
      int h = nw;
      int s = samples;
      bool truncated = false;
      if (nw > r) {
         h = r;
         s = r;
         if (fTrace) {
            printf("split data nw %d, with %d here, %d to follow\n", nw, h, nw-h);
         }
         fSaveChannel = channel;
         fSaveSamples = samples;
         fSaveNw = nw;
         fSavePos = h;
         truncated = true;
      }
      AddSamples(channel, p, s);
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

void PwbChannelAsm::BeginData(const char* ptr, int size, int start_of_data, int end_of_data, uint32_t ts)
{
   fTs = ts;
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
      int s = fSaveSamples - fSavePos;
      bool truncated = false;
      
      if (h > r) {
         h = r;
         s = r;
         truncated = true;
      }

      if (fTrace) {
         printf("AddData: save channel %d, samples %d, nw %d, pos %d, remaining %d, have %d, truncated %d\n", fSaveChannel, fSaveSamples, fSaveNw, fSavePos, r, h, truncated);
      }
      
      AddSamples(fSaveChannel, p, s);
      
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
      printf("PwbChannelAsm::AddData: Error!\n");
      assert(!"this cannot happen");
   }
}

void PwbChannelAsm::EndData()
{
   if (fSaveNw > 0) {
      printf("PwbChannelAsm::EndData: Error: missing some data at the end\n");
   } else {
      if (fTrace) {
         printf("PwbChannelAsm::EndData: ok!\n");
      }
   }

   fSaveChannel = 0;
   fSaveSamples = 0;
   fSaveNw = 0;
   fSavePos = 0;

   if (fCurrent) {
      fOutput.push_back(fCurrent);
      fCurrent = NULL;
   }

   if (fTrace) {
      PrintFeamChannels(fOutput);
   }
}

void PwbChannelAsm::BuildEvent(FeamEvent* e)
{
   if (fError) {
      e->error = true;
   }

   for (unsigned i=0; i<fOutput.size(); i++) {
      if (fOutput[i]) {
         FeamChannel* c = fOutput[i];
         fOutput[i] = NULL;
         
         // FIXME: bad data from PWB!!!
         if (c->sca_readout >= MAX_FEAM_READOUT) {
            printf("PwbChannelAsm::BuildEvent: Error: skipping invalid channel, sca_readout: %d\n", c->sca_readout);
            delete c;
            continue;
         }

         c->sca_chan = PwbPadMap::Map()->channel[c->sca_readout];
         bool scachan_is_pad = (c->sca_chan > 0);
         //bool scachan_is_fpn = (c->sca_chan >= -4) && (c->sca_chan <= -1);
         
         if (scachan_is_pad) {
            c->pad_col = PwbPadMap::Map()->padcol[c->sca][c->sca_chan];
            c->pad_row = PwbPadMap::Map()->padrow[c->sca][c->sca_chan];
         }

         e->hits.push_back(c);
      }
   }
   fOutput.clear();
}

void PwbChannelAsm::AddPacket(PwbUdpPacket* udp, const char* ptr, int size)
{
   if (fLast_CHANNEL_SEQ == 0) {
      fLast_CHANNEL_SEQ = udp->CHANNEL_SEQ;
   } else if (udp->CHANNEL_SEQ != ((fLast_CHANNEL_SEQ + 1)&0xFFFF)) {
      if (fState == PWB_CA_ST_ERROR) {
         printf("PwbChannelAsm::AddPacket: module %d sca %d state %d: misordered or lost UDP packet: CHANNEL_SEQ jump 0x%04x to 0x%04x\n", fModule, fSca, fState, fLast_CHANNEL_SEQ, udp->CHANNEL_SEQ);
         fError = true;
      } else {
         printf("PwbChannelAsm::AddPacket: Error: module %d sca %d state %d: misordered or lost UDP packet: CHANNEL_SEQ jump 0x%04x to 0x%04x\n", fModule, fSca, fState, fLast_CHANNEL_SEQ, udp->CHANNEL_SEQ);
         fLast_CHANNEL_SEQ = udp->CHANNEL_SEQ;
         fCountErrors++;
         fState = PWB_CA_ST_ERROR;
         fError = true;
      }
   } else {
      fLast_CHANNEL_SEQ++;
   }

   if (fState == PWB_CA_ST_INIT || fState == PWB_CA_ST_LAST || fState == PWB_CA_ST_ERROR) {
      if (udp->CHUNK_ID == 0) {
         PwbEventHeader* eh = new PwbEventHeader(ptr, size);
         uint32_t ts = eh->TriggerTimestamp1;
         fFormatRevision = eh->FormatRevision;
         if (fTrace) {
            eh->Print();
         }
         if (eh->fError) {
            printf("PwbChannelAsm::AddPacket: Error: module %d sca %d state %d: error in event header\n", fModule, fSca, fState);
            fCountErrors++;
            fState = PWB_CA_ST_ERROR;
            fError = true;
         } else {
            fState = PWB_CA_ST_DATA;
            fError = false;
            BeginData(ptr, size, eh->start_of_data, udp->end_of_payload, ts);
         }
         delete eh;
      } else {
         printf("PwbChannelAsm::AddPacket: module %d sca %d state %d: Ignoring UDP packet with CHUNK_ID 0x%02x while waiting for an event header\n", fModule, fSca, fState, udp->CHUNK_ID);
      }
   } else if (fState == PWB_CA_ST_DATA) {
      AddData(ptr, size, udp->start_of_payload, udp->end_of_payload);
      if (udp->FLAGS & 1) {
         fState = PWB_CA_ST_LAST;
         EndData();
      }
   } else {
      printf("PwbChannelAsm::AddPacket: Error: module %d sca %d state %d: invalid state\n", fModule, fSca, fState);
      fCountErrors++;
   }
}

bool PwbChannelAsm::CheckComplete() const
{
   if (fState == PWB_CA_ST_INIT)
      return true;
   if (fState == PWB_CA_ST_LAST)
      return true;
   return false;
}

PwbModuleAsm::PwbModuleAsm(int module, int column, int ring)
{
   fModule = module;
   fColumn = column;
   fRing = ring;
   Reset();
}

PwbModuleAsm::~PwbModuleAsm()
{
   for (unsigned i=0; i<fChannels.size(); i++) {
      if (fChannels[i]) {
         delete fChannels[i];
         fChannels[i] = NULL;
      }
   }

   if (fCountErrors) {
      printf("PwbModuleAsm: module %d: %d errors, lost %d/%d/%d/%d\n", fModule, fCountErrors, fCountLost1, fCountLost2, fCountLost3, fCountLostN);
   }
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

   if (fTrace) {
      udp->Print();
   }

   if (fLast_PKT_SEQ == 0) {
      fLast_PKT_SEQ = udp->PKT_SEQ;
   } else if (udp->PKT_SEQ == fLast_PKT_SEQ + 1) {
      fLast_PKT_SEQ++;
   } else if (udp->PKT_SEQ == fLast_PKT_SEQ + 2) {
      printf("PwbModuleAsm::AddPacket: Error: module %d lost one UDP packet: PKT_SEQ jump 0x%08x to 0x%08x\n", fModule, fLast_PKT_SEQ, udp->PKT_SEQ);
      fCountErrors++;
      fCountLost1++;
      fLast_PKT_SEQ = udp->PKT_SEQ;
   } else if (udp->PKT_SEQ == fLast_PKT_SEQ + 3) {
      printf("PwbModuleAsm::AddPacket: Error: module %d lost two UDP packets: PKT_SEQ jump 0x%08x to 0x%08x\n", fModule, fLast_PKT_SEQ, udp->PKT_SEQ);
      fCountErrors++;
      fCountLost2++;
      fLast_PKT_SEQ = udp->PKT_SEQ;
   } else if (udp->PKT_SEQ == fLast_PKT_SEQ + 4) {
      printf("PwbModuleAsm::AddPacket: Error: module %d lost three UDP packets: PKT_SEQ jump 0x%08x to 0x%08x\n", fModule, fLast_PKT_SEQ, udp->PKT_SEQ);
      fCountErrors++;
      fCountLost3++;
      fLast_PKT_SEQ = udp->PKT_SEQ;
   } else {
      int skip = udp->PKT_SEQ - fLast_PKT_SEQ - 1;
      printf("PwbModuleAsm::AddPacket: Error: module %d misordered or lost UDP packets: PKT_SEQ jump 0x%08x to 0x%08x, skipped %d\n", fModule, fLast_PKT_SEQ, udp->PKT_SEQ, skip);
      fCountErrors++;
      fCountLostN++;
      fLast_PKT_SEQ = udp->PKT_SEQ;
   }

   int s = udp->CHANNEL_ID;

   if (s < 0 || s > 4) {
      printf("PwbModuleAsm::AddPacket: Error: invalid CHANNEL_ID 0x%08x\n", udp->CHANNEL_ID);
      fCountErrors++;
   } else {
      while (s >= fChannels.size()) {
         fChannels.push_back(NULL);
      }

      if (!fChannels[s]) {
         fChannels[s] = new PwbChannelAsm(fModule, fColumn, fRing, s);
      }

      fChannels[s]->AddPacket(udp, ptr, size);
   }
   
   delete udp;
}

bool PwbModuleAsm::CheckComplete() const
{
   for (unsigned i=0; i<fChannels.size(); i++) {
      if (fChannels[i]) {
         if (!fChannels[i]->CheckComplete())
            return false;
      }
   }
   return true;
}

void PwbModuleAsm::BuildEvent(FeamEvent* e)
{
   fTs = 0;
   for (unsigned i=0; i<fChannels.size(); i++) {
      if (fChannels[i]) {
         fChannels[i]->BuildEvent(e);
         if (fTs == 0) {
            fTs = fChannels[i]->fTs;
         } else {
            if (fChannels[i]->fTs != fTs) {
               printf("PwbModuleAsm::BuildEvent: Error: channel %d timestamp mismatch 0x%08x should be 0x%08x\n", i, fChannels[i]->fTs, fTs);
               fCountErrors++;
               e->error = true;
            }
         }
      }
   }

   double ts_freq = 125000000.0; // 125 MHz

   if (e->counter == 1) {
      fTsFirstEvent = fTs;
      fTsLastEvent = 0;
      fTsEpoch = 0;
      fTimeFirstEvent = fTs/ts_freq;
   }
   
   if (fTs < fTsLastEvent)
      fTsEpoch++;
   
   fTime = fTs/ts_freq - fTimeFirstEvent + fTsEpoch*2.0*0x80000000/ts_freq;
   fTimeIncr = fTime - fTimeLastEvent;
   
   fTsLastEvent = fTs;
   fTimeLastEvent = fTime;
}

PwbAsm::PwbAsm() // ctor
{
   // empty
}

PwbAsm::~PwbAsm() // dtor
{
   int errors = 0;
   int lost1 = 0;
   int lost2 = 0;
   int lost3 = 0;
   int lostN = 0;

   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i]) {
         errors += fModules[i]->fCountErrors;
         lost1  += fModules[i]->fCountLost1;
         lost2  += fModules[i]->fCountLost2;
         lost3  += fModules[i]->fCountLost3;
         lostN  += fModules[i]->fCountLostN;
         delete fModules[i];
         fModules[i] = NULL;
      }
   }

   if (errors > 0) {
      printf("PwbAsm: %d errors, lost %d/%d/%d/%d\n", errors, lost1, lost2, lost3, lostN);
   }
}

void PwbAsm::Reset()
{
   fCounter = 0;
   fLastTime = 0;
   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i]) {
         fModules[i]->Reset();
      }
   }
}

void PwbAsm::AddPacket(int module, int column, int ring, const char* ptr, int size)
{
   while (module >= fModules.size()) {
      fModules.push_back(NULL);
   }

   if (!fModules[module]) {
      fModules[module] = new PwbModuleAsm(module, column, ring);
   }

   fModules[module]->AddPacket(ptr, size);
}

bool PwbAsm::CheckComplete() const
{
   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i]) {
         if (!fModules[i]->CheckComplete())
            return false;
      }
   }
   return true;
}

void PwbAsm::BuildEvent(FeamEvent* e)
{
   e->complete = true;
   e->error = false;
   e->counter = ++fCounter;
   e->time = 0;
   e->timeIncr = 0;

   bool first_ts = true;

   for (unsigned i=0; i<fModules.size(); i++) {
      if (fModules[i]) {
         fModules[i]->BuildEvent(e);
         if (first_ts) {
            first_ts = false;
            e->time = fModules[i]->fTime;
         } else {
            if (fabs(e->time - fModules[i]->fTime) > 1e9) {
               printf("PwbModuleAsm::BuildEvent: Error: module %d event time mismatch %f should be %f sec, diff %f ns\n", i, fModules[i]->fTime, e->time, (fModules[i]->fTime - e->time)*1e9);
               e->error = true;
            }
         }
      }
   }

   e->timeIncr = e->time - fLastTime;
   fLastTime = e->time;
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */


