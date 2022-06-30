// unpack_cb.cxx

#include "unpack_cb.h"

#include <stdio.h>
#include <iostream>

CbUnpack::CbUnpack(uint32_t num_inputs) // ctor
{
   fNumInputs = num_inputs;
   fChanLastTimestamp.resize(fNumInputs);
   fChanEpoch.resize(fNumInputs);
}

CbUnpack::~CbUnpack() // dtor
{

}

void CbUnpack::SaveScalers(CbScalers* scalers)
{
   //printf("save scalers!\n");
   
   if (fScalers.empty()) {
      fScalers.resize(fNumInputs+1);
   }
   
   fScalersEpoch.resize(fScalers.size());
   fScalersIncr.resize(fScalers.size());
   
   for (size_t i=0; i<fScalers.size(); i++) {
      if (fScalersPrev.size() > 0) {
         if (fScalers[i] < fScalersPrev[i]) {
            // 32-bit overflow
            fScalersEpoch[i] += 1;
         }
         fScalersIncr[i] = fScalers[i] - fScalersPrev[i]; // NB: this unsigned 32-bit arithmetic handles 32-bit overflow correctly. K.O.
      }
      //fScalersSum[i] = fScalers[i] + fScalersEpoch[i]*(double)0x100000000;
   }
   
   //double ext_clock_freq = 10e6; // 10 MHz timestamp clock
   double sys_clock_freq = 100e6; // 100 MHz timestamp clock
   
   double ext_clock_incr = fScalersIncr[fNumInputs-1];
   double sys_clock_incr = fScalersIncr[fNumInputs];
   
   double dt = sys_clock_incr/sys_clock_freq;
   
   double ext_clock = ext_clock_incr/dt;
   double sys_clock = sys_clock_incr/dt;

   if (fVerbose) {
      printf("clock incr %12.0f cnt %12.1f Hz, ext %12.0f cnt %12.1f Hz, dt %6.3f\n", sys_clock_incr, sys_clock, ext_clock_incr, ext_clock, dt);
   }

   CbLatchedScalers ls;
      
   for (size_t i=0; i<fScalers.size(); i++) {
      CbScaler s;

      s.raw   = fScalers[i];
      s.epoch = fScalersEpoch[i];
      s.sum = fScalers[i] + fScalersEpoch[i]*(double)0x100000000;
      if (dt > 0) {
         s.rate = fScalersIncr[i]/dt;
      } else {
         s.rate = 0;
      }

      ls.push_back(s);
   }
   
   if (scalers) {
      scalers->push_back(ls);
   }
   
   fScalersPrev = fScalers;
   fScalers.clear();
}

void CbUnpack::Unpack(const uint32_t* fifo_data, size_t nwords, CbHits* hits, CbScalers* scalers)
{
   if (fFailed)
      return;
  
   for (size_t i=0; i<nwords; i++) {
      uint32_t v = fifo_data[i];

      if (fVerbose) {
         printf("read %3zu: 0x%08x: ", i, v);
      }

      if (!fInScalersPacket && v == 0 && fKludge==1) {
         v = 0xFF000000;
         if (fVerbose) {
            printf("kludge1: corrected data word is 0x%08x: ", v);
         }
      }
      
      if (fInScalersPacket) {
         if (fVerbose) {
            printf(" scaler %d", (int)fScalers.size());
         }

         fScalers.push_back(v);

         if (fScalers.size() == fNumScalers) {
            fInScalersPacket = false;
            fNumScalers = 0;
            if (fVerbose) {
               printf(" (last)");
            }
            SaveScalers(scalers);
         }
      } else if ((v & 0xFF000000) == 0xFF000000) { // overflow marker
         int ts_top_bit = (v >> 23) & 1;
         uint32_t ts_counter = v & 0x7FFFFF;
         uint32_t ts_epoch = ts_counter/2;
         if (fKludge==1) {
            if (ts_top_bit==0 && (ts_counter&1))
               ts_epoch += 1;
         }
         if (fVerbose) {
            printf(" overflow 0x%06x, ts top bit %d, this_epoch %d", (int)ts_counter, ts_top_bit, (int)ts_epoch);
         }
         if (1) {
            if (ts_top_bit==0) {
               fCurrentEpoch = ts_epoch;
               fCurrentTsRangeMin = 0x800000;
               fCurrentTsRangeMax = 0xFFFFFF;
               for (size_t i=0; i<fChanEpoch.size(); i++) {
                  fChanEpoch[i] = fCurrentEpoch;
                  if (fWaitingForData) {
                     fChanLastTimestamp[i] = 0x700000;
                  } else if (fChanLastTimestamp[i] < 0x700000) {
                     fChanLastTimestamp[i] = 0x700000;
                  } else if (fChanLastTimestamp[i] > 0x900000) {
                     fChanLastTimestamp[i] = 0x700000;
                  }
               }
               if (fWaitingForData) {
                  if (fWaitForEpoch0 && ts_epoch != 0) {
                     if (fVerbose) {
                        printf(", waiting for epoch 0");
                     }
                  } else {
                     fWaitingForData = false;
                     if (fVerbose) {
                        printf(", start of data!");
                     }
                  }
               }
            } else {
               if (fWaitingForData) {
                  if (fVerbose) {
                     printf(", waiting for start of data");
                  }
               } else {
                  fCurrentTsRangeMin = 0x000000;
                  fCurrentTsRangeMax = 0x800000;
                  for (size_t i=0; i<fChanEpoch.size(); i++) {
                     if (fChanLastTimestamp[i] < 0xE00000)
                        fChanLastTimestamp[i] = 0xE00000;
                  }
               }
            }
         }
      } else if ((v & 0xFF000000) == 0xFE000000) { // start of scalers block
         fNumScalers = v & 0xFFFF;
         fInScalersPacket = true;
         if (fVerbose) {
            printf(" packet of %d scalers", fNumScalers);
         }
         
         if (fNumScalers != fNumInputs + 1) {
            fprintf(stderr, "CbUnpack::Unpack:: unexpected number of scalers %d vs expected %d!\n", fNumScalers, fNumInputs + 1);
            fFailed = true;
            return;
         }
         
      } else { // regular hit
         if (hits) {
            CbHit hit;
            //if (fCbEpochFromReset) {
            //   hit.epoch = fEpoch;
            //} else {
            //   hit.epoch = fEpochCounter;
            //}
            hit.timestamp = (v & 0x00FFFFFF);
            hit.channel   = (v & 0x7F000000)>>24;
            hit.flags = 0;
            if (v&1) {
               hit.flags |= CB_HIT_FLAG_TE;
            }

            if (fVerbose) {
               printf(" hit chan %2d, timestamp 0x%06x,",  hit.channel, hit.timestamp);
            }

            if (hit.channel < fNumInputs) {

               if (fWaitingForData) {
                  if (fVerbose) {
                     printf(" waiting for start of data");
                  }
               } else {
                  bool wrap = false;
                  if (hit.timestamp < fChanLastTimestamp[hit.channel]) {
                     fChanEpoch[hit.channel]++;
                     wrap = true;
                  }
                  hit.epoch = fChanEpoch[hit.channel];
                  hit.time = hit.timestamp/fCbTsFreq + 1.0/fCbTsFreq*hit.epoch*0x01000000;
                  
                  if (fVerbose) {
                     printf(" epoch %d, time %.6f sec", hit.epoch, hit.time);
                     if (wrap) {
                        printf(", wrap from timestamp 0x%06x", fChanLastTimestamp[hit.channel]);
                     }
                  }

                  hits->push_back(hit);

                  fChanLastTimestamp[hit.channel] = hit.timestamp;
               }
            } else {
               hit.epoch = 0;
               hit.time = 0;
               if (fVerbose) {
                  printf(" invalid channel number, epoch %d, time %.6f sec",  hit.epoch, hit.time);
               }
            }

            //fLastTimestamp = hit.timestamp;

            //if (hit.timestamp < fCurrentTsRangeMin || hit.timestamp > fCurrentTsRangeMax) {
            //   printf("BAD hit chan %2d, timestamp 0x%06x, epoch %d, time %.6f sec, 0x%08x, 0x%08x\n",  hit.channel, hit.timestamp, hit.epoch, hit.time, fCurrentTsRangeMin, fCurrentTsRangeMax);
            //}
         }
      }

      if (fVerbose) {
         printf("\n");
      }
   }
}

void PrintCbHits(const CbHits& hits)
{
   printf("PrintCbHits: %d hits\n", (int)hits.size());
   for (size_t i=0; i<hits.size(); i++) {
      printf("hit %3d: time %.6f sec, %d+%d, channel %d, %s\n", (int)i, hits[i].time, hits[i].timestamp, hits[i].epoch, hits[i].channel, (hits[i].flags&CB_HIT_FLAG_TE)?"te":"le");
   }
}

void PrintCbScalers(const CbScalers& scalers)
{
   printf("PrintCbScalers: %d sets of scalers\n", (int)scalers.size());
   for (size_t i=0; i<scalers.size(); i++) {
      printf("PrintCbScalers: set %d with %d scalers\n", (int)i, (int)scalers[i].size());
      for (size_t j=0; j<scalers[i].size(); j++) {
         printf("  chan %2d: raw 0x%08x (%d), epoch %.0f, sum %.0f, rate %g\n", (int)j, scalers[i][j].raw, scalers[i][j].raw, scalers[i][j].epoch, scalers[i][j].sum, scalers[i][j].rate);
      }
   }
}

/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
