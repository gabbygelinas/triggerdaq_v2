//
// ko_limits.h - common histogram limits and ranges, etc
//
// K.Olchanski
//

#ifndef KO_LIMITS_H
#define KO_LIMITS_H

// number of BSC bars (64 bottom, 64 top)
#define NUM_BSC 128

// number of TPC anode wires
#define NUM_AW 512

// TPC anode wire pulse or hit histogram limit
#define MAX_AW_AMP 64000

// General time limit, ns
#define MAX_TIME 8200

// TPC pad pulse or hit histogram limit
#define MAX_PAD_AMP 4100

// TPC field wire pulser time, ns

#define AW_PULSER_TIME_100 1800
#define AW_PULSER_TIME_625 1800
#define PAD_PULSER_TIME_625 7200

// number of TPC pad columns
#define NUM_PC (8*4)

// number of TPC pad rows
#define NUM_PR (18*4*8)

#endif
/* emacs
 * Local Variables:
 * tab-width: 8
 * c-basic-offset: 3
 * indent-tabs-mode: nil
 * End:
 */
