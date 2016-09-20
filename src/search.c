
#include <stdint.h>
#include <stdio.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"


// this searches for spoke patterns by testing for successful lacing
// patterns (against a rim modulo the pattern length).  duplicates
// are avoided using a simple bit sieve.

// in retrospect, for the limits used, a direct enumeration would
// likely have been fast enough.  but the search used here is
// significantly more efficient and could, with a different filter
// for duplicates, be used to search a larger space (eg all patterns
// for a given wheel size).


// two representations of the lacing are used: offsets and pattern.
// offsets is an array of integers, each integer is a spoke.
// pattern is a single (long) integer where every OFFSET_BITS is a spoke.
// they share the same bit-level representation, which is "signed"
// inside the OFFSET_BITS, where -0 is UNUSED_OFFSET (for extra
// confusion, but simplified code, the value -1 in the signed
// type used to represent offsets is used to initialise offsets
// before search - initial incrementation moves the value to 0).

// the "new" search, for any particular length of pattern,
// generates offsets that vary most rapidly in the highest index.
// since we want to generate patterns in lexical order we read
// this left (index 0) to right (index length-1).

// for groups A and B reflection (and negation) is necessary.
// this is about the "right hand end" (index length-1).

// when converting to a binary pattern, for simplicity we add
// and shift as we scan the array, so offset[0] is towards hsb
// while offset[length-1] is towards lsb.

// worked example:
// group A, after 1,-3 we backtrack to 2,0 (offsets).
// this is pattern 2,0A.  reflect to 2,0,-2
// binary pattern 010 000 110 (hex 86)


// these can be changed, but going much deeper requires too much
// memory for the sieve (length 10 requires 128MB and gives 10387
// patterns; length 12 requires 8GB and gives 72532 patterns).
#define OFFSET_BITS 3
#define MAX_LENGTH 6

#define PATTERN_FILE "patterns.txt"

// derived sizes
#define PATTERN_BITS (OFFSET_BITS * MAX_LENGTH)
#define UNUSED_OFFSET (1L << (OFFSET_BITS - 1))
#define MAX_OFFSET (UNUSED_OFFSET - 1)
#define OFFSET_SIGN UNUSED_OFFSET
#define OFFSET_VALUE MAX_OFFSET
#define OFFSET_MASK (OFFSET_SIGN | OFFSET_VALUE)
#define OFFSET_LIMIT (1L << OFFSET_BITS)
#define NEG(o) (o ? (o ^ OFFSET_SIGN) : o)
#define LENGTH_A ((MAX_LENGTH + 1) / 2)
#define LENGTH_B (MAX_LENGTH / 2)
#define LENGTH_C MAX_LENGTH
#define PATTERN_RIGHT_MASK (OFFSET_LIMIT - 1)
#define LEFT_ROTATION OFFSET_BITS

// all these are likely best as uint64 for fast access
#define SIEVE_T uint64_t
#define OFFSET_T int64_t   // minimum length OFFSET_BITS+1 must be signed
#define PATTERN_T uint64_t  // minimum length PATTERN_BITS
#define HOLES_T uint64_t  // minimum length MAX_LENGTH
// maybe the rim bit pattern type should be here too

// global state.  for a single-minded, math-intensive program it's
// pointless to pass these around as arguments
SIEVE_T *sieve = NULL;
lulog *dbg = NULL;
FILE *out = NULL;

// more derived sizes
#define SIEVE_WIDTH (8 * sizeof(*sieve))
#define SIEVE_LEN_BITS (1L << PATTERN_BITS)
#define SIEVE_LEN ((SIEVE_LEN_BITS + (SIEVE_WIDTH - 1)) / SIEVE_WIDTH)
#define SIEVE_LEN_BYTES (SIEVE_LEN_BITS / 8)

#define SIEVE_INDEX(n) (n / SIEVE_WIDTH)
#define SIEVE_SHIFT(n) (n - SIEVE_WIDTH * SIEVE_INDEX(n))
#define SET_SIEVE(n) (sieve[SIEVE_INDEX(n)] |= (1L << SIEVE_SHIFT(n)))
#define GET_SIEVE(n) (1 & (sieve[SIEVE_INDEX(n)] >> SIEVE_SHIFT(n)))

// the only mathematical insight is here.  that we can consider a pattern of length L
// as if it is laced to a tiny wheel with L holes (on one side) and so use modular
// arithmetic.
// length is repeated because % is remainder, not modulus, so we need positive values.
#define RIM_INDEX(offset, index, length) (1L << ((index + (offset & OFFSET_SIGN ? -1 : 1) * (offset & OFFSET_VALUE) + 2 * length) % length))


// equivalent to macros above - uncomment and lowercase call for debugging
//int get_sieve(int n) {
//    int index = SIEVE_INDEX(n);
//    int shift = SIEVE_SHIFT(n);
//    int s = GET_SIEVE(n);
//    ludebug(dbg, "Pattern %d -> sieve %d at %d/%d", n, s, index, shift);
//    return s;
//}
//
//void set_sieve(int n) {
//    int index = SIEVE_INDEX(n);
//    int shift = SIEVE_SHIFT(n);
//    SET_SIEVE(n);
//    ludebug(dbg, "Pattern %d -> sieve set at %d/%d", n, index, shift);
//}
//
//int rim_index(OFFSET_T offset, int index, int length) {
//    int sign = offset & OFFSET_SIGN;
//    int value = offset & OFFSET_VALUE;
//    ludebug(dbg, "Offset %d -> value %d sign %d", offset, value, sign);
//    int rim = 1L << ((index + (sign ? -1 : 1) * value + 2 * length) % length);
//    ludebug(dbg, "Offset %d at index %d -> rim %d", sign ? -value : value, index, rim);
//    return rim;
//}

// returns true if lacing is possible.  used to check padded patterns.
// could be used for a (much simpler) scan approach.
int check_lacing(PATTERN_T pattern, int length) {
    int rim = 0;
    for (int i = 0; i < length; ++i) {
        // unpacking from right first
        int addition = RIM_INDEX(pattern & PATTERN_RIGHT_MASK, length - i, length);
        pattern >>= OFFSET_BITS;
        if (rim & addition) {ludebug(dbg, "Bad lace %x / %x", rim, addition); return 0;}
        rim |= addition;
    }
    ludebug(dbg, "Laced ok, rim %x", rim);
    return 1;
}

void set_sieve_all_rotn(PATTERN_T pattern, int length) {
    PATTERN_T left_mask = ((1L << (length * OFFSET_BITS)) - 1) ^ PATTERN_RIGHT_MASK;
    int right_rotation = (length - 1) * OFFSET_BITS;
    for (int i = 0; i < length; ++i) {
        ludebug(dbg, "Setting %x, length %d", pattern, length);
        SET_SIEVE(pattern);
        pattern = ((pattern & left_mask) >> LEFT_ROTATION) | ((pattern & PATTERN_RIGHT_MASK) << right_rotation);
    }
}

void set_sieve_all(PATTERN_T pattern, int length) {
    PATTERN_T repeated = 0;
    int repeated_length = 0;
    while (repeated_length + length <= MAX_LENGTH) {
        repeated = (repeated << (length * OFFSET_BITS)) | pattern;
        repeated_length += length;
        set_sieve_all_rotn(repeated, repeated_length);
    }
}

void write_pattern(OFFSET_T *offsets, int length, char group, int padding, int full_length) {

    char buffer[3*MAX_LENGTH+3], *p;

    p = buffer;
    for (int i = 0; i < length; ++i) {
        if (i) p += sprintf(p, ",");
        int offset = offsets[i];
        if (offset > UNUSED_OFFSET) {
            p += sprintf(p, "-%d", NEG(offset));
        } else {
            p += sprintf(p, "%d", offset);
        }
    }
    *(p++) = group;
    if (padding) p += sprintf(p, "%d", padding);
    *p = '\0';

    luinfo(dbg, "Writing %s", buffer);
    fprintf(out, "%s %d\n", buffer, full_length);
}

int candidate_a(OFFSET_T *offsets, int length) {

    int count = 0;
    int half = (length + 1) / 2;
    PATTERN_T pattern = 0;

    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= offsets[i];}
    for (int i = 1; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= NEG(offsets[half - 1 - i]);}

    ludebug(dbg, "Candidate A length %d offsets %d %d %d -> %x", length, offsets[0], offsets[1], offsets[2], pattern);

    if (GET_SIEVE(pattern)) {
        ludebug(dbg, "Pattern %x already exists", pattern);
    } else if (length > 1 && !offsets[0]) {   // allow A0
        ludebug(dbg, "Skipping zero leading offset");
    } else if (offsets[0] & OFFSET_SIGN) {
        ludebug(dbg, "Skipping negative leading offset");
    } else if (offsets[half-1]) {
        ludebug(dbg, "Skipping non-radial central spoke");
    } else {
        int unbalanced = length == 1 && offsets[0];
        if (unbalanced) ludebug(dbg, "Unbalanced %d %d", length, pattern);
        for (int i = 0; i < MAX_LENGTH - length + 1; ++i) {
            PATTERN_T padded = pattern << (i * OFFSET_BITS);
            if (!unbalanced && !GET_SIEVE(padded) && check_lacing(padded, length + i)) {
                write_pattern(offsets, half, 'A', i, length + i);
                count++;
                set_sieve_all(padded, length + i);
            }
        }
    }

    return count;
}

void search_a() {

    luinfo(dbg, "Searching for A group patterns");

    OFFSET_T offsets[LENGTH_A] = {0};  // zeroed for nice display only
    int count = 0, length = 0;

    for (length = 1; length <= MAX_LENGTH; length += 2) {
        ludebug(dbg, "Looking for patterns of length %d", length);
        int half = (length + 1) / 2, middle = half - 1;
        for (int i = 0; i < half; ++i) offsets[i] = -1;
        int spoke = 0, rim = 0;
        while (spoke >= 0) {
            // at this point, spoke we are adjusting is not in rim
            int offset = offsets[spoke] + 1;
            if (offset == UNUSED_OFFSET) offset++;
            ludebug(dbg, "New offset for spoke %d is %d", spoke, offset);
            if (offset == OFFSET_LIMIT) {
                offsets[spoke] = -1;
                spoke--;
                // remove spoke we're backtracking to
                if (spoke >= 0) {
                    rim ^= RIM_INDEX(offsets[spoke], spoke, length);
                    if (spoke != middle) rim ^= RIM_INDEX(NEG(offsets[spoke]), -1 - spoke, length);
                }
                ludebug(dbg, "Out of options, so backtrack to spoke %d (rim %d)", spoke, rim);
            } else {
                offsets[spoke] = offset;
                int addition1 = RIM_INDEX(offset, spoke, length);
                int addition2 = RIM_INDEX(NEG(offset), -1 - spoke, length);
                if ((spoke == middle && !(rim & addition1)) || (spoke != middle && addition1 != addition2 && !((rim & addition1) | (rim & addition2)))) {
                    if (spoke == middle) {
                        ludebug(dbg, "Spoke(s) made rim complete");
                        count += candidate_a(offsets, length);
                        // we never added rim to spoke, so just continue
                    } else {
                        rim |= (addition1 | addition2);  // we're not at middle, so use both
                        spoke++;
                        ludebug(dbg, "Spoke(s) fits (rim %d), move to spoke %d", rim, spoke);
                    }
                }
            }
        }
    }

    luinfo(dbg, "Found %d A group patterns", count);
}

int candidate_b(OFFSET_T *offsets, int length) {

    int count = 0;
    int half = length / 2;
    PATTERN_T pattern = 0;

    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= offsets[i];}
    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= NEG(offsets[half - 1 - i]);}

    ludebug(dbg, "Candidate B length %d offsets %d %d %d -> %x", length, offsets[0], offsets[1], offsets[2], pattern);

    if (GET_SIEVE(pattern)) {
        ludebug(dbg, "Pattern %x already exists", pattern);
    } else if (!offsets[0]) {
        ludebug(dbg, "Skipping zero leading offset");
    } else if (offsets[0] & OFFSET_SIGN) {
        ludebug(dbg, "Skipping negative leading offset");
    } else {
        for (int i = 0; i < MAX_LENGTH - length + 1; ++i) {
            PATTERN_T padded = pattern << (i * OFFSET_BITS);
            if (!GET_SIEVE(padded) && check_lacing(padded, length + i)) {
                write_pattern(offsets, half, 'B', i, length + i);
                count++;
                set_sieve_all(padded, length + i);
            }
        }
    }

    return count;
}

void search_b() {

    luinfo(dbg, "Searching for B group patterns");

    OFFSET_T offsets[LENGTH_B] = {0};  // zeroed for nice display only
    int count = 0, length = 0;

    for (length = 2; length <= MAX_LENGTH; length += 2) {
        ludebug(dbg, "Looking for patterns of length %d", length);
        int half = length / 2;
        for (int i = 0; i < half; ++i) offsets[i] = -1;
        int spoke = 0, rim = 0;
        while (spoke >= 0) {
            // at this point, spoke we are adjusting is not in rim
            int offset = offsets[spoke] + 1;
            if (offset == UNUSED_OFFSET) offset++;
            ludebug(dbg, "New offset for spoke %d is %d", spoke, offset);
            if (offset == OFFSET_LIMIT) {
                offsets[spoke] = -1;
                spoke--;
                // remove spoke we're backtracking to
                if (spoke >= 0) {
                    rim ^= RIM_INDEX(offsets[spoke], spoke, length);
                    rim ^= RIM_INDEX(NEG(offsets[spoke]), -1 - spoke, length);
                }
                ludebug(dbg, "Out of options, so backtrack to spoke %d (rim %d)", spoke, rim);
            } else {
                offsets[spoke] = offset;
                int addition1 = RIM_INDEX(offset, spoke, length);
                int addition2 = RIM_INDEX(NEG(offset), -1 - spoke, length);
                if (!((rim & addition1) | (rim & addition2) | addition1 == addition2)) {
                    if (spoke == half-1) {
                        ludebug(dbg, "Spokes made rim complete");
                        count += candidate_b(offsets, length);
                        // we never added rim to spoke, so just continue
                    } else {
                        rim |= (addition1 | addition2);
                        spoke++;
                        ludebug(dbg, "Spokes fits (rim %d), move to spoke %d", rim, spoke);
                    }
                }
            }
        }
    }

    luinfo(dbg, "Found %d B group patterns", count);
}

int candidate_c(OFFSET_T *offsets, int length) {

    int count = 0;
    PATTERN_T pattern1 = 0, pattern2 = 0;

    for (int i = 0; i < length; ++i) {pattern1 <<= OFFSET_BITS; pattern1 |= offsets[i];}
    // negate and reverse to give second equivalent pattern
    // (does not apply to A/B because symmetric)
    for (int i = 0; i < length; ++i) {pattern2 <<= OFFSET_BITS; pattern2 |= NEG(offsets[length - 1 - i]);}

    ludebug(dbg, "Candidate C length %d offsets %d %d %d %d %d %d -> %x, %x", length,
            offsets[0], offsets[1], offsets[2], offsets[3], offsets[4], offsets[5], pattern1, pattern2);

    int pos = 0, neg = 0;
    for (int i = 0; i < length; ++i) {
        if (offsets[i]) {
            if (offsets[i] & OFFSET_SIGN) {
                neg = 1;
            } else {
                pos = 1;
            }
        }
    }

    if (!(pos & neg)) {
        ludebug(dbg, "All in one direction");
    } else if (GET_SIEVE(pattern1)) {
        ludebug(dbg, "Pattern %x already exists", pattern1);
    } else if (GET_SIEVE(pattern2)) {
        ludebug(dbg, "Pattern %x already exists", pattern2);
    } else if (!offsets[0]) {
        ludebug(dbg, "Skipping zero leading offset");
    } else if (offsets[0] & OFFSET_SIGN) {
        ludebug(dbg, "Skipping negative leading offset");
    } else {
        write_pattern(offsets, length, 'C', 0, length);
        set_sieve_all(pattern1, length);
        set_sieve_all(pattern2, length);
    }

    return count;
}

void search_c() {

    luinfo(dbg, "Searching for C group patterns");

    OFFSET_T offsets[LENGTH_C] = {0};  // zeroed for nice display only
    int count = 0, length = 0;

    for (length = 1; length <= MAX_LENGTH; ++length) {
        ludebug(dbg, "Looking for patterns of length %d", length);
        for (int i = 0; i < length; ++i) offsets[i] = -1;
        int spoke = 0, rim = 0;
        while (spoke >= 0) {
            // at this point, spoke we are adjusting is not in rim
            int offset = offsets[spoke] + 1;
            if (offset == UNUSED_OFFSET) offset++;
            ludebug(dbg, "New offset for spoke %d is %d", spoke, offset);
            if (offset == OFFSET_LIMIT) {
                offsets[spoke] = -1;
                spoke--;
                // remove spoke we're backtracking to
                if (spoke >= 0) rim ^= RIM_INDEX(offsets[spoke], spoke, length);
                ludebug(dbg, "Out of options, so backtrack to spoke %d (rim %d)", spoke, rim);
            } else {
                offsets[spoke] = offset;
                int addition = RIM_INDEX(offset, spoke, length);
                if (!(rim & addition)) {
                    if (spoke == length-1) {
                        ludebug(dbg, "Spoke made rim complete");
                        count += candidate_c(offsets, length);
                        // we never added rim to spoke, so just continue
                    } else {
                        rim |= addition;
                        spoke++;
                        ludebug(dbg, "Spoke fits (rim %d), move to spoke %d", rim, spoke);
                    }
                }
            }
        }
    }

    luinfo(dbg, "Found %d C group patterns", count);
}

void usage(const char *progname) {
    luinfo(dbg, "Search for spoke patterns (max offset %d, max length %d)", MAX_OFFSET, MAX_LENGTH);
    luinfo(dbg, "%s -h     display this message", progname);
    luinfo(dbg, "%s        run a search (output to %s)", progname, PATTERN_FILE);
}

// error handling is for lulib routines; don't bother elsewhere.
int main(int argc, char** argv) {

    LU_STATUS
    lulog_mkstdout(&dbg, lulog_level_debug);

    if (argc != 1) {
        usage(argv[0]);
    } else {

        luinfo(dbg, "Maximum spoke offset %d; Maximum pattern length %d", MAX_OFFSET, MAX_LENGTH);
        luinfo(dbg, "Sieve size %ldkB (%ld entries)", SIEVE_LEN_BYTES / 1024, SIEVE_LEN);
        LU_ALLOC(dbg, sieve, SIEVE_LEN)
        LU_ASSERT(!lufle_exists(dbg, PATTERN_FILE), LU_ERR_IO, dbg, "Output file %s already exists", PATTERN_FILE)
        lufle_open(dbg, PATTERN_FILE, "w", &out);

        search_a();
        search_b();
        search_c();

    }

LU_CLEANUP
    free(sieve);
    if (out) fclose(out);
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}
