
#include <stdint.h>
#include <stdio.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"


// this searches for spoke patterns by testing for successful lacing
// patterns.  duplicates are avoided using a simple bit sieve.

// in retrospect, for the limits used, a direct enumeration would
// likely have been fast enough.  but the search used here is
// significantly more efficient and could, with a different filter
// for duplicates, be used to search a larger space (eg all patterns
// for a given wheel size).


// two representations of the lacing are used: offsets and pattern.
// offsets is an array of integers, each integer is a spoke.
// pattern is a single (long) integer where every OFFSET_BITS is a spoke.
// they share the same bit-level representation, which is signed
// inside the OFFSET_BITS, where -0 is UNUSED_OFFSET.
// the 0 index offset is the "least significant" spoke in a search.
// for A and B patterns, that is the "centermost" spoke(s).
// for C patterns it is the "rightmost" spoke.
// in both cases, the textual name is generated by writing the
// values from right (0 index) to left (length(offsets)-1 index).
// the pattern, when visualised as a binary sequence, with lsb to
// the right, has the same order as the name, but is fully expanded
// (group A and B are reflected).

// worked example 1:
// spoke pattern: 2, 0, -2, 0 (going clockwise round wheel)
// this is extended crow's foot.  group A with padding. 2,0A1
// encoded as offset values (if OFFSET_BITS = 3): 2, 0, 4, 0
// encoded as offset array [0, 2]
// encoded as binary pattern ...010000100 plus 1 padding
// or as ...010000100000 with padding included

// worked example 2:
// spoke pattern 1, -3, 2 (random mess - may not actually lace)
// name is 1,-3,2C
// encoded as offset values: 1,7,2
// encoded as offset array [2, 7, 1]
// encoded as binary pattern: ...001111010

// as can be seen, the least significant bits of the pattern correspond
// to the 0 index of the offset array for C, and to the reflected spokes
// (which are negated) for A and B.  so "reading" (l-r) the pattern gives
// the same order as moving clockwise round the wheel.


// these can be changed, but going much deeper requires too much
// memory for the sieve.
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


//int get_sieve(int n) {
//    int index = SIEVE_INDEX(n);
//    int shift = SIEVE_SHIFT(n);
//    int s = GET_SIEVE(n);
//    ludebug(dbg, "Pattern %d -> sieve %d at %d/%d", n, s, index, shift);
//    return s;
//}
//
void set_sieve(int n) {
    int index = SIEVE_INDEX(n);
    int shift = SIEVE_SHIFT(n);
    SET_SIEVE(n);
    ludebug(dbg, "Pattern %d -> sieve set at %d/%d", n, index, shift);
}
//
//int rim_index(OFFSET_T offset, int index, int length) {
//    int sign = offset & OFFSET_SIGN;
//    int value = offset & OFFSET_VALUE;
//    ludebug(dbg, "Offset %d -> value %d sign %d", offset, value, sign);
//    int rim = 1L << ((index + (sign ? -1 : 1) * value + 2 * length) % length);
//    ludebug(dbg, "Offset %d at index %d -> rim %d", sign ? -value : value, index, rim);
//    return rim;
//}

// this shouldn't be necessary as all candidates are constructed, but
// it does allow enumeration of possible values from sieve gaps.
void flag_unused() {

    luinfo(dbg, "Setting sieve entries with unused offsets");

    OFFSET_T offsets[MAX_LENGTH] = {0};
    int i = 0, count = 0;
    OFFSET_T offset = 0;

    // run through all possible patterns
    while (i < MAX_LENGTH) {  // overflow

        // increment pattern
        for (i = 0, offset = 0; i < MAX_LENGTH && !offset; ++i) {
            offset = offsets[i] + 1;  // increment current digit
            if (offset == OFFSET_LIMIT) offset = 0;  // carry
            offsets[i] = offset;
        }

        // check for unused
        int unused = 0;
        for (int j = 0; j < MAX_LENGTH && !unused; ++j) unused = offsets[j] == UNUSED_OFFSET;

        if (unused) {
            // build and set pattern
            PATTERN_T pattern = 0;
            for (int k = 0; k < MAX_LENGTH; ++k) {
                pattern <<= OFFSET_BITS;
                pattern |= offsets[MAX_LENGTH - k - 1];
            }
            count++;
            SET_SIEVE(pattern);
        }

    }

    luinfo(dbg, "Set %d (/%d = %2.0f%%)  entries", count, SIEVE_LEN_BITS, count * 100.0 / SIEVE_LEN_BITS);
}

// returns true if lacing is possible.
// could be used for a (much simpler) scan approach.
int check_lacing(PATTERN_T pattern, int length) {
    int rim = 0;
    for (int i = 0; i < length; ++i) {
        // -i because we're unpacking pattern "backwards"
        int addition = RIM_INDEX(pattern & PATTERN_RIGHT_MASK, -i, length);
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
        set_sieve(pattern);
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
        if (offsets[length - i - 1] > UNUSED_OFFSET) {
            p += sprintf(p, "-%d", NEG(offsets[length - i - 1]));
        } else {
            p += sprintf(p, "%d", offsets[length - i - 1]);
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

    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= offsets[half - i - 1];}
    for (int i = 1; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= NEG(offsets[i]);}

    ludebug(dbg, "Candidate A length %d offsets %d %d %d -> %x", length, offsets[0], offsets[1], offsets[2], pattern);

    if (GET_SIEVE(pattern)) {
        ludebug(dbg, "Pattern %x already exists", pattern);
    } else if (half > 1 && !offsets[half-1]) {   // allow A0
        ludebug(dbg, "Skipping zero leading offset");
    } else if (offsets[half-1] > UNUSED_OFFSET) {
        ludebug(dbg, "Skipping negative leading offset");
    } else if (offsets[0]) {
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
        int half = (length + 1) / 2;
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
                    // spoke negated here because increasing i goes left
                    rim ^= RIM_INDEX(offsets[spoke], -spoke, length);
                    // add 1 to spoke here because unlike A we are not symmetric about 0
                    if (spoke) rim ^= RIM_INDEX(NEG(offsets[spoke]), spoke, length);
                }
                ludebug(dbg, "Out of options, so backtrack to spoke %d (rim %d)", spoke, rim);
            } else {
                offsets[spoke] = offset;
                int addition1 = RIM_INDEX(offset, -spoke, length);
                int addition2 = RIM_INDEX(NEG(offset), spoke, length);
                if ((!spoke && !(rim & addition1)) || (spoke && !((rim & addition1) | (rim & addition2)))) {
                    if (spoke == half-1) {
                        ludebug(dbg, "Spoke(s) made rim complete");
                        count += candidate_a(offsets, length);
                        // we never added rim to spoke, so just continue
                    } else {
                        if (spoke) rim |= (addition1 | addition2);
                        else rim |= addition1;
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

    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= offsets[half - i - 1];}
    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= NEG(offsets[i]);}

    ludebug(dbg, "Candidate B length %d offsets %d %d %d -> %x", length, offsets[0], offsets[1], offsets[2], pattern);

    if (GET_SIEVE(pattern)) {
        ludebug(dbg, "Pattern %x already exists", pattern);
    } else if (!offsets[half-1]) {
        ludebug(dbg, "Skipping zero leading offset");
    } else if (offsets[half-1] > UNUSED_OFFSET) {
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
                    // spoke negated here because increasing i goes left
                    rim ^= RIM_INDEX(offsets[spoke], -spoke, length);
                    // add 1 to spoke here because unlike A we are not symmetric about 0
                    rim ^= RIM_INDEX(NEG(offsets[spoke]), spoke + 1, length);
                }
                ludebug(dbg, "Out of options, so backtrack to spoke %d (rim %d)", spoke, rim);
            } else {
                offsets[spoke] = offset;
                int addition1 = RIM_INDEX(offset, -spoke, length);
                int addition2 = RIM_INDEX(NEG(offset), spoke + 1, length);
                if (!((rim & addition1) | (rim & addition2))) {
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
    PATTERN_T pattern = 0;

    // remove trailing 0
//    while (!offsets[0]) {
//        offsets++;
//        length--;
//    }

    for (int i = 0; i < length; ++i) {pattern <<= OFFSET_BITS; pattern |= offsets[length - i - 1];}

    ludebug(dbg, "Candidate C length %d offsets %d %d %d %d %d %d -> %x", length,
            offsets[0], offsets[1], offsets[2], offsets[3], offsets[4], offsets[5], pattern);

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
    } else if (GET_SIEVE(pattern)) {
        ludebug(dbg, "Pattern %x already exists", pattern);
    } else if (!offsets[length-1]) {
        ludebug(dbg, "Skipping zero leading offset");
    } else if (offsets[length-1] > UNUSED_OFFSET) {
        ludebug(dbg, "Skipping negative leading offset");
    } else {
        for (int i = 0; i < MAX_LENGTH - length + 1; ++i) {
            PATTERN_T padded = pattern << (i * OFFSET_BITS);
            if (!GET_SIEVE(padded) && check_lacing(padded, length + i)) {
                write_pattern(offsets, length, 'C', i, length + i);
                count++;
                set_sieve_all(padded, length + i);
            }
        }
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
                // spoke negated here because increasing i goes left
                if (spoke >= 0) rim ^= RIM_INDEX(offsets[spoke], -spoke, length);
                ludebug(dbg, "Out of options, so backtrack to spoke %d (rim %d)", spoke, rim);
            } else {
                offsets[spoke] = offset;
                int addition = RIM_INDEX(offset, -spoke, length);
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
    luinfo(dbg, "%s -h     display this message", progname, PATTERN_FILE);
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
        luinfo(dbg, "Sieve size %dkB (%d entries)", SIEVE_LEN_BYTES / 1024, SIEVE_LEN);
        LU_ALLOC(dbg, sieve, SIEVE_LEN)
        LU_ASSERT(!lufle_exists(dbg, PATTERN_FILE), LU_ERR_IO, dbg, "Output file %s already exists", PATTERN_FILE)
        lufle_open(dbg, PATTERN_FILE, "w", &out);

        flag_unused();
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
