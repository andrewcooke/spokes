
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


// these can be changed, but going much deeper requires too much
// memory for the sieve.
#define OFFSET_BITS 3
#define MAX_LENGTH 6

#define PATTERN_FILE "patterns.txt"

// derived sizes
#define PATTERN_BITS (OFFSET_BITS * MAX_LENGTH)
#define UNUSED_OFFSET (1L << (OFFSET_BITS - 1))
#define OFFSET_SIGN UNUSED_OFFSET
#define MAX_OFFSET (UNUSED_OFFSET - 1)
#define OFFSET_VALUE MAX_OFFSET
#define OFFSET_LIMIT (1L << OFFSET_BITS)
#define LENGTH_A ((MAX_LENGTH + 1) / 2)
#define LENGTH_B (MAX_LENGTH / 2)
#define LENGTH_C MAX_LENGTH
#define PATTERN_RIGHT_MASK (OFFSET_LIMIT - 1)
#define PATTERN_LEFT_MASK (((1L << PATTERN_BITS) - 1) ^ PATTERN_RIGHT_MASK)
#define LEFT_ROTATION OFFSET_BITS
#define RIGHT_ROTATION (PATTERN_BITS - LEFT_ROTATION)

// all these are likely best as uint64 for fast access
#define SIEVE_T uint64_t
#define OFFSET_T uint64_t   // minimum length OFFSET_BITS
#define PATTERN_T uint64_t  // minimum length PATTERN_BITS
#define HOLES_T uint64_t  // minimum length MAX_LENGTH

// global state.  for a single-minded, math-intensive program it's
// pointless to pass these around as arguments
SIEVE_T *sieve = NULL;
lulog *dbg = NULL;
FILE *out = NULL;

// more derived sizes
#define SIEVE_WIDTH (8 * sizeof(*sieve))
#define SIEVE_LEN_BITS (1L << (PATTERN_BITS - 1))
#define SIEVE_LEN ((SIEVE_LEN_BITS + (SIEVE_WIDTH - 1)) / SIEVE_WIDTH)
#define SIEVE_LEN_BYTES (SIEVE_LEN * 8)

#define SIEVE_INDEX(n) (n / SIEVE_WIDTH)
#define SIEVE_SHIFT(n) (n - SIEVE_WIDTH * SIEVE_INDEX(n))
#define SET_SIEVE(n) sieve[SIEVE_INDEX(n)] |= (1L << SIEVE_SHIFT(n))
#define GET_SIEVE(n) 1 & (sieve[SIEVE_INDEX(n)] >> SIEVE_SHIFT(n))

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
//void set_sieve(int n) {
//    int index = SIEVE_INDEX(n);
//    int shift = SIEVE_SHIFT(n);
//    SET_SIEVE(n);
//    ludebug(dbg, "Pattern %d -> sieve set at %d/%d", n, index, shift);
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

void set_sieve_all(PATTERN_T pattern) {
    for (int i = 0; i < MAX_LENGTH; ++i) {
        SET_SIEVE(pattern);
        pattern = ((pattern & PATTERN_LEFT_MASK) >> LEFT_ROTATION) | ((pattern & PATTERN_RIGHT_MASK) << RIGHT_ROTATION);
    }
}

void write_pattern_a(OFFSET_T *offsets, int length, int padding) {

    char buffer[3*MAX_LENGTH+3], *p;

    p = buffer;
    for (int i = 0; i < length; ++i) {
        if (i) p += sprintf(p, ",");
        if (offsets[length - i - 1] > UNUSED_OFFSET) {
            p += sprintf(p, "-%d", offsets[length - i - 1] - UNUSED_OFFSET);
        } else {
            p += sprintf(p, "%d", offsets[length - i - 1]);
        }
    }
    *(p++) = 'A';
    if (padding) p += sprintf(p, "%d", padding);
    *p = '\0';

    luinfo(dbg, "Writing %s", buffer);
    fprintf(out, "%s\n", buffer);
}

int candidate_a(OFFSET_T *offsets, int length) {

    int count = 0;
    int half = (length + 1) / 2;
    PATTERN_T pattern = 0;

    ludebug(dbg, "Candidate length %d offsets %d %d %d", length, offsets[0], offsets[1], offsets[2]);

    for (int i = 0; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= offsets[half - i - 1];}
    for (int i = 1; i < half; ++i) {pattern <<= OFFSET_BITS; pattern |= (offsets[i] | OFFSET_SIGN);}

    if (GET_SIEVE(pattern)) {
        ludebug(dbg, "Pattern %x already exists", pattern);
    } else if (offsets[half-1] > UNUSED_OFFSET) {
        ludebug(dbg, "Skipping negative leading offset");
    } else if (offsets[0]) {
        ludebug(dbg, "Skipping non-radial central spoke");
    } else {
        if (length == 1 && offsets[0]) {
            ludebug(dbg, "Unbalanced %d %d", length, pattern);
        } else {
            for (int i = 0; i < MAX_LENGTH - length + 1; ++i) {
                write_pattern_a(offsets, half, i);
                count++;
            }
        }
        set_sieve_all(pattern);
    }

    return count;
}

void search_a() {

    luinfo(dbg, "Searching for A group patterns");

    // start with single radial spoke
    OFFSET_T offsets[LENGTH_A] = {0};
    int i = 0, offset = 0, length = 1, count = 0;
    PATTERN_T pattern = 0;
    HOLES_T rim = 1, hub = rim;

    // run through all possible patterns
    while (length < MAX_LENGTH) {

        count += candidate_a(offsets, length);

        // remove current spoke(s) from rim
        offset = offsets[i];
        // i negated here because increasing i goes left
        rim ^= RIM_INDEX(offset, -i, length);
        if (i) rim ^= RIM_INDEX(-offset, i, length);
        ludebug(dbg, "Rim after removal %d", rim);

        // search for next lacing
        do {

            // increment
            offset = offsets[i] + 1;
            if (offset == UNUSED_OFFSET) offset++;
            if (offset == OFFSET_LIMIT) offset = 0;   // carry
            offsets[i] = offset;
            ludebug(dbg, "Increment at index %d to %d", i, offset);

            if (!offset) {

                // if we carried, then we need to increment the next level
                // up before we start testing spokes
                i++;
                ludebug(dbg, "Index increased to %d", i);
                if (2 * i + 1 > length) {
                    // this may mean that we are now considering a longer length
                    if (rim != 0) {luerror(dbg, "Non-zero rim!"); return;}
                    length = 2 * i + 1;
                    hub = (hub << 2) | 3;
                    ludebug(dbg, "New length %d, rim %d, hub %d", length, rim, hub);
                } else {
                    // otherwise, we need to reset lower spokes
                    for (int j = i; j > 0; --j) offsets[j-1] = 0;
                    // and remove the one we will increment
                    offset = offsets[i];
                    ludebug(dbg, "Index %d offset %d", i, offset);
                    rim ^= RIM_INDEX(offset, -i, length);
                    if (i) rim ^= RIM_INDEX(-offset, i, length);
                    ludebug(dbg, "Rim after removal (new index) %d", rim);
                }

            } else {

                // if we can lace new spokes do so, until all laced or we
                // have a new point to increment from
                int ok = 1;
                while (ok) {
                    offset = offsets[i];
                    // test if lacing ok
                    int test = rim, addition = RIM_INDEX(offset, -i, length);
                    if (test & addition) {
                        ok = 0;
                    } else {
                        test |= addition;
                        if (i) {
                            addition = RIM_INDEX(-offset, i, length);
                            if (test & addition) {
                                ok = 0;
                            } else {
                                test |= addition;
                            }
                        }
                    }
                    if (ok) {
                        rim = test;
                        ludebug(dbg, "Successful lace at index %d, rim %d", i, rim);
                        if (i) {
                            ludebug(dbg, "Decrementing index");
                            i--;
                        } else {
                            ok = 0;
                        }
                    } else {
                        ludebug(dbg, "Failed to lace");
                    }
                }

            }

        } while (rim != hub && length < MAX_LENGTH);

    }

    luinfo(dbg, "Found %d A group patterns", count);
}

void usage(const char *progname) {
    luinfo(dbg, "Search for spoke patterns (max offset %d, max length %d)", MAX_OFFSET, MAX_LENGTH);
    luinfo(dbg, "%s -h     display this message", progname, PATTERN_FILE);
    luinfo(dbg, "%s        run a search (output to %s)", progname, PATTERN_FILE);
}

// error handling is for lulib routines; don't bother elsewhere.
int main(int argc, char** argv) {

    LU_STATUS

    lulog_mkstderr(&dbg, lulog_level_debug);

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

    }

LU_CLEANUP
    free(sieve);
    if (out) fclose(out);
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}
