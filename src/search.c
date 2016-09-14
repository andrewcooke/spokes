
#include <stdint.h>
#include <stdio.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"


// these can be changed, but i doubt going any deeper is going to be
// efficient - time is probably exponential in their product
#define OFFSET_BITS 3
#define MAX_LENGTH 6

#define PATTERN_FILE "patterns.txt"

// derived sizes
#define PATTERN_BITS (OFFSET_BITS * MAX_LENGTH)
#define UNUSED_OFFSET (1 << (OFFSET_BITS - 1))
#define MAX_OFFSET (UNUSED_OFFSET - 1)
#define OFFSET_LIMIT (1 << OFFSET_BITS)
#define LENGTH_A ((MAX_LENGTH + 1) / 2)
#define LENGTH_B (MAX_LENGTH / 2)
#define LENGTH_C MAX_LENGTH

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
#define SIEVE_LEN_BITS (1 << (PATTERN_BITS - 1))
#define SIEVE_LEN ((SIEVE_LEN_BITS + (SIEVE_WIDTH - 1)) / SIEVE_WIDTH)
#define SIEVE_LEN_BYTES (SIEVE_LEN * 8)

#define SIEVE_INDEX(n) (n / SIEVE_WIDTH)
#define SIEVE_SHIFT(n) (n - SIEVE_WIDTH * SIEVE_INDEX(n))
#define SET_SIEVE(n) sieve[SIEVE_INDEX(n)] |= 1 << SIEVE_SHIFT(n)
#define GET_SIEVE(n) 1 & (sieve[SIEVE_INDEX(n)] >> SIEVE_SHIFT(n))

#define RIM_INDEX(offset, index, length) (1 << ((index + offset + length + length) % length))


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
//            ludebug(dbg, "Setting %x (%d)", pattern, count);
            SET_SIEVE(pattern);
            count++;
        }

    }

    luinfo(dbg, "Set %d (/%d = %2.0f%%)  entries", count, SIEVE_LEN_BITS, count * 100.0 / SIEVE_LEN_BITS);
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

        // TODO - sieve
        // TODO - write output here
        ludebug(dbg, "Pattern? length %d offsets %d %d %d rim %d", length, offsets[0], offsets[1], offsets[2], rim);

        // remove current spoke(s) from rim
        offset = offsets[i];
        rim ^= RIM_INDEX(offset, i, length);
        if (i) rim ^= RIM_INDEX(-offset, -i, length);
        ludebug(dbg, "rim after removal %d", rim);

        // search for next lacing
        do {

            // increment
            offset = offsets[i] + 1;
            if (offset == UNUSED_OFFSET) offset++;
            if (offset == OFFSET_LIMIT) offset = 0;  // carry
            offsets[i] = offset;
            ludebug(dbg, "Increment at index %d to %d", i, offset);

            if (!offset) {

                // if we carried, then we need to increment the next level
                // up before we start testing spokes
                i++;
                ludebug(dbg, "Index increased to %d", i);
                // this may mean that we are now considering a longer length
                if (2 * i + 1 > length) {
                    if (rim != 0) {luerror(dbg, "Non-zero rim!"); return;}
                    if (offsets[i] != 0) {luerror(dbg, "Non-zero offset!"); return;}
                    length = 2 * i + 1;
                    hub = (hub << 2) | 3;
                    ludebug(dbg, "New length %d, rim %d, hub %d", length, rim, hub);
                } else {
                    offset = offsets[i];
                    rim ^= RIM_INDEX(offset, i, length);
                    if (i) rim ^= RIM_INDEX(-offset, -i, length);
                    ludebug(dbg, "Rim after removal %d", rim);
                }

            } else {

                // if we can lace new spokes do so, until all laced or we
                // have a new point to increment from
                int ok = 1;
                while (ok) {
                    offset = offsets[i];
                    // test if lacing ok
                    int test = rim, addition = RIM_INDEX(offset, i, length);
                    if (test & addition) {
                        ok = 0;
                    } else {
                        test |= addition;
                        if (i) {
                            addition = RIM_INDEX(-offset, -i, length);
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
                            ludebug(dbg, "Exiting: rim %d, hub %d", rim, hub);
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

// error handling is for lulib routines; don't bother elsewhere.
int main(int argc, char** argv) {

    LU_STATUS

    lulog_mkstderr(&dbg, lulog_level_debug);
    luinfo(dbg, "Maximum spoke offset %d; Maximum pattern length %d", MAX_OFFSET, MAX_LENGTH);
    luinfo(dbg, "Sieve size %dkB (%d entries)", SIEVE_LEN_BYTES / 1024, SIEVE_LEN);
    LU_ALLOC(dbg, sieve, SIEVE_LEN)
    LU_ASSERT(!lufle_exists(dbg, PATTERN_FILE), LU_ERR_IO, dbg, "Output file %s already exists", PATTERN_FILE)
    lufle_open(dbg, PATTERN_FILE, "w", &out);

    flag_unused();
    search_a();

LU_CLEANUP
    free(sieve);
    if (out) fclose(out);
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}

//
//
//
//
//
//
//
//
//
//
//#define LENGTH_A ((MAX_LENGTH - 1) / 2)
//#define LENGTH_B (MAX_LENGTH / 2)
//#define LENGTH_C MAX_LENGTH
//
//#define UNUSED_OFFSET (1 << (OFFSET_BITS - 1))
//#define MAX_OFFSET (UNUSED_OFFSET - 1)
//#define OFFSET_LIMIT (1 << OFFSET_BITS)
//#define PATTERN_BITS (OFFSET_BITS * MAX_LENGTH)
//
//// if you change this, also change 63, 64 below
//#define SIEVE uint64_t
//#define SIEVE_BITS 64
//
//#define SIEVE_LEN_BITS (1 << (PATTERN_BITS - 1))
//#define SIEVE_LEN ((SIEVE_LEN_BITS + (SIEVE_BITS - 1)) / SIEVE_BITS)
//#define SIEVE_LEN_BYTES (SIEVE_LEN * 8)
//
//#define SIEVE_INDEX(pattern) (pattern / SIEVE_BITS)
//#define SIEVE_MASK(pattern) (1 << (pattern - SIEVE_INDEX(pattern) * SIEVE_BITS))
//
//#define SIEVE_RIGHT_MASK (OFFSET_LIMIT - 1)
//
//// this must contain at least OFFSET_BITS bits (unsigned)
//#define OFFSET uint64_t
//#define OFFSET_MASK ((1 << OFFSET_BITS) - 1)
//#define OFFSET_SIGN UNUSED_OFFSET
//#define OFFSET_VALUE (OFFSET_MASK ^ OFFSET_SIGN)
//#define GET_OFFSET(pattern, index) ((pattern >> (index * OFFSET_BITS)) & OFFSET_MASK)
//#define SET_OFFSET(pattern, index, offset) (pattern | (offset << (index * OFFSET_BITS)))
//
//// this must contain at least PATTERN_BITS bits (unsigned)
//#define PATTERN uint64_t
//
//// this must contain at least MAX_LENGTH bits (unsigned)
//#define HOLES uint64_t
//
//#define PATTERN_FILE "patterns.txt"
//
//
//int check_sieve(SIEVE *sieve, PATTERN pattern) {
//    int index = SIEVE_INDEX(pattern);
//    SIEVE mask = SIEVE_MASK(pattern);
//    return sieve[index] & mask;
//}
//
//// we update the sieve for all possible rotations.  this implies that we
//// must generate the "least padded" pattern first and then write all
//// related paddings.
//void update_sieve(SIEVE *sieve, PATTERN pattern, int length) {
//    int right_rotn = OFFSET_BITS * (length - 1);
//    PATTERN left_mask = ((1 << length) - 1) ^ SIEVE_RIGHT_MASK;
//    for (int i = 0; i < length; ++i) {
//        int index = SIEVE_INDEX(pattern);
//        SIEVE mask = SIEVE_MASK(pattern);
//        sieve[index] |= mask;
//        pattern = ((pattern & left_mask) >> OFFSET_BITS) | ((pattern & SIEVE_RIGHT_MASK) << right_rotn);
//    }
//}
//
//void print_result(lulog *log, FILE *out, PATTERN pattern, const char *name) {
//    ludebug(log, "%08x %s", pattern, name);
//    fprintf(out, "%s\n", name);
//}
//
//char *sprint_offset(char *p, OFFSET offset) {
//    if (offset & OFFSET_SIGN) sprintf(p++, "-");
//    sprintf(p++, "%d", offset & OFFSET_VALUE);
//    return p;
//}
//
//void print_result_a(lulog *log, FILE *out, PATTERN pattern, int padding) {
//    char name[2*MAX_LENGTH+3] = {0};
//    char *p = name;
//    PATTERN mask = OFFSET_MASK << OFFSET_BITS;
//    int half = 1;
//    while (pattern & mask) {half++; mask <<= OFFSET_BITS;}
//    for (int i = half; i > 0; --i) {
//        p = sprint_offset(p, GET_OFFSET(pattern, i-1));
//    }
//    sprintf(p++, "A");
//    if (padding) sprintf(p, "%d", padding);
//    print_result(log, out, pattern, name);
//}
//
//void radial(lulog *log, SIEVE *sieve, FILE *out) {
//    update_sieve(sieve, 0, 1);
//    print_result_a(log, out, 0, 0);
//}
//
//void search_a_half(lulog *log, SIEVE *sieve, FILE *out, int half) {
//
//    // start with full radial and increment "outer" spoke
//    PATTERN pattern = 0;
//    HOLES rim = (1 << (2 * half + 1)) - 1;
//    HOLES hub = rim;
//    int index = half - 1;
//
//    while (1) {
//
//        // increment with carry
//        OFFSET offset = 0;
//        while (index < half && !offset) {
//            offset = GET_OFFSET(pattern, index) + 1;
//            if (offset == UNUSED_OFFSET) offset++;
//            offset = OFFSET & OFFSET_MASK;
//            pattern = SET_OFFSET(pattern, index, offset);
//            index++;
//        }
//
//        // overflow
//        if (index == half) return;
//
//        // expand to full pattern here
//        if (!check_sieve(sieve, pattern)) {
//            update_sieve_a(sieve, pattern);
//            print_result_a(log, out, pattern);
//        }
//
//        // start increment from "inner" spoke after initial bootstrap
//        index = 0;
//    }
//}
//
//void search_a(lulog *log, SIEVE *sieve, FILE *out) {
//    for (int i = 0; i < LENGTH_A; ++i) {
//        search_a_half(log, sieve, out, i+1);
//    }
//}
//
//int search_b(lulog *log, SIEVE *sieve) {
//    return LU_OK;
//}
//
//int search_c(lulog *log, SIEVE *sieve) {
//    return LU_OK;
//}
//
