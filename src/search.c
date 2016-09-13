
#include <stdint.h>
#include <stdio.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#define OFFSET_BITS 3
#define MAX_LENGTH 6

#define LENGTH_A ((MAX_LENGTH - 1) / 2)
#define LENGTH_B (MAX_LENGTH / 2)
#define LENGTH_C MAX_LENGTH

#define UNUSED_OFFSET (1 << (OFFSET_BITS - 1))
#define MAX_OFFSET (UNUSED_OFFSET - 1)
#define OFFSET_LIMIT (1 << OFFSET_BITS)
#define PATTERN_BITS (OFFSET_BITS * MAX_LENGTH)

// if you change this, also change 63, 64 below
#define SIEVE uint64_t
#define SIEVE_BITS 64

#define SIEVE_LEN_BITS (1 << (PATTERN_BITS - 1))
#define SIEVE_LEN ((SIEVE_LEN_BITS + (SIEVE_BITS - 1)) / SIEVE_BITS)
#define SIEVE_LEN_BYTES (SIEVE_LEN * 8)

#define SIEVE_INDEX(pattern) (pattern / SIEVE_BITS)
#define SIEVE_MASK(pattern) (1 << (pattern - SIEVE_INDEX(pattern) * SIEVE_BITS))

#define SIEVE_RIGHT_MASK (OFFSET_LIMIT - 1)

// this must contain at least OFFSET_BITS bits (unsigned)
#define OFFSET uint64_t
#define OFFSET_MASK ((1 << OFFSET_BITS) - 1)
#define OFFSET_SIGN UNUSED_OFFSET
#define OFFSET_VALUE (OFFSET_MASK ^ OFFSET_SIGN)
#define GET_OFFSET(pattern, index) ((pattern >> (index * OFFSET_BITS)) & OFFSET_MASK)

// this must contain at least PATTERN_BITS bits (unsigned)
#define PATTERN uint64_t

#define PATTERN_FILE "patterns.txt"


int check_sieve(SIEVE *sieve, PATTERN pattern) {
    int index = SIEVE_INDEX(pattern);
    SIEVE mask = SIEVE_MASK(pattern);
    return sieve[index] & mask;
}

// we update the sieve for all possible rotations.  this implies that we
// must generate the "least padded" pattern first and then write all
// related paddings.
void update_sieve(SIEVE *sieve, PATTERN pattern, int length) {
    int right_rotn = OFFSET_BITS * (length - 1);
    PATTERN left_mask = ((1 << length) - 1) ^ SIEVE_RIGHT_MASK;
    for (int i = 0; i < length; ++i) {
        int index = SIEVE_INDEX(pattern);
        SIEVE mask = SIEVE_MASK(pattern);
        sieve[index] |= mask;
        pattern = ((pattern & left_mask) >> OFFSET_BITS) | ((pattern & SIEVE_RIGHT_MASK) << right_rotn);
    }
}

void print_result(lulog *log, FILE *out, PATTERN pattern, const char *name) {
    ludebug(log, "%08x %s", pattern, name);
    fprintf(out, "%s\n", name);
}

char *sprint_offset(char *p, OFFSET offset) {
    if (offset & OFFSET_SIGN) sprintf(p++, "-");
    sprintf(p++, "%d", offset & OFFSET_VALUE);
    return p;
}

void print_result_a(lulog *log, FILE *out, PATTERN pattern, int padding) {
    char name[2*MAX_LENGTH+3] = {0};
    char *p = name;
    PATTERN mask = OFFSET_MASK << OFFSET_BITS;
    int half = 1;
    while (pattern & mask) {half++; mask <<= OFFSET_BITS;}
    for (int i = half; i > 0; --i) {
        p = sprint_offset(p, GET_OFFSET(pattern, i-1));
    }
    sprintf(p++, "A");
    if (padding) sprintf(p, "%d", padding);
    print_result(log, out, pattern, name);
}

void radial(lulog *log, SIEVE *sieve, FILE *out) {
    update_sieve(sieve, 0, 1);
    print_result_a(log, out, 0, 0);
}

void search_a(lulog *log, SIEVE *sieve) {
    PATTERN pattern = 0;
    while (1) {
        int index = 0;
        do {
            OFFSET offset = GET_OFFSET(pattern, index) + 1;
            if (offset == UNUSED_OFFSET) offset++;

        } while (1);
        if (!check_sieve(sieve, pattern)) {

        }
    }
}

int search_b(lulog *log, SIEVE *sieve) {
    return LU_OK;
}

int search_c(lulog *log, SIEVE *sieve) {
    return LU_OK;
}

int main(int argc, char** argv) {

    LU_STATUS
    lulog *log = NULL;
    FILE *out = NULL;
    SIEVE *sieve = NULL;

    lulog_mkstderr(&log, lulog_level_debug);
    luinfo(log, "Maximum spoke offset %d; Maximum pattern length %d", MAX_OFFSET, MAX_LENGTH);
    luinfo(log, "Sieve size %dkB", SIEVE_LEN_BYTES / 1024);
    LU_ALLOC(log, sieve, SIEVE_LEN)
    LU_ASSERT(!lufle_exists(log, PATTERN_FILE), LU_ERR_IO, log, "Output file %s already exists", PATTERN_FILE)
    lufle_open(log, PATTERN_FILE, "w", &out);

    radial(log, sieve, out);
    search_a(log, sieve);
    LU_CHECK(search_b(log, sieve))
    LU_CHECK(search_c(log, sieve))

LU_CLEANUP
    if (out) fclose(out);
    if (log) status = log->free(&log, status);
    return status;
}
