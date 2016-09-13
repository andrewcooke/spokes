
#include <stdint.h>
#include <stdio.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#define OFFSET_BITS 3
#define MAX_DEPTH 6

#define UNUSED_OFFSET (1 << (OFFSET_BITS - 1))
#define MAX_OFFSET (UNUSED_OFFSET - 1)
#define OFFSET_LIMIT (1 << OFFSET_BITS)
#define PATTERN_BITS (OFFSET_BITS * MAX_DEPTH)

// if you change this, also change 63, 64 below
#define SIEVE uint64_t
#define SIEVE_BITS 64

#define SIEVE_LEN_BITS (1 << (PATTERN_BITS - 1))
#define SIEVE_LEN ((SIEVE_LEN_BITS + (SIEVE_BITS - 1)) / SIEVE_BITS)
#define SIEVE_LEN_BYTES (SIEVE_LEN * 8)

#define SIEVE_INDEX(pattern) (pattern / SIEVE_BITS)
#define SIEVE_MASK(pattern) (1 << (pattern - SIEVE_INDEX(pattern) * SIEVE_BITS))

#define SIEVE_RIGHT_MASK (OFFSET_LIMIT - 1)
#define SIEVE_LEFT_MASK (((1 << PATTERN_BITS) - 1) ^ SIEVE_RIGHT_MASK)
#define SIEVE_RIGHT_ROTN (PATTERN_BITS - OFFSET_BITS)
#define SIEVE_LEFT_ROTN OFFSET_BITS

// this must contain at least OFFSET_BITS bits (unsigned)
#define SPOKE uint64_t

// this must contain at least OFFSET_BITS * MAX_DEPTH bits (unsigned)
#define PATTERN uint64_t

#define PATTERN_FILE "patterns.txt"

#define K 1024
#define M (K * K)


int check_sieve(SIEVE *sieve, PATTERN pattern) {
    int index = SIEVE_INDEX(pattern);
    SIEVE mask = SIEVE_MASK(pattern);
    return sieve[index] & mask;
}

void update_sieve(SIEVE *sieve, PATTERN pattern) {
    for (int i = 0; i < MAX_OFFSET; ++i) {
        int index = SIEVE_INDEX(pattern);
        SIEVE mask = SIEVE_MASK(pattern);
        sieve[index] |= pattern;
        pattern = ((pattern & SIEVE_LEFT_MASK) >> SIEVE_LEFT_ROTN) | ((pattern & SIEVE_RIGHT_MASK) << SIEVE_RIGHT_ROTN);
    }
}

void print_result(FILE *out, PATTERN pattern, const char *name) {
    fprintf(out, "%08x %s\n", pattern, name);
}

int radial(lulog *log, SIEVE *sieve, FILE *out) {
    update_sieve(sieve, 0);
    print_result(out, 0, "A0");
    return LU_OK;
}

int search_a(lulog *log, SIEVE *sieve) {
    for (SPOKE i = 0; i < OFFSET_LIMIT; ++i) {
        if (i != UNUSED_OFFSET) {
            ludebug(log, "%d", i);
        }
    }
    return LU_OK;
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
    luinfo(log, "Maximum offset %d.  Maximum depth %d", MAX_OFFSET, MAX_DEPTH);
    ludebug(log, "Sieve %dB (%dkB, %dMB)", SIEVE_LEN_BYTES, SIEVE_LEN_BYTES / K, SIEVE_LEN_BYTES / M);
    LU_ALLOC(log, sieve, SIEVE_LEN)
    LU_ASSERT(!lufle_exists(log, PATTERN_FILE), LU_ERR_IO, log, "Output file %s already exists", PATTERN_FILE)
    lufle_open(log, PATTERN_FILE, "w", &out);

    LU_CHECK(radial(log, sieve, out))
    LU_CHECK(search_a(log, sieve))
    LU_CHECK(search_b(log, sieve))
    LU_CHECK(search_c(log, sieve))

LU_CLEANUP
    if (out) fclose(out);
    if (log) status = log->free(&log, status);
    return status;
}
