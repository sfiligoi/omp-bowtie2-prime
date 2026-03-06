/*
 * Include this file and use its functions in 
 *  SwAligner::alignEnd2EndSseU8
 * befor and after calling
 *  EEU8_alignNucleotides
 * to save the relevant buffers and contants.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

static pthread_mutex_t aee8_file_mutex = PTHREAD_MUTEX_INITIALIZER;
static uint32_t aee8_file_cnt = 0;

static void appendArraysToBinary(
		                        const SSERegI profbuf[], const size_t profbuf_size,
                                        const char   rf[], const size_t rfd, const size_t mat_size,
                                        const size_t iter, const size_t colstride, const size_t lastWordIdx,
                                        const int32_t minsc, const size_t nrow,
                                        const int8_t refGapOpen, const int8_t refGapExtend, const int8_t readGapOpen, const int8_t readGapExtend) {
    char filename[256];
    sprintf(filename,"EEU8_an_in_%i.bin",int(sizeof(SSERegI)));
    // Open file in "ab" mode: Append and Binary
    FILE *file = fopen(filename, "ab");

    const uint32_t magic = 0x12943ef2;
    const uint32_t isize = sizeof(SSERegI);
    const int8_t gaps[4] = {refGapOpen, refGapExtend, readGapOpen, readGapExtend};
    // ignore errors
    fwrite(&magic,sizeof(uint32_t), 1, file);
    fwrite(&isize,sizeof(uint32_t), 1, file);
    fwrite(&aee8_file_cnt,sizeof(uint32_t), 1, file);
    fwrite(&nrow,sizeof(size_t), 1, file);
    fwrite(&iter,sizeof(size_t), 1, file);
    fwrite(&colstride,sizeof(size_t), 1, file);
    fwrite(&lastWordIdx,sizeof(size_t), 1, file);
    fwrite(&minsc,sizeof(int32_t), 1, file);
    fwrite(&mat_size,sizeof(size_t), 1, file);
    fwrite(&profbuf_size,sizeof(size_t), 1, file);
    fwrite(profbuf, sizeof(SSERegI), profbuf_size, file);
    fwrite(&rfd,sizeof(size_t), 1, file);
    fwrite(rf, sizeof(char), rfd, file);
    fwrite(gaps, sizeof(int8_t), 4, file);
    fclose(file);
}

static void appendRes(
					const int32_t lrmax,  const int32_t btnfilled) {
    char filename[256];
    sprintf(filename,"EEU8_an_res_%i.bin",int(sizeof(SSERegI)));
    FILE *file = fopen(filename, "ab");
    fwrite(&aee8_file_cnt,sizeof(uint32_t), 1, file);
    fwrite(&lrmax,sizeof(int32_t), 1, file);
    fwrite(&btnfilled,sizeof(int32_t), 1, file);
    fclose(file);
}

static void appendArraysToBinary2(
					const int32_t lrmax,  const int32_t btnfilled,
		                        const SSERegI mat[],
                                        const size_t mat_size,
                                        const size_t colstride,
                                        const size_t rfd, const size_t nrow) {
    char filename[256];
    sprintf(filename,"EEU8_an_mat_%i.bin",int(sizeof(SSERegI)));
    // Open file in "ab" mode: Append and Binary
    FILE *file = fopen(filename, "ab");

    const uint32_t magic = 0x12a32ef3;
    const uint32_t isize = sizeof(SSERegI);
    // ignore errors
    fwrite(&magic,sizeof(uint32_t), 1, file);
    fwrite(&isize,sizeof(uint32_t), 1, file);
    fwrite(&aee8_file_cnt,sizeof(uint32_t), 1, file);
    fwrite(&lrmax,sizeof(int32_t), 1, file);
    fwrite(&btnfilled,sizeof(int32_t), 1, file);
    fwrite(&nrow,sizeof(size_t), 1, file);
    fwrite(&rfd,sizeof(size_t), 1, file);
    fwrite(&colstride,sizeof(size_t), 1, file);
    fwrite(&mat_size,sizeof(size_t), 1, file);
    fwrite(mat, sizeof(SSERegI), mat_size, file);
    fclose(file);
}


