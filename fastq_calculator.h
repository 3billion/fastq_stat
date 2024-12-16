#ifndef FASTQSTAT_CALC_H_
#define FASTQSTAT_CALC_H_

#include <atomic>
#include <fstream>
#include <immintrin.h>
#include <iostream>
#include <string.h>
#include <thread>
#include <vector>
struct Stat {
    long linenum = 0;
    long yield = 0;
    long q20_base = 0;
    long q30_base = 0;
};

class FastqCalculator {
  public:
    FastqCalculator() {}
    ~FastqCalculator() {
        if (buffer != nullptr) {
            delete[] buffer;
        }
    }
    bool Calculate(const char *buffer, long buffer_size, int thread_count,
                   Stat &stat_file);
    void Free() {
        if (buffer != nullptr) {
            delete[] buffer;
            buffer = nullptr;
        }
    }
    char *buffer = nullptr;
};

void correct_pos(const char *buffer, long *const cur_end_pos,
                 long *const next_start_pos) {
    while (buffer[*next_start_pos] != '@') {
        (*next_start_pos)++;
    }
    *cur_end_pos = *next_start_pos;
}

void find_qc_pos_job(const char *buffer, long start_pos, long end_pos,
                     Stat &stat) {

    // avx-512 code, seq라인을 읽어서 score line의 length를 읽어들이고, 그 길이만큼 simd명령을 수행한다.
    int line_index =0;
    __m512i q20 = _mm512_set1_epi8(20 + 33);
    __m512i q30 = _mm512_set1_epi8(30 + 33);
    char score_line[64*4];
    int yield =0;

    for(long i = start_pos; i < end_pos; i++) {
        char c = buffer[i];
        if(c == '@') {
            line_index = 0;
            stat.linenum++;
            yield = 0;
            memset(score_line, 0, sizeof(score_line));
        }
        if(c == '\n') {
            line_index++;
            continue;
        } else if(c != 0) {
            if(line_index == 1){
                yield++;
            } 
            if (line_index == 3) {
                memcpy(score_line, buffer + i, yield);
                int yield_64_count = (yield + 63) / 64;
                for(int j = 0; j < yield_64_count; j++) {
                    __m512i score = _mm512_loadu_epi8(score_line + j*64);
                    __mmask64 q20_mask = _mm512_cmpge_epi8_mask(score, q20);
                    __mmask64 q30_mask = _mm512_cmpge_epi8_mask(score, q30);
                    stat.q20_base += _mm_popcnt_u64(q20_mask);
                    stat.q30_base += _mm_popcnt_u64(q30_mask);
                }
                stat.yield += yield;
                i += yield;
            }
        }
    }
    /* 기존 한 글자씩 읽는 코드..
    int line_index = 0;
    for (long i = start_pos; i < end_pos; i++) {
        char c = buffer[i];
        if (c == '@') {
            line_index = 0;
            stat.linenum++;
        }

        if (c == '\n') {
            if (line_index == 3) {
            }
            line_index++;
            continue;
        } else if (c != 0) {
            if (line_index == 3) {
                stat.yield++;
                if (c - 33 >= 20) {
                    stat.q20_base++;
                }
                if (c - 33 >= 30) {
                    stat.q30_base++;
                }
            }
        }
    }
    */
}

bool FastqCalculator::Calculate(const char *buffer, long buffer_size,
                                int thread_count, Stat &stat_file) {
    long *thread_start_pos = new long[thread_count];
    long *thread_end_pos = new long[thread_count];
    long basic_step = buffer_size / thread_count;
    for (int i = 0; i < thread_count; i++) {
        thread_start_pos[i] = i * basic_step;
        thread_end_pos[i] =
            (i == thread_count - 1) ? buffer_size : (i + 1) * basic_step;
    }

    // Correct start and end position
    for (int i = 0; i < thread_count - 1; i++) {
        long before_start_pos = thread_start_pos[i];
        long before_end_pos = thread_end_pos[i];
        correct_pos(buffer, &thread_end_pos[i], &thread_start_pos[i + 1]);
    }

    // Extract QC data from buffer in parallel
    std::thread *threads = new std::thread[thread_count];
    Stat *stat = new Stat[thread_count];
    for (int i = 0; i < thread_count; i++) {
        threads[i] = std::thread(find_qc_pos_job, buffer, thread_start_pos[i],
                                 thread_end_pos[i], std::ref(stat[i]));
    }
    for (int i = 0; i < thread_count; i++) {
        threads[i].join();

        stat_file.linenum += stat[i].linenum;
        stat_file.yield += stat[i].yield;
        stat_file.q20_base += stat[i].q20_base;
        stat_file.q30_base += stat[i].q30_base;
    }

    delete[] thread_start_pos;
    delete[] thread_end_pos;
    delete[] threads;

    return true;
}
#endif // FASTQSTAT_CALC_H_
