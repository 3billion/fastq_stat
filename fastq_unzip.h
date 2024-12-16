#ifndef FASTQSTAT_UNZIP_H_
#define FASTQSTAT_UNZIP_H_
#include "./libdeflate/libdeflate.h"
#include <algorithm>
#include <atomic>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <string.h>
#include <vector>
#include <immintrin.h>

typedef struct gzip_block {
    long index = 0;
    long offset = 0;
    long compressed_size = 0;
    long uncompressed_size = 0;
    char *compressed = 0;
    char *uncompressed = 0;
    long actual_out_nbytes = 0;
} GZIP_BLOCK;

typedef struct {
    const unsigned char *buffer;
    long start;
    long length;
    std::vector<long> offsets;
} OFFSET_THREAD_PARAM;

typedef struct {
    long start;
    long length;
    std::vector<GZIP_BLOCK> *blocks;
    std::atomic<long> *total_size;
} THREAD_PARAM;

typedef struct {
    std::vector<GZIP_BLOCK> *blocks;
    char *buffer;
    long start;
    long end;
} MERGE_THREAD_PARAM;

/**
 * @class FastqUnzip
 * @brief A class to handle the unzipping of FASTQ files.
 * 
 * This class provides functionality to unzip FASTQ files using a specified number of threads.
 * It maintains an uncompressed buffer and tracks the total output size.
 */
class FastqUnzip {
  public:
    FastqUnzip() { total_out_size = 0; };
    ~FastqUnzip() {}
    bool Unzip(const std::string &file_name, int thread_count);
    char *uncompressed_buffer = nullptr;
    std::atomic<long> total_out_size;
};

/**
 * @brief Thread function to find specific patterns in a buffer and store their offsets.
 *
 * This function searches through a buffer for a specific pattern of bytes that 
 * indicate the start of a compressed block (gzip header). When such a pattern 
 * is found, the offset (position) of the pattern in the buffer is stored in 
 * the offsets vector of the parameter structure.
 *
 * @param arg A pointer to an OFFSET_THREAD_PARAM structure containing the 
 *        buffer to search, the start position, the length of the segment to 
 *        search, and the vector to store the offsets.
 * @return Always returns (void *)1.
 */
void *offset_thread(void *arg) {
    OFFSET_THREAD_PARAM *param = (OFFSET_THREAD_PARAM *)arg;
    long i = 0;
    for (i = param->start; i < param->start + param->length; ++i) {
        if (param->buffer[i] == 0x1f && param->buffer[i + 1] == 0x8b &&
            param->buffer[i + 2] == 0x08) {
            if ((param->buffer[i + 3] & 0xE0) ==
                0x0) { /// flags. The most significant 3 bits are reserved for
                       /// 0s.
                if (param->buffer[i + 8] == 0x0 ||
                    param->buffer[i + 8] == 0x2 ||
                    param->buffer[i + 8] ==
                        0x4) { /// extra flags. Only 2 and 4 are valid for
                               /// deflate method (CM=8).
                    if (param->buffer[i + 9] == 0x0 ||
                        param->buffer[i + 9] == 0x3 ||
                        param->buffer[i + 9] == 255) { /// OS.
                        param->offsets.push_back(i);
                    }
                }
            }
        }
    }

    return (void *)1;
}

/**
 * @brief Finds and processes gzip members in a given buffer.
 *
 * This function divides the buffer into chunks and processes each chunk in a separate thread
 * to find gzip members. The offsets of the gzip members are collected and used to create
 * GZIP_BLOCK structures which are stored in the destination vector.
 *
 * @param buffer Pointer to the buffer containing the gzip data.
 * @param file_size Size of the buffer.
 * @param core Number of threads to use for processing.
 * @param dest Pointer to a vector where the found GZIP_BLOCK structures will be stored.
 * @return true if gzip members are found and processed successfully, false otherwise.
 */
bool find_gzip_members(const char *const buffer, const long &file_size,
                       const int &core, std::vector<GZIP_BLOCK> *const dest) {

    long chunk_size = file_size / core;
    long start = 0;
    int i = 0;
    std::vector<OFFSET_THREAD_PARAM *> params;
    std::vector<pthread_t> threads;
    for (i = 0, start = 0; i < core; ++i, start += chunk_size) {
        OFFSET_THREAD_PARAM *p = new OFFSET_THREAD_PARAM;
        p->start = start - 10; /// 3-byte overlap. +7 for more sure.
        if (p->start + chunk_size >= file_size) {
            p->length = file_size - p->start;
        } else {
            p->length = chunk_size;
        }
        p->buffer = (const unsigned char *)buffer;

        pthread_t t;
        int s = pthread_create(&t, NULL, offset_thread, (void *)p);
        if (s != 0) {
            delete p;
        }
        params.push_back(p);
        threads.push_back(t);
    }

    std::vector<pthread_t>::iterator thread;
    for (thread = threads.begin(); thread != threads.end(); ++thread) {
        void *res = NULL;
        pthread_join(*thread, &res);
    }

    std::vector<long> offsets;
    std::vector<OFFSET_THREAD_PARAM *>::iterator param;
    for (param = params.begin(); param != params.end(); ++param) {
        if ((*param)->offsets.empty() == true) {
            delete (*param);
            continue;
        }
        offsets.insert(offsets.begin(), (*param)->offsets.begin(),
                       (*param)->offsets.end());

        delete (*param);
    }
    params.clear();

    //std::cout << "total number of members:" << offsets.size() << std::endl;

    if (offsets.empty() == true) {
        std::cerr << "No offsets found." << std::endl;
        return false;
    }

    if (offsets.size() == 1) {
        long length = file_size - offsets[0];
        unsigned int isize = *(unsigned int *)(buffer + file_size - 4);

        GZIP_BLOCK block;
        block.index = 0;
        block.offset = offsets[0];
        block.compressed_size = length;
        block.uncompressed_size = isize;
        block.compressed = (char *)(buffer + block.offset);
        block.uncompressed = NULL;
        dest->push_back(block);

        return true;
    }

    std::sort(offsets.begin(), offsets.end());
    long index = 0;
    std::vector<long>::const_iterator offset;
    std::vector<long>::const_iterator next_offset;
    for (offset = offsets.begin(), next_offset = offsets.begin() + 1, index = 0;
         next_offset != offsets.end(); ++offset, ++next_offset, ++index) {
        if (*offset == *next_offset)
            continue;

        long length = *next_offset - *offset + 1;
        unsigned int isize = *(unsigned int *)(buffer + *next_offset - 4);

        GZIP_BLOCK block;
        block.index = index;
        block.offset = *offset;
        block.compressed_size = length;
        block.uncompressed_size = isize;
        block.compressed = (char *)(buffer + block.offset);
        block.uncompressed = NULL;
        dest->push_back(block);
    }

    /// Process last member.
    long last_offset = offsets.back();
    if (last_offset != dest->back().offset) {
        long length = file_size - last_offset;
        unsigned int isize = *(unsigned int *)(buffer + file_size - 4);

        GZIP_BLOCK block;
        block.index = index;
        block.offset = last_offset;
        block.compressed_size = length;
        block.uncompressed_size = isize;
        block.compressed = (char *)(buffer + block.offset);
        block.uncompressed = NULL;
        dest->push_back(block);
    }

    return true;
}

/**
 * @brief Thread function to decompress GZIP blocks.
 *
 * This function is executed by a thread to decompress a range of GZIP blocks
 * specified by the THREAD_PARAM argument. It allocates a decompressor, iterates
 * over the specified blocks, and decompresses each block. If the decompressed
 * size is insufficient, it reallocates the buffer with double the size and
 * retries decompression. The total decompressed size is accumulated in the
 * THREAD_PARAM structure.
 *
 * @param arg Pointer to a THREAD_PARAM structure containing the parameters for
 *            the thread, including the range of blocks to decompress and a
 *            pointer to the total size accumulator.
 * @return Always returns (void *)1.
 */
void *job_thread(void *arg) {
    THREAD_PARAM param = *(THREAD_PARAM *)arg;
    delete (THREAD_PARAM *)arg;
    libdeflate_decompressor *decompressor = libdeflate_alloc_decompressor();

    long i = 0;
    std::vector<GZIP_BLOCK>::iterator block;
    for (block = param.blocks->begin() + param.start, i = 0; i < param.length;
         ++block, ++i) {
        if (block->compressed == NULL)
            continue;
        if (block->uncompressed == NULL) {
            block->uncompressed = new char[block->uncompressed_size];
            //std::cout << "block id : " << i
            //          << ", size : " << block->uncompressed_size << std::endl;
        }
        size_t actual_in_nbytes = 0;
        size_t actual_out_nbytes = 0;
        enum libdeflate_result result;
        while (block->compressed_size != 0) {
            result = libdeflate_gzip_decompress_ex(
                decompressor, block->compressed, block->compressed_size,
                block->uncompressed, block->uncompressed_size,
                &actual_in_nbytes,
                &actual_out_nbytes); // &actual_in_bytes, &actual_out_bytes);
            if (result == LIBDEFLATE_INSUFFICIENT_SPACE) {
                block->uncompressed_size *= 2;
                delete[] block->uncompressed;
                block->uncompressed = new char[block->uncompressed_size];
                //std::cout << "block id : " << i << ", allocated_buffer_size : "
                //          << block->uncompressed_size << std::endl;
                continue;
            }
            block->compressed += actual_in_nbytes;
            block->compressed_size -= actual_in_nbytes;
            block->actual_out_nbytes += actual_out_nbytes;
        }
        *param.total_size += block->actual_out_nbytes;
        //std::cout << "block id : " << i
        //          << ", actual_out_nbytes : " << block->actual_out_nbytes
        //          << std::endl;
    }
    libdeflate_free_decompressor(decompressor);
    return (void *)1;
}

/**
 * @brief Thread function to merge decompressed GZIP blocks into a buffer.
 *
 * This function is intended to be run as a thread. It processes a range of GZIP blocks,
 * copying their uncompressed data into a specified buffer at the appropriate offsets.
 * After copying, it deallocates the memory used by the uncompressed data.
 *
 * @param arg Pointer to a MERGE_THREAD_PARAM structure containing the parameters for the thread.
 * @return Always returns (void *)1.
 *
 * The MERGE_THREAD_PARAM structure should contain:
 * - start: The starting index of the blocks to process.
 * - end: The ending index of the blocks to process.
 * - blocks: A pointer to a vector of GZIP_BLOCK structures.
 * - buffer: A pointer to the buffer where the uncompressed data should be copied.
 *
 * Each GZIP_BLOCK structure should contain:
 * - offset: The offset in the buffer where the uncompressed data should be copied.
 * - uncompressed: A pointer to the uncompressed data.
 * - uncompressed_size: The size of the uncompressed data.
 * - actual_out_nbytes: The actual number of bytes to copy from the uncompressed data.
 */
void *merge_thread(void *arg) {
    MERGE_THREAD_PARAM *param = (MERGE_THREAD_PARAM *)arg;
    long i = 0;
    for (i = param->start; i < param->end; ++i) {
        GZIP_BLOCK *block = &(*param->blocks)[i];
        if (block->offset < 0)
            continue;
        if (block->uncompressed == NULL)
            continue;
        if (block->uncompressed_size == 0)
            continue;
        //std::copy(block->uncompressed, block->uncompressed + block->actual_out_nbytes,
         //         param->buffer + block->offset);
        memcpy(param->buffer + block->offset, block->uncompressed,
              block->actual_out_nbytes);
        delete[] block->uncompressed;
    }
    return (void *)1;
}

void *buffer_release_job(void *arg) {
    char *buffer = (char *)arg;
    delete[] buffer;
    return NULL;
}

/**
 * @brief Unzips a given FASTQ file using multiple threads.
 *
 * This function reads a compressed FASTQ file, identifies GZIP blocks,
 * and decompresses them using multiple threads. The decompressed data
 * is then merged into a single buffer.
 *
 * @param file_name The name of the compressed FASTQ file.
 * @param thread_count The number of threads to use for decompression.
 * @return true if the file was successfully decompressed, false otherwise.
 */
bool FastqUnzip::Unzip(const std::string &file_name, int thread_count) {
    std::ifstream ifile(file_name.c_str(), std::ifstream::binary);
    ifile.seekg(0, ifile.end);
    long file_size = ifile.tellg();
    ifile.seekg(0, ifile.beg);

    unsigned char *compressed_buffer = new unsigned char[file_size];
    ifile.read((char *const)compressed_buffer,
               file_size); // 4L*1000*1000*1000);
    ifile.close();

    std::vector<GZIP_BLOCK> blocks;

    if (find_gzip_members((const char *const)compressed_buffer, file_size,
                          thread_count, &blocks) == false) {
        return false;
    }
    if (blocks.size() == 1) {
        blocks[0].uncompressed_size = blocks[0].compressed_size * 6;
    }
    int i = 0;
    long chunk_size = blocks.size() / thread_count;
    long start_index = 0;
    std::vector<pthread_t> threads;
    for (i = 0; i < thread_count; ++i) {
        THREAD_PARAM *p = new THREAD_PARAM;
        p->start = start_index;
        p->length = chunk_size;
        p->blocks = &blocks;
        p->total_size = &total_out_size;
        start_index += chunk_size;
        if (i < blocks.size() % thread_count) {
            p->length++;
            start_index++;
        }

        if (p->start + p->length > blocks.size())
            p->length = blocks.size() - p->start;

        pthread_t t;
        int s = pthread_create(&t, NULL, job_thread, (void *)p);
        threads.push_back(t);
    }
    std::vector<pthread_t>::iterator thread;
    for (thread = threads.begin(); thread != threads.end(); ++thread) {
        void *res = NULL;
        pthread_join(*thread, &res);
    }

    // release buffer
    pthread_t buffer_release_thread;
    pthread_create(&buffer_release_thread, NULL, buffer_release_job,
                   (void *)compressed_buffer);

    long offset = 0;
    std::vector<GZIP_BLOCK>::iterator block;
    for (block = blocks.begin(); block != blocks.end(); ++block) {
        block->offset = offset;
        offset += block->actual_out_nbytes;
    }
    if (blocks.size() > 1){
        uncompressed_buffer = new char[total_out_size];
        std::vector<MERGE_THREAD_PARAM *> merge_params;
        std::vector<pthread_t> merge_threads;
        chunk_size = blocks.size() / thread_count;
        start_index = 0;
        int end_index =0;
        for (i = 0; i < thread_count; ++i) {
            MERGE_THREAD_PARAM *p = new MERGE_THREAD_PARAM;
            p->buffer = uncompressed_buffer;
            p->start = end_index;
            p->end = start_index + chunk_size;
            end_index = p->end;
    
            if (i == thread_count - 1) {
                p->end = blocks.size();
            }

            if (p->start + p->end > blocks.size())
            {
                p->end = blocks.size() - p->start;
            }

            if(p->start == p->end)
                continue;

            p->blocks = &blocks;

            pthread_t t;
            int s = pthread_create(&t, NULL, merge_thread, (void *)p);
            if (s != 0) {
                std::cerr << "thread failed: [ " << block->index << " ]."
                        << std::endl;
                continue;
            }
            merge_params.push_back(p);
            merge_threads.push_back(t);
        }

        for (thread = merge_threads.begin(); thread != merge_threads.end();
            ++thread) {
            void *res = NULL;
            pthread_join(*thread, &res);
        }
    }
    else{
        uncompressed_buffer = blocks[0].uncompressed;
    }
    return true;
}
#endif // FASTQSTAT_UNZIP_H_
