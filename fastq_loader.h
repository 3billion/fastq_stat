#ifndef FASTQSTAT_LOADER_H_
#define FASTQSTAT_LOADER_H_
#include "fastq_unzip.h"
/**
 * @class FastqLoader
 * @brief A class to load and process FASTQ files.
 *
 * The FastqLoader class provides functionality to load FASTQ files
 * and manage the associated resources.
 */
class FastqLoader {
  public:
    FastqLoader() {};
    ~FastqLoader() {
        if (decomper != nullptr) {
            delete decomper;
        }
        if (buffer != nullptr) {
            delete[] buffer;
        }
    };
    bool Load(const std::string &file_name, int thread_count);
    char *buffer = nullptr;
    long buffer_size = 0;

  private:
    FastqUnzip *decomper = nullptr;
};

/**
 * @brief Loads a FASTQ file, either compressed (.gz) or uncompressed.
 *
 * This function attempts to load a FASTQ file specified by the file_name parameter.
 * If the file is compressed (with a .gz extension), it will be decompressed using
 * the FastqUnzip class. If the file is uncompressed, it will be read directly.
 *
 * @param file_name The name of the FASTQ file to load.
 * @param thread_count The number of threads to use for decompression (if applicable).
 * @return true if the file was successfully loaded, false otherwise.
 */
bool FastqLoader::Load(const std::string &file_name, int thread_count) {
    std::ifstream file;

    if (file_name.find(".gz") != std::string::npos) {
        decomper = new FastqUnzip();
        if (!decomper->Unzip(file_name, thread_count)) {
            std::cerr << "Failed to decompress file: " << file_name
                      << std::endl;
            return false;
        }
        buffer = decomper->uncompressed_buffer;
        buffer_size = decomper->total_out_size;
    } else {
        file.open(file_name);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << file_name << std::endl;
            return false;
        }

        // Get file size
        file.seekg(0, std::ios::end);
        buffer_size = file.tellg();
        file.seekg(0, std::ios::beg);

        // Allocate buffer
        buffer = new char[buffer_size];
        if (buffer == nullptr) {
            std::cerr << "Failed to allocate buffer" << std::endl;
            return false;
        }
        // Read file
        file.read(buffer, buffer_size);
        file.close();
    }
    //std::cout << "load finish. file_name : " << file_name
    //          << ", load_size : " << buffer_size << std::endl;
    return true;
}
#endif // FASTQSTAT_LOADER_H_