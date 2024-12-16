#include "fastq_calculator.h"
#include "fastq_loader.h"
#include <atomic>
#include <iostream>
#include <thread>


void thread_job(FastqLoader &loader, const std::string &file_name,
                int thread_count, Stat &stat_file) {
    if (loader.Load(file_name, thread_count) == false) {
        std::cerr << "Failed to load file: " << file_name << std::endl;
        return;
    }
    FastqCalculator calculator;

    if (calculator.Calculate(loader.buffer, loader.buffer_size,
                                thread_count, stat_file)) {
        calculator.Free();
    }
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0]
                  << " <fastq_file1> <fastq_file2> <thread_count>" << std::endl;
        return -1;
    }
    const int file_count = 2;
    std::string file_name[file_count];
    std::string thread_count_str(argv[3]);
    int thread_count = std::stoi(thread_count_str);
    int thread_count_per_file = (thread_count + 0.5) / file_count;
    std::thread *threads[file_count];
    FastqLoader loader[file_count];
    Stat stat_file[file_count];
    Stat stat_total;

    for(int i = 0; i<file_count; i++)
    {
        file_name[i] = std::string(argv[i + 1]);
        threads[i] = new std::thread(thread_job, std::ref(loader[i]), file_name[i], thread_count_per_file, std::ref(stat_file[i]));
    }
    for (int i = 0; i < file_count; i++) {
        threads[i]->join();
        stat_total.linenum += stat_file[i].linenum;
        stat_total.yield += stat_file[i].yield;
        stat_total.q20_base += stat_file[i].q20_base;
        stat_total.q30_base += stat_file[i].q30_base;
    }

    double read_length = (double)stat_total.yield / (double)stat_total.linenum;
    double Q20 = (double)stat_total.q20_base / (double)stat_total.yield * 100;
    double Q30 = (double)stat_total.q30_base / (double)stat_total.yield * 100;
    std::cout << "{\n";
    std::cout << "  \"Total Yield\": " << stat_total.yield << ",\n";
    std::cout << "  \"Total reads\": " << stat_total.linenum << ",\n";
    std::cout << "  \"Average read length\": " << read_length << ",\n";
    std::cout << "  \"Q20(%)\": " << Q20 << ",\n";
    std::cout << "  \"Q30(%)\": " << Q30 << "\n";
    std::cout << "}" << std::endl;

    return 0;
}
