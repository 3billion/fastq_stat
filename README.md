# Fastq Stat

Fastq Stat is a high-performance tool designed to process and analyze FASTQ files, which are commonly used in bioinformatics to store nucleotide sequences and their corresponding quality scores.

## Prerequisites

Before building and running the Fastq Stat program, ensure you have the following prerequisites installed:

- **Intel oneAPI DPC++/C++ Compiler**: Required for compiling the program. You can download it from the [Intel oneAPI Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html).
- **libdeflate library**: Used for efficient decompression of gzip-compressed FASTQ files. You can clone and build it from the [libdeflate GitHub repository](https://github.com/ebiggers/libdeflate).
- **large size memory**: It loads two fastq.gz files into memory and decompresses them. So, it loads two fastq.gz and two fastq files into memory at the same time.
- **avx512bw**: cpu that supports avx512bw


## Architecture

The program is written in C++ and leverages the Intel oneAPI DPC++/C++ Compiler for optimal performance. It uses the libdeflate library for efficient decompression of gzip-compressed FASTQ files. The program is designed to take advantage of modern CPU features, such as AVX-512, to accelerate computation.

## Build

To build the Fastq Stat program, use the following command:

```sh
/opt/intel/oneapi/compiler/latest/bin/icpx -o fastq_stat fastq_stat.cpp ./libdeflate/build/libdeflate.a -lpthread -std=c++20 -O3 -mavx512bw
```

This command compiles the `fastq_stat.cpp` source file, links it with the libdeflate library, and optimizes the code for AVX-512 capable processors.

## Execute

To execute the Fastq Stat program, use the following command:

```sh
./fastq_stat {input1.fastq.gz} {input2.fastq.gz} {thread_count}
```

- `{input1.fastq.gz}`: Path to the first input FASTQ file.
- `{input2.fastq.gz}`: Path to the second input FASTQ file.
- `{thread_count}`: Number of threads to use for processing.

## Example

Here is an example of how to run the program:

```sh
./fastq_stat ../EPG24-AAMY_1.fastq.gz ../EPG24-AAMY_2.fastq.gz 10
```

This command processes the `EPG24-AAMY_1.fastq.gz` and `EPG24-AAMY_2.fastq.gz` files using 10 threads.