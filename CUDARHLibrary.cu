#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cuda_runtime.h>

// Function to read a FASTA file and return a vector of sequences
std::vector<std::string> readFASTA(const std::string& file_path) {
    std::ifstream file(file_path);
    std::vector<std::string> sequences;
    std::string line, sequence;

    while (std::getline(file, line)) {
        if (line[0] == '>') {  // Header line
            if (!sequence.empty()) {
                sequences.push_back(sequence);  // Save the current sequence
                sequence.clear();
            }
        } else {
            sequence += line;  // Append sequence lines
        }
    }

    if (!sequence.empty()) {
        sequences.push_back(sequence);  // Add the last sequence
    }

    return sequences;
}

// Kernel function to compute the rolling hash values for the input string `s`.
__global__ void computeRollingHash(const char* s, int window_size, int* hash_values, int base, int mod, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx <= n - window_size) {
        int current_hash = 0;
        int base_power = 1;  // Keeps track of base^(window_size-1) % mod

        // Compute the hash value for the initial window starting at `idx`.
        for (int i = 0; i < window_size; ++i) {
            current_hash = (current_hash * base + s[idx + i]) % mod;
            if (i < window_size - 1) {
                base_power = (base_power * base) % mod;
            }
        }

        hash_values[idx] = current_hash;

        // Update the hash value for the subsequent windows
        for (int i = 1; i < n - window_size + 1 - idx; ++i) {
            current_hash = (current_hash * base - s[idx + i - 1] * base_power % mod + mod) % mod;
            current_hash = (current_hash + s[idx + i + window_size - 1]) % mod;
            hash_values[idx + i] = current_hash;
        }
    }
}

// Expose the main functionality to Python
extern "C" void processFASTA(const char* file_path, int window_size, int base, int mod) {
    // Read sequences from FASTA file
    std::vector<std::string> sequences = readFASTA(file_path);

    // Set up CUDA processing for each sequence
    for (const std::string& s : sequences) {
        int n = s.length();
        std::vector<int> hash_values(n - window_size + 1);

        // Device pointers
        char* d_s;
        int* d_hash_values;

        // Allocate memory on GPU
        cudaMalloc(& d_s, n * sizeof(char));
        cudaMalloc(& d_hash_values, (n - window_size + 1) * sizeof(int));

        // Copy data to GPU
        cudaMemcpy(d_s, s.c_str(), n * sizeof(char), cudaMemcpyHostToDevice);

        int threadsPerBlock = 256;
        int blocksPerGrid = (n - window_size + 1 + threadsPerBlock - 1) / threadsPerBlock;

        // Launch kernel
        computeRollingHash <<< blocksPerGrid, threadsPerBlock >>> (d_s, window_size, d_hash_values, base, mod, n);

        // Check for any kernel errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl;
        }

        cudaDeviceSynchronize();

        // Copy results back to host
        cudaMemcpy(hash_values.data(), d_hash_values, (n - window_size + 1) * sizeof(int), cudaMemcpyDeviceToHost);

        // Free GPU memory
        cudaFree(d_s);
        cudaFree(d_hash_values);

        // Print the sequence and hash values
        std::cout << "Hash values for sequence: " << s << std::endl;
        for (int hash : hash_values) {
            std::cout << hash << " ";
        }
        std::cout << std::endl;
    }
}