#include <iostream>
#include <vector>
#include <fstream>
#include <string>

// Function to read a FASTA file and return a vector of sequences
std::vector<std::string> readFASTA(const std::string& file_path) {
    std::ifstream file(file_path);
    std::vector<std::string> sequences;
    std::string line, sequence;

    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }

    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }

    return sequences;
}

// Kernel function to compute the rolling hash values for the input string `s`.
__global__ void computeRollingHash(const char* s, int window_size, int* hash_values, int base, int mod, int* power, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx <= n - window_size) {
        int current_hash = 0;

        // Compute the hash value for the initial window starting at `idx`.
        for (int i = 0; i < window_size; ++i) {
            current_hash = (current_hash * base + s[idx + i]) % mod;
        }

        hash_values[idx] = current_hash;

        // Update the hash value for the subsequent windows
        for (int i = 1; i < n - window_size + 1 - idx; ++i) {
            current_hash = (current_hash - s[idx + i - 1] * power[window_size - 1]) % mod;
            if (current_hash < 0) current_hash += mod;
            current_hash = (current_hash * base + s[idx + i + window_size - 1]) % mod;
            hash_values[idx + i] = current_hash;
        }
    }
}

int main() {
    // Read sequences from FASTA file
    std::vector<std::string> sequences = readFASTA("/blue/boucher/tyler.pencinger/sequences.fasta");

    // Set up CUDA processing for each sequence
    for (const std::string& s : sequences) {
        int window_size = 2; // As in the example with "world"
        int base = 31; // Prime number as base
        int mod = 1000000007; // Large prime number
        int n = s.length();

        std::vector<int> hash_values(n - window_size + 1);
        std::vector<int> power(n + 1, 1);

        // Precompute the powers of the base modulo the mod
        for (int i = 1; i <= n; i++) {
            power[i] = (power[i - 1] * base) % mod;
        }

        // Device pointers
        char* d_s;
        int* d_hash_values;
        int* d_power;

        // Allocate memory on GPU
        cudaMalloc(&d_s, n * sizeof(char));
        cudaMalloc(&d_hash_values, (n - window_size + 1) * sizeof(int));
        cudaMalloc(&d_power, (n + 1) * sizeof(int));

        // Copy data to GPU
        cudaMemcpy(d_s, s.c_str(), n * sizeof(char), cudaMemcpyHostToDevice);
        cudaMemcpy(d_power, power.data(), (n + 1) * sizeof(int), cudaMemcpyHostToDevice);

        int threadsPerBlock = 256;
        int blocksPerGrid = (n - window_size + 1 + threadsPerBlock - 1) / threadsPerBlock;

        // Launch kernel
        computeRollingHash<<<blocksPerGrid, threadsPerBlock>>>(d_s, window_size, d_hash_values, base, mod, d_power, n);

        cudaDeviceSynchronize();

        // Copy results back to host
        cudaMemcpy(hash_values.data(), d_hash_values, (n - window_size + 1) * sizeof(int), cudaMemcpyDeviceToHost);

        // Free GPU memory
        cudaFree(d_s);
        cudaFree(d_hash_values);
        cudaFree(d_power);

        // Print the hash values
        std::cout << "Hash values for sequence: " << s << std::endl;
        for (int hash : hash_values) {
            std::cout << hash << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}