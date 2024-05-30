from pathlib import WindowsPath
import struct
import Bio 
from Bio import SeqIO
import logging
from collections import defaultdict

# Setup for basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MersenneKarpRabinHash:
    #Implements a Mersenne Karp-Rabin rolling hash for windowed sequence hashing
    def __init__(self, w):
        self.w = w
        self.current_hash = 0

    def initialize(self, data):
        #Initializes the hash with the first window of data
        self.current_hash = sum(data)

    def update(self, prev, next):
        #Update the hash for the rolling window
        self.current_hash += next - prev

    def get_hash(self):
        #Return the current hash value
        return self.current_hash

    def reset(self):
        #Reset the hash for a new sequence
        self.current_hash = 0

class ParserFasta:
    #A parser class for FASTA files using BioPython for reading and custom hash for processing
    def __init__(self, params, out_file_prefix):
        self.params = params
        self.out_file_prefix = out_file_prefix
        self.out_file_name = f"{out_file_prefix}.parse"
        self.out_file = None
        self.dictionary = defaultdict(int)  # Stores hashes of phrases and their frequencies
        self.trigger_strings = {}  # Stores positions of trigger strings
        self.parse_size = 0
        self.closed = False

    def init(self):
        #Open the output file and prepare for writing parsed data
        self.out_file = open(self.out_file_name, 'wb')
        logging.info(f"Output file {self.out_file_name} opened for writing.")

    def process(self, in_file_path):
        #Process each sequence from a FASTA file using BioPython's SeqIO and custom methods
        sequences = self.read_fasta(in_file_path)  # Using the new integrated read_fasta method
        for record in sequences:
            self.process_record(str(record.seq))

    def process_record(self, sequence):
        #Process a single sequence from the FASTA file, applying a rolling hash and storing results
        kr_hash = MersenneKarpRabinHash(self.params['w'])
        phrase = [ord('$')]  # Start each new sequence with a special start symbol

        for i, char in enumerate(sequence):
            char = ord(char)
            if len(phrase) >= self.params['w']:
                kr_hash.update(phrase[-self.params['w']-1], char)
            phrase.append(char)

            if len(phrase) > self.params['w']:
                hash_val = kr_hash.get_hash()
                if hash_val % self.params['p'] == 0:
                    window_tuple = tuple(phrase[-self.params['w']:])
                    self.dictionary[WindowsPath(tuple), ] += 1
                    # Log position of each trigger string
                    self.trigger_strings.setdefault(window_tuple, []).append(i - self.params['w'] + 1)

        # After processing all characters, handle the final window
        self.finalize_parsing(phrase)

    def finalize_parsing(self, phrase):
        #Handle the last phrase after processing all characters
        if len(phrase) >= self.params['w']:
            phrase.extend([ord('$')] * self.params['w'])  # Append ending symbol
            final_hash = self.dictionary[tuple(phrase[-self.params['w']:])]
            self.out_file.write(struct.pack('I', final_hash))

    def close(self):
        #Close the output file and perform cleanup
        if not self.closed and self.out_file:
            self.out_file.close()
            self.closed = True
            logging.info("Output file closed and parser cleanup completed.")

    @staticmethod
    def read_fasta(file_path):
        #Read a FASTA file and return a list of sequences
        sequences = []
        with open(file_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append(record)
        return sequences

    @staticmethod
    def parse_fasta(file_path):
        #Parse a FASTA file and return a dictionary of sequences indexed by their IDs
        sequences = {}
        with open(file_path, 'r') as file:
            current_id = None
            current_seq = []
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    current_id = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
        return sequences

# Usage
params = {'w': 5, 'p': 1}
prefix = 'example_output'
parser = ParserFasta(params, prefix)
parser.init()
parser.process("path/to/your/fasta_file.fasta")  # Directs processing from a FASTA file path
parser.close()

'''
SUMMARY:

Initialization:

The ParserFasta class is instantiated with specific parameters 
(w for window size, and p for modulus value) and an output file prefix.
This setup prepares the parser to write to an output file named according to the provided prefix, 
which will contain the parsed results.

Reading and Processing:

The process method employs the read_fasta function to read sequences from a specified FASTA file. 
This method leverages BioPython's SeqIO.parse to load and iterate through sequences in the file efficiently. 
Each sequence retrieved is then processed individually.

Sequence Processing:

Within the process_record method, each sequence undergoes processing where a rolling hash mechanism 
(Mersenne Karp-Rabin) is applied. As the sequence is read character by character, the hash is continuously updated.
When the hash value of the current window (of length w) modulo p equals zero (trigger condition), 
the window is converted into a tuple and recorded in a dictionary that counts the occurrences of 
each unique window, or "phrase".

Finalizing Parsing:

After all characters in a sequence are processed, finalize_parsing manages the final window 
of characters by appending a special end symbol ($) and writing the final hash value to the output file. 
This method ensures that the parsing results are properly recorded, including the last segment of the sequence.

Cleanup:

The close method is called to ensure the output file is properly closed after all parsing 
activities are completed, marking the end of file processing and securing the written data.

Static Methods:

read_fasta: Reads sequences directly from a FASTA file and returns them as a list, 
which simplifies the handling of sequence data.
parse_fasta: Parses a FASTA file and returns a dictionary of sequences indexed by their headers (sequence IDs). 
This method provides quick access to specific sequences, facilitating efficient data manipulation and analysis.

'''


'''
SOURCES:

https://biopython.org/DIST/docs/tutorial/Tutorial-1.83.pdf
https://github.com/marco-oliva/pfp/blob/master/pfp%2B%2B.cpp
https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
chat.openai.com for some C++ conversions
https://brilliant.org/wiki/rabin-karp-algorithm/#:~:text=The%20Rabin%2DKarp%20algorithm%20is,important%20application%20of%20computer%20science.
https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/


'''