import struct
import logging
from collections import defaultdict

# Setup for basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

SPECIAL_TYPES = {
    "ENDOFDICT": 0,
    "ENDOFWORD": 1,
    "DOLLAR": 2,
    "DOLLAR_SEQUENCE": 4,
    "DOLLAR_PRIME": 5
}

class MersenneKarpRabinHash:
    # Implements a Mersenne Karp-Rabin rolling hash for windowed sequence hashing
    def __init__(self, w):
        self.w = w
        self.current_hash = 0

    def initialize(self, data):
        # Initializes the hash with the first window of data
        self.current_hash = sum(data)

    def update(self, prev, next):
        # Update the hash for the rolling window
        self.current_hash += next - prev

    def get_hash(self):
        # Return the current hash value
        return self.current_hash

    def reset(self):
        # Reset the hash for a new sequence
        self.current_hash = 0

class ParserFasta:
    # A parser class for FASTA files using custom methods for reading and custom hash for processing
    def __init__(self, params, out_file_prefix):
        self.params = params
        self.out_file_prefix = out_file_prefix
        self.out_file_name = f"{out_file_prefix}.parse"
        self.out_file = None
        self.dictionary = {}  # Stores hashes of phrases
        self.trigger_strings = defaultdict(list)  # Stores positions of trigger strings
        self.parse_size = 0
        self.closed = False

    def init(self):
        # Open the output file and prepare for writing parsed data
        self.out_file = open(self.out_file_name, 'wb')
        logging.info(f"Output file {self.out_file_name} opened for writing.")

    def process(self, in_file_path):
        # Process each sequence from a FASTA file using custom methods
        sequences = self.read_fasta(in_file_path)
        kr_hash = MersenneKarpRabinHash(self.params['w'])

        # Initialize the hash with a starting phrase
        phrase = [SPECIAL_TYPES["DOLLAR"]] * (self.params['w'] - 1) + [SPECIAL_TYPES["DOLLAR_SEQUENCE"]]
        kr_hash.initialize(phrase)

        for sequence in sequences:
            for i, char in enumerate(sequence):
                char = ord(char)
                if len(phrase) >= self.params['w']:
                    kr_hash.update(phrase[-self.params['w']], char)
                phrase.append(char)

                if len(phrase) > self.params['w']:
                    hash_val = hash(tuple(phrase)) & 0xFFFFFFFF # Get hash of the whole phrase
                    if hash_val not in self.dictionary:
                        self.dictionary[hash_val] = phrase[:]
                    self.out_file.write(struct.pack('I', hash_val))
                    self.parse_size += 1

                    # Log position of each trigger string
                    self.trigger_strings[hash_val].append(i - self.params['w'] + 1)

                    # Reset phrase and KR Hash
                    phrase = phrase[-self.params['w']:]
                    kr_hash.reset()
                    kr_hash.initialize(phrase)

            # Handle the final window for the current sequence
            self.finalize_parsing(phrase)

        self.out_file.close()

        # Sort the dictionary by phrases
        sorted_dict = sorted(self.dictionary.items(), key=lambda item: item[1])
        sorted_hashes = {hash_val: index for index, (hash_val, _) in enumerate(sorted_dict)}

         # Read the file, replace hashes with positions, and write to a new file
        with open(self.out_file_name, 'rb') as file:
            parse_data = file.read()

        new_out_file_name = f"{self.out_file_prefix}.sorted.parse"
        with open(new_out_file_name, 'wb') as new_file:
            for i in range(0, len(parse_data), 4):
                hash_val, = struct.unpack('I', parse_data[i:i+4])
                new_file.write(struct.pack('I', sorted_hashes[hash_val]))

        logging.info(f"Output file sorted and saved as {new_out_file_name}")

    def finalize_parsing(self, phrase):
    # Handles the last phrase after processing all characters
        if phrase[0] != SPECIAL_TYPES["DOLLAR"] and len(phrase) >= self.params['w']:
            # Append w-1 dollar prime, and one dollar seq at the end of each sequence
            for _ in range(self.params['w'] - 1):
                phrase.append(SPECIAL_TYPES["DOLLAR_PRIME"])
            phrase.append(SPECIAL_TYPES["DOLLAR_SEQUENCE"])

            # Get hash of the final phrase and add to the dictionary if it's not present
            hash_val = hash(tuple(phrase)) & 0xFFFFFFFF
            if hash_val not in self.dictionary:
                self.dictionary[hash_val] = phrase[:]
            self.out_file.write(struct.pack('I', hash_val))
            self.parse_size += 1

            # Resets phrase
            phrase.clear()
            for _ in range(self.params['w'] - 1):
                phrase.append(SPECIAL_TYPES["DOLLAR_PRIME"])
            phrase.append(SPECIAL_TYPES["DOLLAR_SEQUENCE"])
            
            # Resets and reinitializes KR hash
            kr_hash = MersenneKarpRabinHash(self.params['w'])
            kr_hash.reset()
            kr_hash.initialize(phrase)

    def close(self):
        # Close the output file and perform cleanup
        if not self.closed and self.out_file:
            self.out_file.close()
            self.closed = True
            logging.info("Output file closed and parser cleanup completed.")

    @staticmethod
    def read_fasta(file_path):
        # Reads a FASTA file and return a list of sequences
        sequences = []
        with open(file_path, 'r') as file:
            current_seq = []
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(''.join(current_seq))
                        current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                sequences.append(''.join(current_seq))
        return sequences

# Usage
params = {'w': 10, 'p': 3}
prefix = 'example_output'
parser = ParserFasta(params, prefix)    
parser.init()

fasta_file_path = r"c:\Users\tyler\OneDrive\Desktop\VSCODE RESEARCH\sequences.fasta"
parser.process(fasta_file_path)
parser.close()  

print(f"Trigger strings: {parser.trigger_strings}")

"""

Initialization:

The ParserFasta class is instantiated with specific parameters
(w for window size, and p for modulus value) and an output file prefix.
This setup prepares the parser to write to an output file named according to the provided prefix,
which will contain the parsed results.

Reading and Processing:

The process method employs the read_fasta function to read sequences from a specified FASTA file.
This method reads the file line by line, collecting sequences. Each sequence retrieved is
then processed individually.

Sequence Processing:

Within the process method, each sequence undergoes processing where a
rolling hash mechanism (Mersenne Karp-Rabin) is applied. As the sequence is read character by character,
the hash is continuously updated. The hashing now uses the whole phrase
instead of just the last w characters. When the hash value of the current
window (of length w) modulo p equals zero (trigger condition), the window is converted
into a tuple and hashed. This hash is used as a key to store the phrase in a dictionary,
and the hash is written to the output file. Trigger string information is logged.

Phrase Resetting:

After every phrase is added to the dictionary, the phrase variable is reset to contain only the last
w characters of the just added phrase. Additionally, the Karp-Rabin hash is reset and reinitialized
with this updated phrase. This ensures that the rolling window correctly continues with the next part
of the sequence.

Finalizing Parsing:

After all characters in a sequence are processed, finalize_parsing manages the final window of characters
by appending special end symbols and writing the final hash value to the output file.
This method ensures that the parsing results are properly recorded, including the last segment of the sequence.
If the phrase length and content meet specific conditions, the phrase is appended with special symbols,
hashed, and written to the file. The phrase is then reset, and the Karp-Rabin hash is reinitialized
to handle any remaining characters correctly.

Sorting and Replacing:

The dictionary is sorted based on the phrases, not the keys. Each hash in the parse is then replaced
with the position of the corresponding hash in the newly sorted dictionary. For example, if a hash 213
in the parse corresponds to the fifth dictionary item in the sorted dictionary, 213 is replaced with 4
(using 0-based indexing).

Cleanup:

The close method is called to ensure the output file is properly closed after all parsing activities are completed,
marking the end of file processing and securing the written data.

Static Methods:

read_fasta: Reads sequences directly from a FASTA file and returns them as a list,
which simplifies the handling of sequence data.

"""


'''
SOURCES:

https://biopython.org/DIST/docs/tutorial/Tutorial-1.83.pdf
https://github.com/marco-oliva/pfp/blob/master/pfp%2B%2B.cpp
https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
chat.openai.com for some C++ conversions
https://brilliant.org/wiki/rabin-karp-algorithm/#:~:text=The%20Rabin%2DKarp%20algorithm%20is,important%20application%20of%20computer%20science.
https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/


'''