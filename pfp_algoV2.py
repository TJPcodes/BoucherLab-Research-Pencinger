import logging
import os
from collections import defaultdict

# Setup for basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Special character types used for marking sequences in phrases
SPECIAL_TYPES = {
    "ENDOFDICT": '\x00',
    "ENDOFWORD": '\x01',
    "DOLLAR": '$',
    "DOLLAR_SEQUENCE": '#',
    "DOLLAR_PRIME": '&'
}

class MersenneKarpRabinHash:
    # Implements a Mersenne Karp-Rabin rolling hash for windowed sequence hashing
    def __init__(self, w):
        self.w = w  # Window size
        self.current_hash = 0  # Current hash value

    def initialize(self, data):
        # Initializes the hash with the first window of data
        self.current_hash = sum(map(ord, data))

    def update(self, prev, next):
        # Update the hash for the rolling window
        self.current_hash += ord(next) - ord(prev)

    def get_hash(self):
        # Return the current hash value
        return self.current_hash

    def reset(self):
        # Reset the hash for a new sequence
        self.current_hash = 0

class ParserFasta:
    # A parser class for FASTA files using custom methods for reading and custom hash for processing
    def __init__(self, params, out_file_prefix):
        self.params = params  # Parameters including window size
        self.out_file_prefix = out_file_prefix  # Output file prefix
        self.out_file_name = f"{out_file_prefix}.parse"  # Output file name
        self.out_file = None  # File handle for output file
        self.dictionary = {}  # Stores hashes of phrases
        self.trigger_strings = defaultdict(list)  # Stores trigger strings
        self.trigger_strings_list = []  # List to store trigger strings
        self.parse_size = 0  # Tracks the size of the parsed data
        self.closed = False  # Tracks if the file is closed

    def init(self):
        # Check for and delete the existing output file if it exists
        if os.path.exists(self.out_file_name):
            os.remove(self.out_file_name)
            logging.info(f"Existing output file {self.out_file_name} deleted.")
        
        # Open the output file and prepare for writing parsed data
        self.out_file = open(self.out_file_name, 'w')  # Open for text writing
        logging.info(f"Output file {self.out_file_name} opened for writing.")

    def process(self, in_file_path):
        # Process each sequence from a FASTA file using custom methods
        sequences = self.read_fasta(in_file_path)  # Read sequences from FASTA file
        kr_hash = MersenneKarpRabinHash(self.params['w'])  # Initialize the Karp-Rabin hash

        # Initialize the hash with a starting phrase
        phrase = [SPECIAL_TYPES["DOLLAR"]] * (self.params['w'] - 1) + [SPECIAL_TYPES["DOLLAR_SEQUENCE"]]
        kr_hash.initialize(phrase)

        # Iterate over each sequence in the FASTA file
        for sequence in sequences:
            print(sequence)
            print (phrase)
            # Process each character in the sequence
            for i, char in enumerate(sequence):
                if len(phrase) >= self.params['w']:
                    kr_hash.update(phrase[-self.params['w']], char)  # Update the hash for the new character
                phrase.append(char)

                if len(phrase) > self.params['w'] and kr_hash.get_hash() % self.params['p'] == 0:
                    hash_val = hash(tuple(phrase)) & 0xFFFFFFFF  # Get hash of the whole phrase
                    if hash_val not in self.dictionary:
                        self.dictionary[hash_val] = ''.join(phrase[:])  # Store the phrase in the dictionary
                    self.out_file.write(f"{hash_val}\n")  # Write hash as text to the output file
                    self.parse_size += 1

                    # Save the trigger string (last 2 characters) to the list
                    trigger_string = ''.join(phrase[-2:])
                    self.trigger_strings[hash_val].append(trigger_string)
                    self.trigger_strings_list.append(trigger_string)  # Add to the trigger strings list

                    # Reset phrase and KR Hash
                    phrase = phrase[-self.params['w']:]
                    kr_hash.reset()
                    kr_hash.initialize(phrase)

            # Handle the final window for the current sequence
            self.finalize_parsing(phrase)

        self.out_file.close()  # Close the output file

        # Sort dictionary by phrases and ensure no repeated phrases
        sorted_dict = sorted(set(self.dictionary.items()), key=lambda item: item[1])
        sorted_hashes = {hash_val: index for index, (hash_val, _) in enumerate(sorted_dict)}

        # Read the parse file and replace hash values with sorted indices
        with open(self.out_file_name, 'r') as file:
            parse_data = file.readlines()

        new_out_file_name = f"{self.out_file_prefix}.sorted.parse"
        if os.path.exists(new_out_file_name):
            os.remove(new_out_file_name)
            logging.info(f"Existing sorted output file {new_out_file_name} deleted.")
        
        with open(new_out_file_name, 'w') as new_file:
            for line in parse_data:
                hash_val = int(line.strip())
                new_file.write(f"{sorted_hashes[hash_val]}\n")  # Write sorted hash index as text

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
                self.dictionary[hash_val] = ''.join(phrase[:])
            self.out_file.write(f"{hash_val}\n")  # Write hash as text to the output file
            self.parse_size += 1

            # Save the trigger string (last 2 characters) to the list
            trigger_string = ''.join(phrase[-2:])
            self.trigger_strings[hash_val].append(trigger_string)
            self.trigger_strings_list.append(trigger_string)  # Add to the trigger strings list

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

    def save_dictionary_to_file(self, output_file_name):
        # Check for and delete the existing output file if it exists
        if os.path.exists(output_file_name):
            os.remove(output_file_name)
            logging.info(f"Existing output file {output_file_name} deleted.")
        
        # Save the dictionary of phrases to a text file
        with open(output_file_name, 'w') as file:
            sorted_phrases = sorted(set(self.dictionary.values()), key=lambda phrase: [phrase])
            for phrase in sorted_phrases:
                file.write(f"Phrase: {phrase}\n")

            # Save trigger strings
            file.write("\nTrigger Strings:\n")
            for trigger_string in self.trigger_strings_list:
                file.write(f"{trigger_string}\n")

    def save_trigger_strings_to_file(self, output_file_name):
        # Save the trigger strings to a text file
        with open(output_file_name, 'w') as file:
            for hash_val, positions in self.trigger_strings.items():
                file.write(f"Hash: {hash_val}, Trigger Strings: {positions}\n")

# Usage
params = {'w': 2, 'p': 5}  # Window size parameter
prefix = 'trial_1'  # Output file prefix
parser = ParserFasta(params, prefix)    
parser.init()

# Path to the input FASTA file
fasta_file_path = r"c:\Users\tyler\OneDrive\Desktop\VSCODE RESEARCH\sequences.fasta"
parser.process(fasta_file_path)

# Save the dictionary to a text file
output_file_name = "dictionary_output.txt"
parser.save_dictionary_to_file(output_file_name)

# Close the parser and perform cleanup
parser.close()



"""

Initialization:
The ParserFasta class is instantiated with specific parameters 
(w for window size) and an output file prefix. This setup prepares the parser
to write to an output file named according to the provided prefix, which will 
contain the parsed results. The MersenneKarpRabinHash class is also initialized
with a window size to handle the rolling hash calculations.

Reading and Processing:
The process method employs the read_fasta function to read sequences
from a specified FASTA file. This method reads the file line by line,
collecting sequences. Each sequence retrieved is then processed individually.
A rolling hash mechanism (Mersenne Karp-Rabin) is applied to each sequence
as it is read character by character, updating the hash continuously.
When the phrase length exceeds the window size (w), the hash value of 
the phrase is calculated and stored in the dictionary if it is not already present. 
The hash value is written to the output file as text, and trigger string positions are logged.

Phrase Resetting:
After every phrase is added to the dictionary, 
the phrase variable is reset to contain only the last w characters of 
the just added phrase. Additionally, the Karp-Rabin hash is reset and 
reinitialized with this updated phrase. This ensures that the rolling 
window correctly continues with the next part of the sequence.

Finalizing Parsing:
After all characters in a sequence are processed,
finalize_parsing manages the final window of characters 
by appending special end symbols and writing the final hash value to the output file. 
This method ensures that the parsing results are properly recorded, 
including the last segment of the sequence. If the phrase length 
and content meet specific conditions, the phrase is appended with 
special symbols, hashed, and written to the file. The phrase is then reset, 
and the Karp-Rabin hash is reinitialized to handle any remaining characters correctly.

Sorting and Replacing:
The dictionary is sorted based on the phrases, not the keys. 
Each hash in the parse file is then replaced with the position of the corresponding 
hash in the newly sorted dictionary. For example, if a hash 213 in the parse 
corresponds to the fifth dictionary item in the sorted dictionary, 213 is 
replaced with 4 (using 0-based indexing). The sorted indices are written to a new output file as text.

Cleanup:
The close method is called to ensure the output file is properly 
closed after all parsing activities are completed, marking the end
of file processing and securing the written data.

Static Methods:

read_fasta: Reads sequences directly from a FASTA file and returns 
them as a list, which simplifies the handling of sequence data.
Usage:
The ParserFasta class is initialized with the given parameters and 
output file prefix. The process method processes the input FASTA file, 
and the save_dictionary_to_file method saves the dictionary of phrases 
to a text file. Finally, the close method ensures that the output file 
is properly closed and the parser is cleaned up.

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