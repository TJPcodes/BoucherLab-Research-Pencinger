
#include <pfp_algo.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <spdlog/spdlog.h>

from inspect import _void


def vcfbwt (
    class pfp (
        class Mersenne_KarpRabinHash ()
        public:
            Mersenne_KarpRabinHash(int w) {}
            void initialize(const char* data, int w) {}
            void update(char prev, char next) {}
            int get_hash() { return 0; }
            void reset() {}
    );

        class ReferenceParse {
        public:
            void init(const std::string& reference) {}
            std::vector<int> parse;
            std::vector<int> trigger_strings_position;
            std::vector<char> dictionary;
            std::vector<char> acgt_only_table;
            struct Params {
                int w;
                bool acgt_only;
                int p;
            } params;
        };

        class ParserVCF {
        public:
            void init(const ReferenceParse::Params& params, const std::string& prefix, ReferenceParse& rp, int t) {}
            void operator()(const Sample& sample) {}
            void close() {}
            std::vector<int> parse;
            std::vector<int> trigger_strings_position;
            std::vector<int> samples_processed;
            std::vector<int> registered_workers;
            std::vector<char> dictionary;
            std::vector<char> acgt_only_table;
            std::string out_file_prefix;
            std::string out_file_name;
            std::string tmp_out_file_name;
            std::ofstream out_file;
            ReferenceParse* reference_parse;
            ReferenceParse::Params params;
            int w;
            int p;
            int tags;
            int parse_size;
            bool closed;
            char DOLLAR;
            char DOLLAR_PRIME;
            char DOLLAR_SEQUENCE;
            char ENDOFWORD;
            char ENDOFDICT;
            char MAIN;
            char WORKER;
            char COMPRESSED;
            char EXT;
            char EXIT_FAILURE;
            char acgt_only;
            char in_file_path;
            char working_genotype;
        };

        class ParserFasta {
        public:
            void init(const ReferenceParse::Params& params, const std::string& prefix) {}
            void operator()() {}
            void close() {}
            std::vector<int> parse;
            std::vector<int> trigger_strings_position;
            std::vector<int> sequences_processed;
            std::vector<char> dictionary;
            std::vector<char> acgt_only_table;
            std::string out_file_prefix;
            std::string out_file_name;
            std::string tmp_out_file_name;
            std::ofstream out_file;
            ReferenceParse::Params params;
            int w;
            int p;
            int parse_size;
            bool closed;
            char DOLLAR;
            char DOLLAR_PRIME;
            char DOLLAR_SEQUENCE;
            char ENDOFWORD;
            char ENDOFDICT;
            char acgt_only;
            char in_file_path;
        };
    )

_void vcfbwt::pfp::ReferenceParse::init(const std::string& reference) {
    std::vector<char> phrase;
    spdlog::info("Parsing reference");

    #Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);

    #Reference as first sample, just one dollar to be compatible with Giovanni's pscan.cpp
    phrase.emplace_back(DOLLAR);

    for (std::size_t ref_it = 0; ref_it < reference.size(); ref_it++) {
        char c = reference[ref_it];
        if (params.acgt_only) { c = acgt_only_table[c]; }

        phrase.push_back(c);
        if (phrase.size() == params.w) { kr_hash.initialize(phrase.data(), params.w); }
        else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }

        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0)) {
            int hash = this->dictionary.check_and_add(phrase);

            this->parse.push_back(hash);
            this->trigger_strings_position.push_back(ref_it - this->params.w + 1);

            phrase.erase(phrase.begin(), phrase.end() - this->params.w); # Keep the last w chars

            kr_hash.reset(); kr_hash.initialize(phrase.data(), params.w);
        }

        #Review 'std' and 'this' and look into making faster via python shortcuts
        '''def is_valid_dna_sequence(sequence):
        valid_chars = {'A', 'C', 'G', 'T'}
    for char in sequence:
        if char not in valid_chars:
            return False
    return True

# Example usage
sequence = "ACGTAGTCAGT"
print(is_valid_dna_sequence(sequence))  # Output: True

invalid_sequence = "ACGTBXG"
print(is_valid_dna_sequence(invalid_sequence))  # Output: False'''

    }


    # Last phrase
    if (phrase.size() >= this->params.w) {
        // Append w-1 dollar prime and 1 dollar sequence at the end, reference as the first sample
        for (std::size_t j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
        phrase.emplace_back(DOLLAR_SEQUENCE);

        int hash = this->dictionary.check_and_add(phrase);

        this->parse.push_back(hash);
        this->trigger_strings_position.push_back(reference.size() - 1);
    }
    else { spdlog::error("The reference doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}

void vcfbwt::pfp::ParserVCF::init(const ReferenceParse::Params& params, const std::string& prefix, ReferenceParse& rp, int t) {
    this->w = params.w; this->out_file_prefix = prefix; this->p = params.p; this->tags = t; this->parse_size = 0;
    if (not ((tags & MAIN) or (tags & WORKER))) { spdlog::error("A parser must be either the main parser or a worker"); std::exit(EXIT_FAILURE); }

    if (tags & MAIN) {this->out_file_name = out_file_prefix + EXT::PARSE; }
    if ((tags & MAIN) and params.compress_dictionary) { tags = tags | COMPRESSED; }
    this->tmp_out_file_name = TempFile::getName("parse");
    this->out_file.open(tmp_out_file_name, std::ios::binary);
    this->reference_parse = &rp;
    this->dictionary = &this->reference_parse->dictionary;

    this->params = params;
}

void vcfbwt::pfp::ParserVCF::operator()(const Sample& sample) {
    Sample::iterator sample_iterator(sample, this->working_genotype);
    this->samples_processed.push_back(sample.id());

    std::vector<char> phrase;

    // Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);

    // Every sample starts with w-1 dollar prime and one dollar seq
    for (std::size_t j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
    phrase.emplace_back(DOLLAR_SEQUENCE);
    kr_hash.initialize(phrase.data(), params.w);

    // Shorthands
    std::vector<int>& tsp = reference_parse->trigger_strings_position;

    std::size_t start_window = 0, end_window = 0;
    while (not sample_iterator.end()) {
        // Compute where we are on the reference
        std::size_t pos_on_reference = sample_iterator.get_ref_it();

        if ( not ((sample_iterator.get_var_it() > 0) and (sample_iterator.prev_variation_end() > (pos_on_reference - (2 * this->w))))) {
            // Set start postion to the position in the reference parse after the last computed phrase
            if (params.use_acceleration and ((phrase.size() == this->w) and ((pos_on_reference != 0) and (phrase[0] != DOLLAR_PRIME)))){
                start_window = end_window;
                while ((tsp[start_window] + this->w) <= pos_on_reference and (start_window < tsp.size() - 2))
                { start_window++; }

                // Iterate over the parse up to the next variation
                while (tsp[end_window + 1] < (sample_iterator.next_variation() - (this->w + 1))) { end_window++; }

                // If the window is not empty
                if ((start_window < end_window - 1) and (tsp[end_window] > pos_on_reference)) {
                    spdlog::debug("------------------------------------------------------------");
                    spdlog::debug("from {}", sample.get_reference().substr(tsp[start_window - 1], this->w));
                    spdlog::debug("copied from {} to {}", tsp[start_window], tsp[end_window] + this->w);
                    spdlog::debug("next variation: {}", sample_iterator.next_variation());
                    spdlog::debug("skipped phrases: {}", end_window - start_window);

                    // copy from parse[start_window : end_window]
                    out_file.write((char*) &(this->reference_parse->parse[start_window]), sizeof(int) * (end_window - start_window + 1));
                    this->parse_size += end_window - start_window + 1;

                    // move iterators and re-initialize phrase
                    sample_iterator.go_to(tsp[end_window]);
                    phrase.clear();
                    for (std::size_t i = 0; i < this->w; i++) {
                        ++sample_iterator;
                        char next_char  = (params.acgt_only) ? acgt_only_table[*sample_iterator] : *sample_iterator;
                        phrase.push_back(next_char);
                    }

                    kr_hash.reset(); kr_hash.initialize(phrase.data(), params.w);

                    ++sample_iterator;
                    spdlog::debug("New phrase [{}]: {}", phrase.size(), std::string((char*) phrase.data(), phrase.size()));
                    spdlog::debug("------------------------------------------------------------");
                }
            }
        }

        # Next phrase should contain a variation so parse as normal, also if we don't
        # want to use the acceleration we should always end up here
        char next_char  = (params.acgt_only) ? acgt_only_table[*sample_iterator] : *sample_iterator;
        phrase.push_back(next_char);
        kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]);
        ++sample_iterator;

        if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0)) {
            int hash = this->dictionary->check_and_add(phrase);

            out_file.write((char*) (&hash), sizeof(int)); this->parse_size += 1;

            if (phrase[0] != DOLLAR_PRIME) {
                spdlog::debug("------------------------------------------------------------");
                spdlog::debug("Parsed phrase [{}] {}", phrase.size(), std::string((char*) phrase.data(), phrase.size()));
                spdlog::debug("------------------------------------------------------------");
            }

            phrase.erase(phrase.begin(), phrase.end() - this->w); // Keep the last w chars

            kr_hash.reset(); kr_hash.initialize(phrase.data(), params.w);
        }
    }

    # Last phrase
    if (phrase.size() >= this->params.w) {
        // append w-1 dollar prime and a dollar sequence at the end of each sample
        phrase.insert(phrase.end(), params.w - 1, DOLLAR_PRIME);
        phrase.emplace_back(DOLLAR_SEQUENCE);

        // write down last phrase of last sequence
        int hash = this->dictionary->check_and_add(phrase);
        out_file.write((char*) (&hash), sizeof(int));   this->parse_size += 1;

        // if this is the last sample, add w dollars at the end
        if (sample.last(this->working_genotype)) {
            phrase.erase(phrase.begin(), phrase.end() - this->w); // keep the last w chars
            phrase.insert(phrase.end(), params.w, DOLLAR);

            // write down last phrase
            int hash_l = this->dictionary->check_and_add(phrase);
            out_file.write((char*) (&hash_l), sizeof(int));   this->parse_size += 1;
        }
    }
    else { spdlog::error("A sample doesn't have w dollar prime at the end!"); std::exit(EXIT_FAILURE); }
}

void vcfbwt::pfp::ParserVCF::close() {
    if (closed) return; closed = true;

    if ((tags & MAIN) or (tags & WORKER)) {
        vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
        this->out_file.close();
    }

    // Output parse, substitute hash with rank
    if (tags & MAIN) {
        // close all the registered workers and merge their dictionaries
        spdlog::info("Main parser: closing all registered workers");
        for (auto worker : registered_workers) { worker.get().close(); }

        spdlog::info("Main parser: Replacing hash values with ranks in MAIN, WORKERS and reference");

        // repeat for reference
        if (not this->reference_parse->parse.empty()) {
            for (std::size_t i = 0; i < this->reference_parse->parse.size(); i++) {
                int rank = this->dictionary->hash_to_rank(this->reference_parse->parse[i]);
                this->reference_parse->parse[i] = rank;
            }
        }

        // MAIN read hash values and substitute with ranks
        if (this->parse_size != 0) {
            std::ofstream out_ranks(tmp_out_file_name + ".ranks");
            if (not out_ranks.is_open()) { spdlog::error("Can't open {}", tmp_out_file_name + ".ranks"); std::exit(EXIT_FAILURE); }

            std::ifstream in_hash(tmp_out_file_name);
            if (not in_hash.is_open()) { spdlog::error("Can't open {}", tmp_out_file_name); std::exit(EXIT_FAILURE); }

            for (std::size_t i = 0; i < this->parse_size; i++) {
                int hash;
                in_hash.read((char*) &hash, sizeof(int));
                int rank = this->dictionary->hash_to_rank(hash);
                out_ranks.write((char*) &rank, sizeof(int));
            }
            in_hash.close();
            vcfbwt::DiskWrites::update(out_ranks.tellp());
            out_ranks.close();
        }

        // repeat for every worker
        for (auto worker : registered_workers) {
            if (worker.get().parse_size != 0) {
                std::ofstream out_ranks(worker.get().tmp_out_file_name + ".ranks");
                if (not out_ranks.is_open()) { spdlog::error("Can't open {}", worker.get().tmp_out_file_name + ".ranks"); std::exit(EXIT_FAILURE); }

                std::ifstream in_hash(worker.get().tmp_out_file_name);
                if (not in_hash.is_open()) { spdlog::error("Can't open {}", worker.get().tmp_out_file_name); std::exit(EXIT_FAILURE); }

                for (std::size_t i = 0; i < worker.get().parse_size; i++) {
                    int hash;
                    in_hash.read((char*) &hash, sizeof(int));
                    int rank = this->dictionary->hash_to_rank(hash);
                    out_ranks.write((char*) &rank, sizeof(int));
                }
                in_hash.close();
                vcfbwt::DiskWrites::update(out_ranks.tellp());
                out_ranks.close();
            }
        }

        // Merging files together
        spdlog::info("Main parser: concatenating parsings from workers and reference, reference as first");
        std::ofstream merged(out_file_name, std::ios_base::binary);

        // Reference
        std::size_t out_parse_size = 0;
        out_parse_size += this->reference_parse->parse.size();
        for (auto& e : this->reference_parse->parse) {
            int out_e = e;
            merged.write((char*) &out_e, sizeof(int));
        }

        // Main
        if ((this->parse_size) != 0) {
            std::ifstream main_parse(this->tmp_out_file_name + ".ranks");
            merged << main_parse.rdbuf();
            out_parse_size += this->parse_size;
        }

        // Workers
        for (auto worker : registered_workers) {
            if (worker.get().parse_size != 0) {
                out_parse_size += worker.get().parse_size;
                std::ifstream worker_parse(worker.get().tmp_out_file_name + ".ranks");
                merged << worker_parse.rdbuf();
            }
        }
        vcfbwt::DiskWrites::update(merged.tellp()); // Disk Stats
        merged.close();

        this->parse_size = out_parse_size;

        // Print dicitionary on disk
        spdlog::info("Main parser: writing dictionary to disk NOT COMPRESSED");
        std::string dict_file_name = out_file_prefix + EXT::DICT;
        std::ofstream dict(dict_file_name);

        for (std::size_t i = 0; i < this->dictionary->size(); i++) {
            dict.write((char*) this->dictionary->sorted_entry_at(i).data(), this->dictionary->sorted_entry_at(i).size());
            dict.put(ENDOFWORD);
        }
        dict.put(ENDOFDICT);

        vcfbwt::DiskWrites::update(dict.tellp()); // Disk Stats
        dict.close();

        vcfbwt::pfp::PropertiesWriter<char> properties_out(this->out_file_prefix, this->params);
        properties_out.write();

        spdlog::info("Main parser: closed");
    }
}

void vcfbwt::pfp::ParserFasta::init(const ReferenceParse::Params& params, const std::string& prefix) {
    this->params = params;
    this->w = this->params.w;
    this->p = this->params.p;
    this->parse_size = 0;
    this->out_file_prefix = prefix;
    this->out_file_name = prefix + EXT::PARSE;
    this->tmp_out_file_name = TempFile::getName("parse");
    this->out_file.open(tmp_out_file_name, std::ios::binary);
}

void vcfbwt::pfp::ParserFasta::operator()() {
    # Open input file with kseq
    gzFile fp; kseq_t *record;
    fp = gzopen(this->in_file_path.c_str(), "r");
    if (fp == 0) {
        spdlog::error("Failed to open input file {}", in_file_path);
        exit(EXIT_FAILURE);
    }

    std::vector<char> phrase;
    spdlog::info("Parsing sequence");

    # Karp Robin Hash Function for sliding window
    Mersenne_KarpRabinHash kr_hash(this->params.w);

    # First sequence start with one dollar
    phrase.emplace_back(DOLLAR);


    #write function for fasta file to return name, comment and sequence
    '''
    from Bio import SeqIO

# Function to read a FASTA file and return a list of sequences
def read_fasta(file_path):
    sequences = []
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record)
    return sequences

# Example usage
file_path = "your_fasta_file.fasta"
sequences = read_fasta(file_path)

# Print sequences
for seq_record in sequences:
    print(f"ID: {seq_record.id}")
    print(f"Sequence: {seq_record.seq}\n")

def parse_fasta(file_path):
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

# Example usage
file_path = 'path/to/your/fasta_file.fasta'
sequences = parse_fasta(file_path)

for ref_id, sequence in sequences.items():
    print(f"Reference ID: {ref_id}")
    print(f"Sequence: {sequence}")
    '''



    record = kseq_init(fp);
    while(kseq_read(record) >= 0) {
        std::string sequence_name("<error reading sequence name>"), sequence_comment("<error reading sequence comment>");
        if (record->name.s != NULL) { sequence_name = record->name.s; }
        if (record->comment.s != NULL) { sequence_comment = record->comment.s; }
        this->sequences_processed.push_back(sequence_name + " " + sequence_comment);
        spdlog::info("Parsed:\t{}", sequence_name + " " + sequence_comment);

        # Previous last phrase
        if (phrase[0] != DOLLAR and phrase.size() >= this->params.w) {
            // Append w-1 dollar prime, and one dollar seq at the end of each sequence
            for (std::size_t j = 0; j < this->params.w - 1; j++) { phrase.emplace_back(DOLLAR_PRIME); }
            phrase.emplace_back(DOLLAR_SEQUENCE);

            int hash = this->dictionary.check_and_add(phrase);

            this->parse.push_back(hash);
            this->trigger_strings_position.push_back(record->seq.l - 1);

            phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars

            kr_hash.reset(); kr_hash.initialize(phrase.data(), params.w);
        }

        for (std::size_t seq_it = 0; seq_it < record->seq.l; seq_it++) {
            char c = record->seq.s[seq_it];

            if (c <= DOLLAR_PRIME) {
                spdlog::error("Input may not contain bytes with integer value less than or equal to 5!");
                std::exit(EXIT_FAILURE);
            }
            if (params.acgt_only) { c = acgt_only_table[c]; }

            phrase.push_back(c);
            if (phrase.size() == params.w) { kr_hash.initialize(phrase.data(), params.w); }
            else if (phrase.size() > params.w) { kr_hash.update(phrase[phrase.size() - params.w - 1], phrase[phrase.size() - 1]); }

            if ((phrase.size() > this->params.w) and ((kr_hash.get_hash() % this->params.p) == 0)) {
                int hash = this->dictionary.check_and_add(phrase);

                out_file.write((char*) (&hash), sizeof(int)); this->parse_size += 1;

                phrase.erase(phrase.begin(), phrase.end() - this->params.w); // Keep the last w chars

                kr_hash.reset(); kr_hash.initialize(phrase.data(), params.w);
            }
        }
    }

    # Last phrase
    if (phrase.size() >= this->params.w) {
        // append w-1 dollar prime and a dollar sequence at the end of each sample
        phrase.insert(phrase.end(), params.w - 1, DOLLAR_PRIME);
        phrase.emplace_back(DOLLAR_SEQUENCE);

        // write down last phrase of last sequence
        int hash = this->dictionary.check_and_add(phrase);
        out_file.write((char*) (&hash), sizeof(int));   this->parse_size += 1;

        // last sequence, add w dollars at the end
        phrase.erase(phrase.begin(), phrase.end() - this->params.w); // keep the last w chars
        phrase.insert(phrase.end(), this->params.w, DOLLAR);

        // write down last phrase
        int hash_l = this->dictionary.check_and_add(phrase);
        out_file.write((char*) (&hash_l), sizeof(int));   this->parse_size += 1;
    }
    else { spdlog::error("Missing w DOLLAR at the end!"); std::exit(EXIT_FAILURE); }

    kseq_destroy(record);
    gzclose(fp);
}

void vcfbwt::pfp::ParserFasta::close() {
    if (closed) return; closed = true;

    vcfbwt::DiskWrites::update(out_file.tellp()); // Disk Stats
    this->out_file.close();

    spdlog::info("Main parser: Replacing hash values with ranks.");

    // mmap file and substitute
    if (this->parse_size != 0) {
        std::ofstream out_ranks(this->out_file_name);
        if (not out_ranks.is_open()) { spdlog::error("Can't open {}", this->out_file_name); std::exit(EXIT_FAILURE); }

        std::ifstream in_hash(tmp_out_file_name);
        if (not in_hash.is_open()) { spdlog::error("Can't open {}", tmp_out_file_name); std::exit(EXIT_FAILURE); }

for i in range(0, self.parse_size):
    hash = hash_type()
    in_hash.readinto(hash)
    rank = self.dictionary.hash_to_rank(hash)
    out_ranks.write(rank)

in_hash.close()
vcfbwt.DiskWrites.update(out_ranks.tell())
out_ranks.close()


# Print dictionary on disk
spdlog.info("Main parser: writing dictionary on disk NOT COMPRESSED")
dict_file_name = out_file_prefix + EXT.DICT
with open(dict_file_name, 'wb') as dict_file:
    for i in range(len(self.dictionary)):
        dict_file.write(self.dictionary.sorted_entry_at(i))
        dict_file.write(ENDOFWORD.to_bytes(1, byteorder='big'))
    dict_file.write(ENDOFDICT.to_bytes(1, byteorder='big'))

vcfbwt.DiskWrites.update(dict_file.tell()) # Disk Stats

vcfbwt.pfp.PropertiesWriter(self.out_file_prefix, self.params).write()

spdlog.info("Main parser: closed")python
def init(self, params, prefix):
    self.params = params
    self.w = self.params.w
    self.p = self.params.p
    self.parse_size = 0
    self.out_file_prefix = prefix
    self.out_file_name = prefix + EXT.PARSE
    self.tmp_out_file_name = TempFile.getName("parse")
    self.out_file = open(self.tmp_out_file_name, 'wb')

def __call__(self):
    # Open input file with gzip
    with gzip.open(self.in_file_path, 'rb') as fp:
        phrase = [DOLLAR]
        spdlog.info(f"Parsing {self.in_file_path}")
        
        # Karp Robin Hash Function for sliding window
        kr_hash = Mersenne_KarpRabinHash(self.params.w)
        
        while True:
            c = fp.read(1)
            if not c:
                break
            c = ord(c)
            if c <= DOLLAR_PRIME:
                spdlog.error("Input may not contain bytes with integer value less than or equal to 5!")
                exit(EXIT_FAILURE)
            
            phrase.append(c)
            if len(phrase) == self.params.w:
                kr_hash.initialize(phrase, self.params.w)
            elif len(phrase) > self.params.w:
                kr_hash.update(phrase[len(phrase) - self.params.w - 1], phrase[len(phrase) - 1])
            
            if len(phrase) > self.params.w and (kr_hash.get_hash() % self.params.p) == 0:
                hash_value = self.dictionary.check_and_add(phrase)
                self.out_file.write(hash_value.to_bytes(sizeof(hash_type), byteorder='big'))
                self.parse_size += 1
                
                phrase = phrase[len(phrase) - self.params.w:]
                
                kr_hash.reset()
                kr_hash.initialize(phrase, self.params.w)
        
        # Last phrase
        if len(phrase) >= self.params.w:
            # Append w dollar at the end
            phrase.extend([DOLLAR] * self.params.w)
            
            hash_value = self.dictionary.check_and_add(phrase)
            
            self.out_file.write(hash_value.to_bytes(sizeof(hash_type), byteorder='big'))
            self.parse_size += 1
        else:
            spdlog.error("A sequence doesn't have w DOLLAR at the end!")
            exit(EXIT_FAILURE)

def close(self):
    if self.closed:
        return
    self.closed = True
    
    vcfbwt.DiskWrites.update(self.out_file.tell()) # Disk Stats
    self.out_file.close()
    
    # Occurrences
    spdlog.info("Main parser: Replacing hash values with ranks.")
    
    # mmap file and substitute
    if self.parse_size != 0:
        with open(self.out_file_name, 'wb') as out_ranks, open(self.tmp_out_file_name, 'rb') as in_hash:
            for i in range(self.parse_size):
                hash_value = int.from_bytes(in_hash.read(sizeof(hash_type)), byteorder='big')
                rank = self.dictionary.hash_to_rank(hash_value)
                out_ranks.write(rank.to_bytes(sizeof(size_type), byteorder='big'))
        vcfbwt.DiskWrites.update(out_ranks.tell())
    
    # Print dictionary on disk
    spdlog.info("Main parser: writing dictionary on disk NOT COMPRESSED")
    dict_file_name = self.out_file_prefix + EXT.DICT
    with open(dict_file_name, 'wb') as dict_file:
        for i in range(len(self.dictionary)):
            dict_file.write(self.dictionary.sorted_entry_at(i))
            dict_file.write(ENDOFWORD.to_bytes(1, byteorder='big'))
        dict_file.write(ENDOFDICT.to_bytes(1, byteorder='big'))
    
    vcfbwt.DiskWrites.update(dict_file.tell()) # Disk Stats
    
    vcfbwt.pfp.PropertiesWriter(self.out_file_prefix, self.params).write()

    spdlog.info("Main parser: closed")python
def init(self, params, prefix):
    self.params = params
    self.w = self.params.w
    self.p = self.params.p
    self.parse_size = 0
    self.out_file_prefix = prefix
    self.out_file_name = prefix + EXT.PARSE
    self.tmp_out_file_name = TempFile.getName("parse")
    self.out_file = open(self.tmp_out_file_name, 'wb')

def __call__(self):
    # Open input file with gzip
    with gzip.open(self.in_file_path, 'rb') as fp:
        phrase = [DOLLAR]
        spdlog.info(f"Parsing {self.in_file_path}")
        
        # Karp Robin Hash Function for sliding window
        kr_hash = Mersenne_KarpRabinHash4(self.params.w * 4)
        
        while True:
            c = fp.read(4)
            if not c:
                break
            c = int.from_bytes(c, byteorder='big') + self.params.integers_shift
            phrase.append(c)
            if len(phrase) == self.params.w:
                kr_hash.initialize(bytes(phrase), len(phrase) * sizeof(uint32_t))
            elif len(phrase) > self.params.w:
                kr_hash.update(bytes([phrase[len(phrase) - self.params.w - 1]]), bytes([phrase[len(phrase) - 1]]))
            
            if len(phrase) > self.params.w and (kr_hash.get_hash() % self.params.p) == 0:
                hash_value = self.dictionary.check_and_add(phrase)
                self.out_file.write(hash_value.to_bytes(sizeof(hash_type), byteorder='big'))
                self.parse_size += 1
                
                phrase = phrase[len(phrase) - self.params.w:]
                
                kr_hash.reset()
                kr_hash.initialize(bytes(phrase), len(phrase) * sizeof(uint32_t))
        
        # Last phrase
        if len(phrase) >= self.params.w:
            # Append w dollar at the end
            phrase.extend([DOLLAR] * self.params.w)
            
            hash_value = self.dictionary.check_and_add(phrase)
            
            self.out_file.write(hash_value.to_bytes(sizeof(hash_type), byteorder='big'))
            self.parse_size += 1
        else:
            spdlog.error("A sequence doesn't have w DOLLAR at the end!")
            exit(EXIT_FAILURE)

def close(self):
    if self.closed:
        return
    self.closed = True
    
    vcfbwt.DiskWrites.update(self.out_file.tell()) # Disk Stats
    self.out_file.close()
    
    spdlog.info("Main parser: Replacing hash values with ranks.")
    
    # mmap file and substitute
    if self.parse_size != 0:
        with open(self.out_file_name, 'wb') as out_ranks, open(self.tmp_out_file_name, 'rb') as in_hash:
            for i in range(self.parse_size):
                hash_value = int.from_bytes(in_hash.read(sizeof(hash_type)), byteorder='big')
                rank = self.dictionary.hash_to_rank(hash_value)
                out_ranks.write(rank.to_bytes(sizeof(size_type), byteorder='big'))
        vcfbwt.DiskWrites.update(out_ranks.tell())
    
    # Print dictionary on disk
    spdlog.info("Main parser: writing dictionary on disk NOT COMPRESSED")
    dict_file_name = self.out_file_prefix + EXT.DICT
    with open(dict_file_name, 'wb') as dict_file:
        for i in range(len(self.dictionary)):
            dict_file.write(bytes(self.dictionary.sorted_entry_at(i)))
            dict_file.write(ENDOFWORD.to_bytes(sizeof(uint32_t), byteorder='big'))
        dict_file.write(ENDOFDICT.to_bytes(sizeof(uint32_t), byteorder='big'))
    
    vcfbwt.DiskWrites.update(dict_file.tell()) # Disk Stats
    
    vcfbwt.pfp.PropertiesWriter(self.out_file_prefix, self.params).write()

    spdlog.info("Main parser: closed")

import gzip
import logging
import struct
import os

class ParserIntegers:
    def __init__(self, in_file_path, out_file_name, tmp_out_file_name, out_file_prefix, params, dictionary):
        self.in_file_path = in_file_path
        self.out_file_name = out_file_name
        self.tmp_out_file_name = tmp_out_file_name
        self.out_file_prefix = out_file_prefix
        self.params = params
        self.dictionary = dictionary
        self.parse_size = 0
        self.closed = False
        self.out_file = open(tmp_out_file_name, 'wb')

    def __call__(self):
        # Open input file with gzip
        try:
            fp = gzip.open(self.in_file_path, 'rb')
        except OSError:
            logging.error(f"Failed to open input file {self.in_file_path}")
            exit(1)

        phrase = []
        logging.info(f"Parsing {self.in_file_path}")

        # Karp Robin Hash Function for sliding window
        kr_hash = Mersenne_KarpRabinHash4(self.params.w * 4)

        # First sequence start with one dollar
        phrase.append(DOLLAR)

        while True:
            c = fp.read(4)
            if not c:
                break
            c = struct.unpack('I', c)[0]
            phrase.append(c + self.params.integers_shift)
            if len(phrase) == self.params.w:
                kr_hash.initialize(bytearray(phrase))
            elif len(phrase) > self.params.w:
                kr_hash.update(bytearray(phrase[-self.params.w-1:]), bytearray(phrase[-1:]))

            if len(phrase) > self.params.w and (kr_hash.get_hash() % self.params.p) == 0:
                hash_value = self.dictionary.check_and_add(phrase)
                self.out_file.write(struct.pack('I', hash_value))
                self.parse_size += 1

                phrase = phrase[-self.params.w:]  # Keep the last w chars

                kr_hash.reset()
                kr_hash.initialize(bytearray(phrase))

        # Last phrase
        if len(phrase) >= self.params.w:
            # Append w dollar at the end
            phrase.extend([DOLLAR] * self.params.w)

            hash_value = self.dictionary.check_and_add(phrase)
            self.out_file.write(struct.pack('I', hash_value))
            self.parse_size += 1
        else:
            logging.error("A sequence doesn't have w DOLLAR at the end!")
            exit(1)

        fp.close()

    def close(self):
        if self.closed:
            return
        self.closed = True

        vcfbwt.DiskWrites.update(self.out_file.tell())
        self.out_file.close()

        logging.info("Main parser: Replacing hash values with ranks.")

        # mmap file and substitute
        if self.parse_size != 0:
            with open(self.out_file_name, 'wb') as out_ranks:
                with open(self.tmp_out_file_name, 'rb') as in_hash:
                    for _ in range(self.parse_size):
                        hash_value = struct.unpack('I', in_hash.read(4))[0]
                        rank = self.dictionary.hash_to_rank(hash_value)
                        out_ranks.write(struct.pack('I', rank))

                vcfbwt.DiskWrites.update(out_ranks.tell())

        # Print dictionary on disk
        logging.info("Main parser: writing dictionary on disk NOT COMPRESSED")
        dict_file_name = self.out_file_prefix + EXT.DICT
        with open(dict_file_name, 'wb') as dict_file:
            for i in range(self.dictionary.size()):
                entry = self.dictionary.sorted_entry_at(i)
                dict_file.write(bytearray(entry))
                dict_file.write(struct.pack('I', ENDOFWORD))
            dict_file.write(struct.pack('I', ENDOFDICT))

        vcfbwt.DiskWrites.update(dict_file.tell())

        properties_out = vcfbwt.pfp.PropertiesWriter(self.out_file_prefix, self.params)
        properties_out.write()

        logging.info("Main parser: closed")



#The Python code translates the C++ code for writing the dictionary to disk and parsing text/integer files. It uses the `open` function to create file objects, and the `write` method to write data to files. The `bytes` and `int.from_bytes` functions are used to convert between byte strings and integers. The `gzip` module is used to read compressed input files. The code also imports any necessary libraries like `gzip`.

'''
  SUMMARY:
  This script is part of a bioinformatics or computational biology software toolkit,
  designed to manage and parse genomic data efficiently using custom hashing algorithms.

  NAMESPACE AND CLASSES:
  - Namespace `vcfbwt::pfp` includes several classes for handling genomic data parsing tasks.
  - `Mersenne_KarpRabinHash` class implements a hashing mechanism for sliding window analysis in genomic sequences.
  - `ReferenceParse` and `ParserVCF` classes manage the parsing of reference genomes or variant call data.

  FUNCTIONALITY:
  - `init()`: Prepares the environment for parsing operations by setting parameters and initializing output files.
  - `operator()()`: Processes parsing operations, applying hash-based techniques to manage and detect sequences.
  - `close()`: Finalizes the operations, ensuring data integrity and performing clean-up tasks.

  HASHING AND DATA MANAGEMENT:
  - Utilizes a Karp-Rabin hashing algorithm to efficiently identify and manage repeated sequences or patterns.
  - Employs dynamic data structures for storing intermediate and final parsed data, enhancing memory management.

  LOGGING AND ERROR HANDLING:
  - Integrates `spdlog` for detailed logging, aiding in troubleshooting and operational transparency.
  - Implements robust error checking to handle potential data format and operational issues.

  FILE OPERATIONS:
  - Manages both input and output file operations, including reading genomic data and writing outputs to disk.
  - Uses temporary files for handling intermediate data, ensuring data consistency and efficient disk usage.

  PURPOSE:
  - Aimed at accelerating and optimizing the parsing of genomic data, particularly for applications involving large datasets.
  - Helps in identifying genomic variants, sequence repeats, and other significant features with high computational efficiency.

  NOTE:
  - The provided code is conceptual and requires specific adaptations to be directly executable in a real-world application.
  - Combines C++-style structure with Python-like annotations, reflecting an instructional or prototype development approach.
'''


#ALL NEEDS TO BE REVIEWED

'''
Sources

https://github.com/marco-oliva/pfp/blob/master/pfp%2B%2B.cpp
https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm
chat.openai.com for some C++ conversions
https://brilliant.org/wiki/rabin-karp-algorithm/#:~:text=The%20Rabin%2DKarp%20algorithm%20is,important%20application%20of%20computer%20science.
https://www.geeksforgeeks.org/rabin-karp-algorithm-for-pattern-searching/

'''