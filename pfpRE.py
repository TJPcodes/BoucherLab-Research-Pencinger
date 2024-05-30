import argparse
import os
import spdlog
from vcfbwt import Version, TempFile, VCF, pfp

def main():
    parser = argparse.ArgumentParser(description="PFP++")
    
    parser.add_argument("-v", "--vcf", type=str, nargs='+', help="List of comma ',' separated vcf files. Assuming in genome order!")
    parser.add_argument("-r", "--ref", type=str, nargs='+', help="List of comma ',' separated reference files. Assuming in genome order!")
    parser.add_argument("-f", "--fasta", type=str, help="Fasta file to parse.")
    parser.add_argument("-i", "--int32t", type=str, help="Integers file to parse.")
    parser.add_argument("--int-shift", type=int, default=0, help="Each integer i in int32t input is interpreted as (i + int-shift).", choices=range(0, 201))
    parser.add_argument("-H", "--haplotype", type=str, default="1", help="Haplotype: [1,2,12].")
    parser.add_argument("-t", "--text", type=str, help="Text file to parse.")
    parser.add_argument("-o", "--out-prefix", type=str, help="Output prefix.")
    parser.add_argument("-m", "--max", type=int, default=0, help="Max number of samples to analyze.")
    parser.add_argument("-S", "--samples", type=str, help="File containing the list of samples to parse.")
    parser.add_argument("-w", "--window-size", type=int, default=0, help="Sliding window size.", choices=range(3, 201))
    parser.add_argument("-p", "--modulo", type=int, default=0, help="Modulo used during parsing.", choices=range(5, 20001))
    parser.add_argument("-j", "--threads", type=int, default=1, help="Number of threads.")
    parser.add_argument("--tmp-dir", type=str, help="Temporary files directory.")
    parser.add_argument("-c", "--compress-dictionary", action='store_true', help="Also output compressed the dictionary.")
    parser.add_argument("--use-vcf-acceleration", action='store_true', help="Use reference parse to avoid re-parsing.")
    parser.add_argument("--print-statistics", action='store_true', help="Print out csv containing stats.")
    parser.add_argument("--output-occurrences", action='store_true', help="Output count for each dictionary phrase.")
    parser.add_argument("--output-sai", action='store_true', help="Output sai array.")
    parser.add_argument("--output-last", action='store_true', help="Output last array.")
    parser.add_argument("--acgt-only", action='store_true', help="Convert all non ACGT characters from a VCF or FASTA file to N.")
    parser.add_argument("--verbose", action='store_true', help="Verbose output.")
    parser.add_argument("--version", action='version', version=Version.print())
    
    args = parser.parse_args()
    
    if args.verbose:
        spdlog.set_level(spdlog.level.debug)
    
    vcfs_file_names = [f for f in args.vcf if f]
    refs_file_names = [f for f in args.ref if f]
    
    spdlog.info("Current Configuration:\n{}", vars(args))
    
    if args.tmp_dir:
        TempFile.setDirectory(args.tmp_dir)
    
    params = pfp.Params(
        integers_shift=args.int_shift,
        w=args.window_size,
        p=args.modulo,
        compress_dictionary=args.compress_dictionary,
        use_acceleration=args.use_vcf_acceleration,
        print_out_statistics_csv=args.print_statistics,
        output_occurrences=args.output_occurrences,
        output_sai=args.output_sai,
        output_last=args.output_last,
        acgt_only=args.acgt_only
    )
    
    if args.fasta:
        out_prefix = args.out_prefix or args.fasta
        main_parser = pfp.ParserFasta(params, args.fasta, out_prefix)
        main_parser()
        main_parser.close()
    elif args.text:
        out_prefix = args.out_prefix or args.text
        main_parser = pfp.ParserText(params, args.text, out_prefix)
        main_parser()
        main_parser.close()
    elif args.int32t:
        out_prefix = args.out_prefix or args.int32t
        main_parser = pfp.ParserIntegers(params, args.int32t, out_prefix)
        main_parser()
        main_parser.close()
    else:
        last_genotype = 1 if args.haplotype in ["2", "12"] else 0
        vcf = VCF(refs_file_names, vcfs_file_names, args.samples, args.max, last_genotype)
        reference_parse = pfp.ReferenceParse(vcf.get_reference(), params)
        main_parser = pfp.ParserVCF(params, args.out_prefix, reference_parse)
        
        # Will need to be edited for GPU (line 80 - 101)

        workers = [pfp.ParserVCF(params, "", reference_parse, pfp.ParserVCF.WORKER | pfp.ParserVCF.UNCOMPRESSED) for _ in range(args.threads)]
        for worker in workers:
            main_parser.register_worker(worker)
        
        
        if args.haplotype in ["1", "2"]:
            genotype = 0 if args.haplotype == "1" else 1
            for worker in workers:
                worker.set_working_genotype(genotype)
            for i in range(len(vcf)):
                spdlog.info("Processing sample [{}/{} H{}]: {}", i, len(vcf), args.haplotype, vcf[i].id())
                workers[i % args.threads](vcf[i])
        elif args.haplotype == "12":
            for i in range(len(vcf)):
                spdlog.info("Processing sample [{}/{} H1]: {}", i, len(vcf), vcf[i].id())
                workers[i % args.threads].set_working_genotype(0)
                workers[i % args.threads](vcf[i])
                spdlog.info("Processing sample [{}/{} H2]: {}", i, len(vcf), vcf[i].id())
                workers[i % args.threads].set_working_genotype(1)
                workers[i % args.threads](vcf[i])
        
        main_parser.close()

if __name__ == "__main__":
    main()


"""
SUMMARY:
This Python script provides a command-line interface (CLI) for parsing genomic data from various input formats such as VCF, 
FASTA, plain text, and integer arrays. It utilizes advanced parsing techniques to efficiently process and manage genomic data.

FEATURES:
- Supports multiple input types including VCF, FASTA, plain text, and integer array files.
- Allows specification of haplotypes and various genomic parsing options through command-line arguments.
- Includes advanced functionalities like VCF acceleration, output compression, and detailed statistical outputs.
- Offers extensive configuration options for window size, hash modulus, and temporary directory management.
- Utilizes multi-threading to enhance processing speed and efficiency.

FUNCTIONALITY:
- Parses command-line arguments to configure the parsing operation.
- Depending on the input type specified, it initializes the appropriate parser and processes the input file.
- Supports conditional execution paths based on the input type and provided arguments.
- Outputs parsed data, statistical information, and potentially the suffix array index (SAI) and last array for genomic positions.

IMPLEMENTATION DETAILS:
- Uses argparse for command-line argument parsing.
- Employs spdlog for logging and debug information.
- Handles temporary file management through a custom TempFile class.
- Utilizes custom classes and methods from the 'vcfbwt' module to perform the core parsing and data handling tasks.

USAGE:
- The script is intended to be used as a command-line tool where users can specify the type of genomic data to parse, along with various 
processing parameters.
- It provides detailed logging and error handling to ensure robust operation.

EXAMPLE:
To run the script, use a command in the following format:
python script_name.py --vcf file1.vcf,file2.vcf --ref file1.ref,file2.ref --haplotype 1 --window-size 100 --modulo 101 --out-prefix output_prefix

NOTE:
- This script is part of a larger bioinformatics toolkit and is designed for high-throughput genomic data processing.
- The actual parsing logic and data handling are abstracted within the 'vcfbwt' module's classes and functions.
"""

'''
Sources:

https://github.com/marco-oliva/pfp/blob/master/pfp%2B%2B.cpp
https://github.com/marco-oliva/pfp/blob/master/pfp%2B%2B.cpp
https://en.wikipedia.org/wiki/Rabin%E2%80%93Karp_algorithm#:~:text=The%20Rabin%E2%80%93Karp%20algorithm%20proceeds,full%20comparison%20at%20that%20position
https://almob.biomedcentral.com/articles/10.1186/s13015-019-0148-5
"Recursive Prefix-Free Parsing" Christina Boucher Dept of Computer and 
Information Science and Engineering Herbert Wertheim College of Engineering PowerPoint

'''