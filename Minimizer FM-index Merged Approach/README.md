# NanoMapper
A Mapping tool using gapped minimizer and BWT FM-index

## Header Files:

 - `kseq.h`: This is a library used for faster input output from fasta file.
 - `predefines.h`: All necessary header files and global variables declared through it.
 - `biofunc.h`: General purpose sequence related functions are kept here.
 - `minimizers.h`: Minimizer related functions are written in this file.
 - `ModifiedFM.h`: Modifications in FM index API functions are done here. Most of them are copied from the API source code and then did necessary modification with edited function names.

## Functions:

### `bio_func.h`:

 - `long long pat2num(string s)`: Gets a DNA sequence (with or without gap)and returns corresponding 64-bit signed integer.

 - `string num2pat(long long num, int k)`: Gets a 64-bit signed integer and the length of k-mer and converts it to corresponding DNA sequence without inserting gap.

 - `string num2pat2(long long num, int k)`: Gets a 64-bit signed integer and the length of k-mer and converts it to corresponding DNA sequence with inserting gap.

 - `int nt2num(char &x)`: Gets the address of a base and converts the base to a binary number ranging 0-3.

 - `char num2nt(int &x)`: Gets the address of a base in binary and converts the base to a character base A,T,C or G.

 - `long long rev_comp(long long x, int k)`: Gets a 64-bit signed integer with its length k. The integer represents a k-mer. The function returns a 64-bit signed integer which represents the reverse complement of the given k-mer(integer).

 - `void insert_gap(string &s, const int init = 0)`: It takes a DNA sequence and inserts gaps in every third base. An optional parameter could be supplied for offset adjusting (offset would be adjusted having modulo 3).

 - `void reverse_comp_of_ref()`: Initializes the REV_REFF, the reverse complement of the forward sequence REFF.

 - `string insert_gap_read(string str)`[Unused]: The function was for any previous version. There is no use of this function here.

 ### `minimizers.h`:

 - `void find_minimizers(string s, unordered_map<long long, pair<size_t, size_t> > &ret, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fm_ind)`:

 	- `s`: Reference DNA sequence which should be indexed.
 	- `ret`: `unordered_map` of minimizers which would be produced by this function as a result.
 	- `fm_ind`: FM-index of the reference DNA sequence `s` [first parameter]

 - `void find_kmers_of_read_in_refer(int flag, int offset, int start_in_ref, string read_name, string red_seq)`[Unused]: Kept for safety check. No use of this function in this version.

 ###`ModifiedFM.h`:
 Parameters of the functions under this header file are so long. So, only the name and indication for unique identification is used.

 - `b_search(7 parameters)`: No change, same as *SDSL-lite*'s `backward_search()` function. Kept here for better understanding. It has an overloaded version.
 
 - `b_search(11 parameters)`: Modified in a way so that it could report last non-zero k-mer with it's occurrences.

 - `b_search_lim()`: A clone of the 11 parameterized `b_search()` function. But it iterates `K` length for keeping hash in `unordered_map`.

 - `pair<size_t, size_t> get_hashes(string pat, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fm_ind)`: Takes a k-mer and indexed FM-index. It iterates the FM-index and returns the hash-values for `unorderd_map`.

 - `my_locate_hash()`: A clone of `locate()` function of *SDSL-lite*. But it is modified in a way so that it could iterate taking the hashes from unordered map and returns non-empty set of locations for the largest possible k-mer.

 ###`predefines.h`:

 - `init()`: Initializes the largest possible mask to filter any larger mask.

 ### `Main_minimizer.cpp`:

 - `void takeReference(string &ref_filepath)`: Takes the file path of the reference fasta and stores the reference in REFF after retrieving. It also generates the reverse complement of the reference and stores in REV_REFF. It inserts gaps too.

 - `void takeReads(string &read_filepath, vector<single_read>&reads)`: Takes the file path of the read sequences and an address of a vector. It retrieves the reads and stores them in the vector.

 - `void indexReferenceForEnhancedFM(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_f, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_r)`: It takes two index objects of FM-index. One for indexing the original sequence and other is for reverse complement. It index REFF and REV_REFF and reports the summary.

 - `void processEnhancedFMAprroach(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>&  fm_f, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_r, string location_out, vector<single_read>&reads)`: It takes two FM-index objects for searching in them, some reads in a vector to process and an output file path to store the mapping. It does the mapping process and stores in the file provided.

 - `void buildMinimizerIndex(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fm_ind_fwd, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fm_ind_rev)`: It takes to FM-indexed objects for searching and getting hash values. It goes through the reference and reverse complement of the reference. Stores the hashes of the minimizers in `unordered_map`.

 - `int main(int argc, const char **argv)`: Most ancient and prestigious function which leads all of the process.

## Developer
Enamul Hassan and Md. Khairullah Gaurab under the direction of Ruhul Amin Shajib Sir

## Contact Information
Enamul Hassan : enamsustcse@gmail.com

Md. Khairullah Gaurab : mkgaurabsarkar@gmail.com

