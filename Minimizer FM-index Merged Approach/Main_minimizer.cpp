/// Necessary header files, functions are included in "predefines.h"
#include "predefines.h"

void takeReference(string &ref_filepath)
{
	gzFile fp;
	kseq_t *seq;
	FILE* ttp = fopen((char*)ref_filepath.c_str(), "r");
	fp = gzdopen(fileno(ttp), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0)
	{
		REFF = (string)seq->seq.s;
		string ref_name = (string)seq->name.s;
		ref_len = seq->seq.l;
	}
	kseq_destroy(seq);
	gzclose(fp);
	fclose(ttp);
	cerr << "Length of Reference: " << ref_len << endl;
	cout << "Reference Length: " << ref_len << endl;
	reverse_comp_of_ref(); // must be called before inserting gap
	insert_gap(REFF);
	//cout << REFF << endl;
	insert_gap(REV_REFF); // must be called after creating reverse complement by calling reverse_com_of_ref()
	// cout << REV_REFF << endl;
}

void takeReads(string &read_filepath, vector<single_read>&reads)
{
	gzFile fp;
	kseq_t *seq;
	FILE* ttp = fopen((char*)read_filepath.c_str(), "r");
	fp = gzdopen(fileno(ttp), "r");
	seq = kseq_init(fp);
	int file_cnter = 0;
	while (kseq_read(seq) >= 0)
	{
		reads.push_back(seq);

		file_cnter++;
	}

	kseq_destroy(seq);
	gzclose(fp);
	fclose(ttp);
}

/**
 * fm_f - To index forward reference
 * fm_r - To index backward reference
 */
void indexReferenceForEnhancedFM(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_f, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_r)
{
	reverse(REFF.begin(), REFF.end());
	reverse(REV_REFF.begin(), REV_REFF.end());
	clock_gettime(CLOCK_MONOTONIC, &start);
	construct_im(fm_f, REFF, 1);
	construct_im(fm_r, REV_REFF, 1);
	clock_gettime(CLOCK_MONOTONIC, &finish);
	reverse(REFF.begin(), REFF.end());
	reverse(REV_REFF.begin(), REV_REFF.end());
	int elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
	cout << " >>> INDEX CREATING <<<\n";
	cout << "Indexing of Reverse Reference (FOR ENHANCED FM):\n";
	cout << "Length of Reference: " << REFF.size() << "\n";
	cout << "Reference Indexing Time: " << (double)elapsed << "\n";
	cout << "Memory Consumption in MB:\n";
	cout << "Forward = " << size_in_mega_bytes(fm_f) << "\n";
	cout << "Reverse = " << size_in_mega_bytes(fm_r) << "\n";
	cout << " >>> END - INDEX CREATING <<<\n";
}

/**
 * fm_f - indexed forward reference
 * fm_r - indexed backward reference
 * location_out - file_path of output
 * reads - vector of reads
 */
void processEnhancedFMAprroach(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>&  fm_f, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_r, string location_out, vector<single_read>&reads)
{
	FILE* loc_out = fopen((char*)location_out.c_str(), "w");
	string kmer_with_gap[3], temp_s;
	double elapsed = 0;
	int gap_counter, L, st;
	long long to_num;
	size_t k_length, hash_start, hash_end;
	int_vector<64u>locations;
	/// Processing time in modified procedure
	for (auto x : reads)
	{
		if (x.alignment)
		{
			L = x.seq.size();
			for (int temp_i = 0; temp_i < 3; temp_i++)
			{
				kmer_with_gap[temp_i] = x.seq;
				insert_gap(kmer_with_gap[temp_i], temp_i);
				reverse(kmer_with_gap[temp_i].begin(), kmer_with_gap[temp_i].end());
			}

			fprintf(loc_out, "> %s\n", x.name.c_str());
			gap_counter = 0;
			for (int ii = 0 ; ii < L; ++ii, gap_counter++)
			{
				gap_counter = (ii % 3); // to maintain the index of kmer_with_gap
				st = ii;
				if (x.forward)
				{
					temp_s = string(kmer_with_gap[gap_counter].end() - ii - min_k, kmer_with_gap[gap_counter].end() - ii);
					reverse(temp_s.begin(), temp_s.end());
					// cerr << temp_s << endl;

					to_num = pat2num(temp_s);
					if (retf.count(to_num) == 0) continue;
					/// ATTENTION:  (1) first element should be omitted,
					///             (2) actual location could be found by subtructing from reference length
					k_length = 0;
					hash_start = retf[to_num].first;
					hash_end = retf[to_num].second;
					clock_gettime(CLOCK_MONOTONIC, &start);
					locations = my_locate_hash(fm_f, kmer_with_gap[gap_counter].begin(), kmer_with_gap[gap_counter].end() - ii,
					                           min_k, max_count, k_length, hash_start, hash_end);
					clock_gettime(CLOCK_MONOTONIC, &finish);
					elapsed += (finish.tv_sec - start.tv_sec);
					elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
					if (k_length >= min_k and SZ(locations) != 0)
					{
						ii += k_length;
						--ii;
					}
					// freq_f[occs]++;
				}
				else
				{
					temp_s = string(kmer_with_gap[gap_counter].end() - ii - min_k, kmer_with_gap[gap_counter].end() - ii);
					reverse(temp_s.begin(), temp_s.end());
					// cerr << temp_s << endl;
					to_num = pat2num(temp_s);
					if (retr.count(to_num) == 0) continue;
					/// ATTENTION: See the note in the code of if section.
					k_length = 0;
					hash_start = retf[to_num].first;
					hash_end = retf[to_num].second;
					clock_gettime(CLOCK_MONOTONIC, &start);
					locations = my_locate_hash(fm_r, kmer_with_gap[gap_counter].begin(), kmer_with_gap[gap_counter].end() - ii,
					                           min_k, max_count, k_length, hash_start, hash_end);
					clock_gettime(CLOCK_MONOTONIC, &finish);
					elapsed += (finish.tv_sec - start.tv_sec);
					elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);

					if (k_length >= min_k and SZ(locations) != 0)
					{
						ii += k_length;
						--ii;
					}
				}
				if (k_length >= min_k and SZ(locations) != 0)
				{
					temp_s = string(kmer_with_gap[gap_counter].end() - ii - 1, kmer_with_gap[gap_counter].end() - st);
					// insert_gap(temp_s);
					reverse(temp_s.begin(), temp_s.end());
					fprintf(loc_out, "%d - %d : %s (%d)\n", st, ii, temp_s.c_str(), SZ(locations) );
					for (int i = 0; i < SZ(locations); i++)
						fprintf(loc_out, "%lu ", (x.forward ? ref_len - locations[i] - k_length : locations[i] + k_length));
					fprintf(loc_out, "\n");
				}
			}
		}
	}
	cout << "LOCATIONS querying Time in modified procedure: " << (double)elapsed << "\n";
	fclose(loc_out);
}

void buildMinimizerIndex(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fm_ind_fwd, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> &fm_ind_rev)
{
	clock_gettime(CLOCK_MONOTONIC, &start);
	find_minimizers(REFF, retf, fm_ind_fwd);
	// cerr << "ENDED FORWARD" << endl;
	find_minimizers(REV_REFF, retr, fm_ind_rev);
	clock_gettime(CLOCK_MONOTONIC, &finish);

	double elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
	cout << "Time needed to compute all minimizers of Reference : ";
	cout << (double)elapsed << "\n";
	cout << "Number of Unique Minimizers:\n";
	cout << "Forward: " << SZ(retf) << "\nReverse: " << SZ(retr) << "\n";
}
/**
 * Command Lines:
 * Arguments:
 *      [0] - File name
 *      [1] - reference file path
 *      [2] - read file path
 *      [3] - K, length of k-mer
 *      [4] - W, window size
 *      [5] - output from enhanced fm index
 */
int main(int argc, const char **argv)
{
	//ios::sync_with_stdio(false);

	// ThreadPool pool(1);
	// vector< future<int> > results;
	double elapsed;

	string ref_filepath, read_filepath;
	string loc_enhanced; // to output fm index

	// Argument Handling: less argument supplied than expected
	if (argc <= 2)
	{
		cout << "Please! Give valid file path![K and W]\n";
		return 0;
	}
	if (argc >= 3)
	{
		ref_filepath = string(argv[1]);
		read_filepath = string(argv[2]);
	}
	if (argc >= 4)
	{
		K = atoi(argv[3]);
		min_k = K;
	}
	if (argc >= 5)
	{
		W = atoi(argv[4]);
	}
	if (argc >= 6) {
		loc_enhanced = string(argv[5]);
	}
	cout << ref_filepath << " " << read_filepath << " " << K << " " << W << " " << loc_enhanced << endl;
	ADK = 10;
	init();
	takeReference(ref_filepath);
	// Creating FM Index on Reference
	csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_ff, fm_rr;
	cout << " >>>>>>>>>>>>>>>>>>>> FM INDEX OUTPUT: BEGIN <<<<<<<<<<<<<<<<<<<<<<<<<\n";
	indexReferenceForEnhancedFM(fm_ff, fm_rr);
	cout << " >>>>>>>>>>>>>>>>>>>> FM INDEX OUTPUT:   END <<<<<<<<<<<<<<<<<<<<<<<<<\n";

	// Processing Reference
	buildMinimizerIndex(fm_ff, fm_rr);

	// Processing Reads
	vector<single_read>reads;

	takeReads(read_filepath, reads);

	processEnhancedFMAprroach(fm_ff, fm_rr, loc_enhanced, reads);

	return 0;
}
