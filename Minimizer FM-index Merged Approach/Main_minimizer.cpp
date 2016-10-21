/// Necessary header files, functions are included in "predefines.h"
#include "predefines.h"
/// OUTPUT_LOC - it would process FM-index navie approach
#define OUTPUT_LOC
/// OUTPUT_LOC1 - it would process FM-index enhanced approach
#define OUTPUT_LOC1
#undef OUTPUT_LOC

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
//        for(auto t : strings)
//        {
//            cout<<t<<" ";
//        }
//        cout<<endl;
		// results.emplace_back(pool.enqueue([file_cnter, NN, SS]
		// {
		file_cnter++;
		//return file_cnter;
		// }));
		//cout<<file_cnter<<"\n";
	}

	kseq_destroy(seq);
	gzclose(fp);
	fclose(ttp);
}
/**
 * fm_f - To index forward reference
 * fm_r - To index backward reference
 */
void indexReferenceForNaiveFM(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>&  fm_f, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_r)
{
	clock_gettime(CLOCK_MONOTONIC, &start);
	construct_im(fm_f, REFF, 1);
	construct_im(fm_r, REV_REFF, 1);
	clock_gettime(CLOCK_MONOTONIC, &finish);
	cout << " >>> INDEX CREATING <<<\n";
	cout << "\nIndexing of Reference (FOR NAIVE FM):\n";
	int elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
	cout << "Length of Reference: " << REV_REFF.size() << "\n";
	cout << "Reference Indexing Time: " << (double)elapsed << "\n";
	cout << "Memory Consumption in MB:\n";
	cout << "Forward = " << size_in_mega_bytes(fm_f) << "\n";
	cout << "Reverse = " << size_in_mega_bytes(fm_r) << "\n";
	cout << " >>> END - INDEX CREATING <<<\n";
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
	size_t k_length;
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
					if (retf.count(pat2num(temp_s)) == 0) continue;
					/// ATTENTION:  (1) first element should be omitted,
					///             (2) actual location could be found by subtructing from reference length
					k_length = 0;
					clock_gettime(CLOCK_MONOTONIC, &start);
					locations = my_locate(fm_f, kmer_with_gap[gap_counter].begin(), kmer_with_gap[gap_counter].end() - ii, min_k, max_count, k_length);
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
					if (retr.count(pat2num(temp_s)) == 0) continue;
					/// ATTENTION: See the note in the code of if section.
					k_length = 0;
					clock_gettime(CLOCK_MONOTONIC, &start);
					locations = my_locate(fm_r, kmer_with_gap[gap_counter].begin(), kmer_with_gap[gap_counter].end() - ii, min_k, max_count, k_length);
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

/**
 * fm_f - indexed forward reference
 * fm_r - indexed backward reference
 * location_out - file_path of output
 * reads - vector of reads
 */
void processNaiveFMAprroach(csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>&  fm_f, csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_r, string location_out, vector<single_read>&reads)
{
	FILE* loc_out = fopen((char*)location_out.c_str(), "w");
	string kmer_with_gap[3], temp_s;
	double elapsed = 0;
	int gap_counter, L, tcnt, tcntprev, st;
	size_t k_length;
	int_vector<64u>locations;
	loc_out = fopen((char*)location_out.c_str(), "w");
	/// Processing Time in normal procedure
	elapsed = 0;
	for (auto x : reads)
	{
		if (x.alignment)
		{
			L = x.seq.size();
			for (int temp_i = 0; temp_i < 3; temp_i++) {
				kmer_with_gap[temp_i] = x.seq;
				insert_gap(kmer_with_gap[temp_i], temp_i);
			}

			fprintf(loc_out, "> %s\n", x.name.c_str());

			for (int ii = 0, j; ii < L; ++ii)
			{
				gap_counter = ii % 3; // to maintain the index of kmer_with_gap
				st = ii;
				if (x.forward)
				{
					tcntprev = tcnt = 0;
					if (ii + min_k <= L)
					{
						temp_s = string(kmer_with_gap[gap_counter].begin() + ii, kmer_with_gap[gap_counter].begin() + min_k + ii);
					}
					else break;
					// checks whether the k-mer exists as minimizer in unordered_map
					if (retf.count(pat2num(temp_s)) == 0) continue;
					for (j = min_k; j + ii <= L; j++)
					{
						if (j != min_k)
						{
							temp_s.push_back(*(kmer_with_gap[gap_counter].begin() + j + ii - 1));
						}
						tcntprev = tcnt;
						clock_gettime(CLOCK_MONOTONIC, &start);
						tcnt = count(fm_f, temp_s.begin(), temp_s.end());
						clock_gettime(CLOCK_MONOTONIC, &finish);
						elapsed += (finish.tv_sec - start.tv_sec);
						elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
						if (tcnt == 0) break;
					}
					if (tcnt != 0) tcntprev = tcnt;
					else temp_s.pop_back();
					j--;
					if (tcntprev != 0 && tcntprev <= max_count)
					{
						clock_gettime(CLOCK_MONOTONIC, &start);

						locations = locate(fm_f, temp_s.begin(), temp_s.end());

						clock_gettime(CLOCK_MONOTONIC, &finish);
						elapsed += (finish.tv_sec - start.tv_sec);
						elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
						ii += j;
						ii--;
					}
				}
				else
				{
					tcntprev = tcnt = 0;
					if (ii + min_k <= L)
					{
						temp_s = string(kmer_with_gap[gap_counter].begin() + ii, kmer_with_gap[gap_counter].begin() + min_k + ii);
					}
					else break;
					// checks whether the k-mer exists as minimizer in unordered_map
					if (retr.count(pat2num(temp_s)) == 0) continue;
					for (j = min_k; j + ii <= L; j++)
					{
						if (j != min_k)
						{
							temp_s.push_back(*(kmer_with_gap[gap_counter].begin() + j + ii - 1));
						}
						tcntprev = tcnt;
						clock_gettime(CLOCK_MONOTONIC, &start);

						tcnt = count(fm_r, temp_s.begin(), temp_s.end());
						clock_gettime(CLOCK_MONOTONIC, &finish);
						elapsed += (finish.tv_sec - start.tv_sec);
						elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
						if (tcnt == 0) break;
					}
					// cerr << temp_s << endl;
					if (tcnt != 0) tcntprev = tcnt;
					else temp_s.pop_back();
					j--;
					if (tcntprev != 0 && tcntprev <= max_count)
					{
						clock_gettime(CLOCK_MONOTONIC, &start);
						locations = locate(fm_r, temp_s.begin(), temp_s.end());
						clock_gettime(CLOCK_MONOTONIC, &finish);
						elapsed += (finish.tv_sec - start.tv_sec);
						elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
						ii += j;
						ii--;
					}
				}

				if (tcntprev != 0 && tcntprev <= max_count)
				{
					fprintf(loc_out, "%d - %d : %s (%d)\n", st, ii, string(temp_s.begin(), temp_s.end()).c_str(), SZ(locations) );
					for (int i = 0; i < SZ(locations); i++)
						fprintf(loc_out, "%lu ", (x.forward ? locations[i] : ref_len - locations[i]));
					fprintf(loc_out, "\n");
				}
			}
		}
	}
	cout << "LOCATIONS querying Time in normal procedure: " << (double)elapsed << endl;
	fclose(loc_out);
}

/**
 * Command Lines:
 * Arguments:
 *      [0] - File name
 *      [1] - reference file path
 *      [2] - read file path
 *      [3] - K, length of k-mer
 *      [4] - W, window size
 *		[5] - RUnning mode, 1 - naive, 2 - enhanced, 3 both
 *      [6] - output from normal fm index(optional depanding on [5])
 *      [7] - output from enhanced fm index (optional depending on[5])
 */

void buildMinimizerIndex()
{
	clock_gettime(CLOCK_MONOTONIC, &start);
	find_minimizers(REFF, retf);
	find_minimizers(REV_REFF, retr);
	clock_gettime(CLOCK_MONOTONIC, &finish);

	double elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
	cout << "Time needed to compute all minimizers of Reference : ";
	cout << (double)elapsed << "\n";
	cout << "Number of Unique Minimizers:\n";
	cout << "Forward: " << SZ(retf) << "\nReverse: " << SZ(retr) << "\n";
}

int main(int argc, const char **argv)
{
	//ios::sync_with_stdio(false);

	// ThreadPool pool(1);
	// vector< future<int> > results;
	double elapsed, total;

	total = 0.0;
	string ref_filepath, read_filepath;
	string loc_naive, loc_enhanced; // to output fm index

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
	if (argc >= 6)
	{
		runningMode = atoi(argv[5]);
	}
	if (argc >= 6) {
		if (runningMode != 2)
			loc_naive = string(argv[6]);
		else loc_enhanced = string(argv[6]);
	}
	if (argc >= 7) {
		if (runningMode == 3)loc_enhanced = string(argv[7]);
	}
	cout << ref_filepath << " " << read_filepath << " " << K << " " << W << " " << loc_enhanced << endl;
	ADK = 10;
	init();

	// Processing Reference
	takeReference(ref_filepath);

	buildMinimizerIndex();

	// Creating FM Index on Reference
	csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_f, fm_r, fm_ff, fm_rr;
	cout << " >>>>>>>>>>>>>>>>>>>> FM INDEX OUTPUT: BEGIN <<<<<<<<<<<<<<<<<<<<<<<<<\n";
	if (runningMode & 1)indexReferenceForNaiveFM(fm_f, fm_r);
	if (runningMode & 2)indexReferenceForEnhancedFM(fm_ff, fm_rr);
	cout << " >>>>>>>>>>>>>>>>>>>> FM INDEX OUTPUT:   END <<<<<<<<<<<<<<<<<<<<<<<<<\n";

	// Processing Reads
	vector<single_read>reads;

	takeReads(read_filepath, reads);

	if (runningMode & 1)processNaiveFMAprroach(fm_f, fm_r, loc_naive, reads);
	if (runningMode & 2)processEnhancedFMAprroach(fm_ff, fm_rr, loc_enhanced, reads);

	return 0;
}
