#include "predefines.h"
#define OUTPUT_LOC
#define OUTPUT_LOC1
#undef OUTPUT_LOC
/**
 * Command Lines:
 * Arguments:
 *      [0] - File name
 *      [1] - reference file path
 *      [2] - read file path
 *      [3] - K, length of k-mer
 *      [4] - W, window size
 *      [5] - output from normal fm index
 *      [6] - output from modified fm index
 */

int main(int argc, const char **argv)
{
    //ios::sync_with_stdio(false);
    //freopen("/Users/mkg/Desktop/Read Mapping Testing/gg_reads.txt", "r", stdin);
    //freopen("output.txt","w",stdout);
    size_t min_k = 14;/// 14
    size_t max_count = 20;
    bool show_string = false; /// shows the serached string after constructing
    bool debug = false; /// prints debuggings

    ADK = 10;
    init();
    // ThreadPool pool(1);
    // vector< future<int> > results;
    struct timespec start, finish;
    double elapsed, total;

    total = 0.0;
    string ref_filepath, read_filepath;
    string location_out, location_out1; // to output fm index

    // Argument Handling: less argument supplied than expected
    if (argc <= 2)
    {
        cout << "Please! Give valid file path![optional K and W]\n";
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
        location_out = string(argv[5]);
    }
    if (argc >= 7) {
        // cerr  << "It is defining!" << endl;
        location_out1 = string(argv[6]);
    }
    // Processing Reference
    string ref_name ;
    int ref_len;
    gzFile fp;
    kseq_t *seq;
    FILE* ttp = fopen((char*)ref_filepath.c_str(), "r");
    fp = gzdopen(fileno(ttp), "r");
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0)
    {
        REFF = (string)seq->seq.s;
        ref_name = (string)seq->name.s;
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
    clock_gettime(CLOCK_MONOTONIC, &start);
    find_minimizers(REFF, retf);
    find_minimizers(REV_REFF, retr);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
    total += elapsed;
    cout << "Time needed to compute all minimizers of Reference : ";
    cout << (double)elapsed << "\n";
    cout << "Number of Unique Minimizers:\n";
    cout << "Forward: " << SZ(retf) << "\nReverse: " << SZ(retr) << "\n";


    // Creating FM Index on Reference
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_f, fm_r, fm_ff, fm_rr;
    cout << " >>>>>>>>>>>>>>>>>>>> FM INDEX OUTPUT: BEGIN <<<<<<<<<<<<<<<<<<<<<<<<<\n";
#ifdef OUTPUT_LOC
    // cout << REFF << endl;
    clock_gettime(CLOCK_MONOTONIC, &start);
    construct_im(fm_ff, REFF, 1);
    construct_im(fm_rr, REV_REFF, 1);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    cout << " >>> INDEX CREATING <<<\n";

    cout << "\nIndexing of Reference:\n";


    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
    total += elapsed;
    cout << "Length of Reference: " << REV_REFF.size() << "\n";
    cout << "Reference Indexing Time: " << (double)elapsed << "\n";
    cout << "Memory Consumption in MB:\n";
    cout << "Forward = " << size_in_mega_bytes(fm_ff) << "\n";
    cout << "Reverse = " << size_in_mega_bytes(fm_rr) << "\n";
    cout << " >>> END - INDEX CREATING <<<\n";
#endif // OUTPUT_LOC


#ifdef OUTPUT_LOC1
    reverse(REFF.begin(), REFF.end());
    reverse(REV_REFF.begin(), REV_REFF.end());
    clock_gettime(CLOCK_MONOTONIC, &start);
    construct_im(fm_f, REFF, 1);
    construct_im(fm_r, REV_REFF, 1);
    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
    total += elapsed;

    cout << " >>> INDEX CREATING <<<\n";
    cout << "Indexing of Reverse Reference:\n";
    cout << "Length of Reference: " << REFF.size() << "\n";
    cout << "Reference Indexing Time: " << (double)elapsed << "\n";
    cout << "Memory Consumption in MB:\n";
    cout << "Forward = " << size_in_mega_bytes(fm_f) << "\n";
    cout << "Reverse = " << size_in_mega_bytes(fm_r) << "\n";
    cout << " >>> END - INDEX CREATING <<<\n";
#endif //OUTPUT_LOC1
    cout << " >>>>>>>>>>>>>>>>>>>> FM INDEX OUTPUT:   END <<<<<<<<<<<<<<<<<<<<<<<<<\n";


    // Processing Reads
    ttp = fopen((char*)read_filepath.c_str(), "r");
    fp = gzdopen(fileno(ttp), "r");
    seq = kseq_init(fp);
    int file_cnter = 0;
    vector<single_read>reads;

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

    cout << "Total Reads = " << file_cnter << endl;
    string temp_s;
    int gap_counter = 0;
    int_vector<64u>locations;
    size_t k_length;
    string kkmer;
    int cnt = 0, st, L, tcnt, tcntprev;
    FILE* loc_out, *loc_out1;

    string kmer_with_gap[3];
#ifdef OUTPUT_LOC
    loc_out = fopen((char*)location_out.c_str(), "w");
    /// Processing Time in normal procedure
    elapsed = 0;
    for (auto x : reads)
    {
        cnt++;
        if (x.alignment)
        {
            L = x.seq.size();
            for (int temp_i = 0; temp_i < 3; temp_i++) {
                kmer_with_gap[temp_i] = x.seq;
                insert_gap(kmer_with_gap[temp_i], temp_i);
                // cerr << kmer_with_gap[temp_i] << endl;
                // cerr << temp_i << ": " << kmer_with_gap[temp_i] << endl;
            }

            fprintf(loc_out, "> %s\n", x.name.c_str());

            for (int ii = 0, j; ii < L; ++ii)
            {
                gap_counter = ii % 3; // to maintain the index of kmer_with_gap
                st = ii;
                if (x.forward)
                {
                    /// ATTENTION:  (1) first element should be omitted,
                    tcntprev = tcnt = 0;
                    if (ii + min_k <= L)
                    {
                        temp_s = string(kmer_with_gap[gap_counter].begin() + ii, kmer_with_gap[gap_counter].begin() + min_k + ii);
                        // cerr << ii << ": " << temp_s << endl;
                    }
                    else break;
                    // checks whether the k-mer exists as minimizer in unordered_map
                    if (retf.count(pat2num(temp_s)) == 0) continue;
                    for (j = min_k; j + ii <= L; j++)
                    {
                        // cerr << "j = " << j << endl;
                        if (j != min_k)
                        {
                            temp_s.push_back(*(kmer_with_gap[gap_counter].begin() + j + ii - 1));
                        }
                        // cerr << "temp_s = " << temp_s << endl;
                        tcntprev = tcnt;
                        clock_gettime(CLOCK_MONOTONIC, &start);
                        tcnt = count(fm_ff, temp_s.begin(), temp_s.end());
                        clock_gettime(CLOCK_MONOTONIC, &finish);
                        elapsed += (finish.tv_sec - start.tv_sec);
                        elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
                        // cerr << "tcnt = " << tcnt << ", tcntprev = " << tcntprev << endl;
                        if (tcnt == 0) break;
                    }
                    if (tcnt != 0) tcntprev = tcnt;
                    else temp_s.pop_back();
                    // cerr << "Final temp_s = " << temp_s << endl;
                    j--;
                    // cerr << tcnt << " <ooo> " << tcntprev << " >> M = " << max_count << endl;
                    // cerr << j << " L = " << L << endl;
                    if (tcntprev != 0 && tcntprev <= max_count)
                    {
                        // locations = locate(fm_ff, kkmer.begin()+ii, kkmer.end());
                        clock_gettime(CLOCK_MONOTONIC, &start);

                        locations = locate(fm_ff, temp_s.begin(), temp_s.end());

                        clock_gettime(CLOCK_MONOTONIC, &finish);
                        elapsed += (finish.tv_sec - start.tv_sec);
                        elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
                        ii += j;
                        ii--;
                    }
                    // freq_f[occs]++;
                }
                else
                {
                    /// ATTENTION: See the note in the code of if section.

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

                        tcnt = count(fm_rr, temp_s.begin(), temp_s.end());
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
                        // locations = locate(fm_ff, kkmer.begin()+ii, kkmer.end());
                        clock_gettime(CLOCK_MONOTONIC, &start);
                        locations = locate(fm_rr, temp_s.begin(), temp_s.end());
                        clock_gettime(CLOCK_MONOTONIC, &finish);
                        elapsed += (finish.tv_sec - start.tv_sec);
                        elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
                        ii += j;
                        ii--;
                    }
                    // freq_r[occs]++;
                }

                if (tcntprev != 0 && tcntprev <= max_count)
                {
                    fprintf(loc_out, "%d - %d : %s (%d)\n", st, ii, string(temp_s.begin(), temp_s.end()).c_str(), SZ(locations) );
                    for (int i = 0; i < SZ(locations); i++)
                        fprintf(loc_out, "%lu ", (x.forward ? locations[i] : ref_len - locations[i]));
                    fprintf(loc_out, "\n");
                    /// print locations here
                }
            }
        }
        if (cnt % 1 == 0) cerr << cnt << endl;
    }
    cout << "LOCATIONS querying Time in normal procedure: " << (double)elapsed << endl;
    fclose(loc_out);
#endif
    cnt = 0;
    elapsed = 0;

#ifdef OUTPUT_LOC1
    loc_out1 = fopen((char*)location_out1.c_str(), "w");
    // string kmer_with_gap[3];

    /// Processing time in modified procedure
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (auto x : reads)
    {
        cnt++;
        if (x.alignment)
        {
            L = x.seq.size();
            for (int temp_i = 0; temp_i < 3; temp_i++) {
                kmer_with_gap[temp_i] = x.seq;
                insert_gap(kmer_with_gap[temp_i], temp_i);
                reverse(kmer_with_gap[temp_i].begin(), kmer_with_gap[temp_i].end());
            }

            // cerr << x.seq << " >> " << kkmer << endl;
            // cout<<"Actual Location = "<<x.strt<<"\n";
            fprintf(loc_out1, "> %s\n", x.name.c_str());
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
                    locations = my_locate(fm_f, kmer_with_gap[gap_counter].begin(), kmer_with_gap[gap_counter].end() - ii, min_k, max_count, k_length);
                    // clock_gettime(CLOCK_MONOTONIC, &finish);
                    // elapsed += (finish.tv_sec - start.tv_sec);
                    // elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
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
                    // clock_gettime(CLOCK_MONOTONIC, &start);
                    locations = my_locate(fm_r, kmer_with_gap[gap_counter].begin(), kmer_with_gap[gap_counter].end() - ii, min_k, max_count, k_length);

                    if (k_length >= min_k and SZ(locations) != 0)
                    {
                        ii += k_length;
                        --ii;
                    }
                    // freq_r[occs]++;
                }
                if (k_length >= min_k and SZ(locations) != 0)
                {
                    temp_s = string(kmer_with_gap[gap_counter].end() - ii - 1, kmer_with_gap[gap_counter].end() - st);
                    // insert_gap(temp_s);
                    reverse(temp_s.begin(), temp_s.end());
                    fprintf(loc_out1, "%d - %d : %s (%d)\n", st, ii, temp_s.c_str(), SZ(locations) );
                    for (int i = 0; i < SZ(locations); i++)
                        fprintf(loc_out1, "%lu ", (x.forward ? ref_len - locations[i] - k_length : locations[i] + k_length));
                    fprintf(loc_out1, "\n");
                    /// print locations here
                    // if (show_string)
                    // {
                    //     temp_s = string(kkmer.end() - ii - 1, kkmer.end() - st);
                    //     reverse(temp_s.begin(), temp_s.end());
                    //     cout << "Kmer = " << temp_s << ", Count = " << SZ(locations) << endl;
                    // }

                    // if (debug)cout << "Locations:\n";
                    // for (int i = 0; i < SZ(locations); i++)
                    // {
                    //     if (debug)cout << ref_len - locations[i] - k_length << endl;
                    // }
                }
            }
        }
        if (cnt % 100 == 0) cerr << cnt << endl;
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed += (finish.tv_sec - start.tv_sec);
    elapsed += ((finish.tv_nsec - start.tv_nsec) * 1.0 / 1000000000.0);
    cout << "LOCATIONS querying Time in modified procedure: " << (double)elapsed << "\n";
    fclose(loc_out1);
#endif

    return 0;
}
