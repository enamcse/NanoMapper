void find_minimizers(string s, unordered_map<long long, bool > &ret) // reduced a parameter which was actually an index to count efficiently a minimizer in a range
{
    ret.clear();
    if (SZ(s) < K or SZ(s) < W)
    {
        if (SZ(s) < W) cerr << "Sequence Length is smaller than W after inserting gap(-).\n";
        else cerr << "Sequence Length is smaller than K after inserting gap(-).\n";
        return ;
    }
    int n = SZ(s);
    long long temps, minimizer, prevmin = -1;
    deque<my_data> sliding_window;
    temps = pat2num(s.substr(0, K));
    sliding_window.push_back(my_data(temps, 0));
    int i ;
    //cerr << num2pat(temps, ((K / 3) * 2)) << " " << temps << endl;
    for (i = K; i < n && i < W; ++i)
    {
        if (s[i] != '_')
        {
            temps <<= 2;
            temps &= kmask;
            temps |= nt2num(s[i]);
        }
        if (s[i - K] != '_') continue;
        minimizer = temps;
        //cerr << num2pat(temps, ((K / 3) * 2)) << "++" << num2pat(minimizer, ((K / 3) * 2)) << endl;
        while ((!sliding_window.empty()) && minimizer < sliding_window.back().minim)
            sliding_window.pop_back();
        sliding_window.push_back(my_data(minimizer, (i - K + 1)));
    }
    prevmin = sliding_window.front().minim;
    //cerr << i << "++" << endl;
    for (; i < n; i++)
    {
        if (s[i] != '_')
        {
            temps <<= 2;
            temps &= kmask;
            temps |= nt2num(s[i]);
        }
        if (s[i - K] != '_') continue;
        minimizer = temps;
        while ((!sliding_window.empty()) && sliding_window.front().idx <= (i - W))
        {

            sliding_window.pop_front();
        }
        while ((!sliding_window.empty()) && minimizer < sliding_window.back().minim)
            sliding_window.pop_back();
        sliding_window.push_back(my_data(minimizer, (i - K + 1)));

        if (sliding_window.front().minim != prevmin)
        {

            if (ret.count(prevmin) == 0) ret[ prevmin ] = true;
            prevmin = sliding_window.front().minim;
        }

    }


    if (ret.count(sliding_window.front().minim) == 0)
    {
        ret[ ret.count(sliding_window.front().minim) ] = true;
    }
    return ;
}

void find_kmers_of_read_in_refer(int flag, int offset, int start_in_ref, string read_name, string red_seq)
{
    string kkmer;
    int L = SZ(red_seq);
    int end_in_ref = start_in_ref + offset  + L, counter;
    start_in_ref += offset;
    int lidx, ridx;
    if (flag)
    {
        lidx = lower_bound(all(for_index), start_in_ref) - for_index.begin();
        ridx = upper_bound(all(for_index), end_in_ref) - for_index.begin();
        //cerr<<lidx<<" and "<<ridx<<" ++ "<<start_in_ref<<" -- "<<end_in_ref<<" :: "<<SZ(for_index)<<endl;
        if (lidx >= ridx) counter = 0;
        else counter = (ridx - lidx + 1);
    }
    else
    {
        lidx = lower_bound(all(rev_index), start_in_ref) - rev_index.begin();
        ridx = upper_bound(all(rev_index), end_in_ref) - rev_index.begin();
        //cerr<<lidx<<" or "<<ridx<<" ++ "<<start_in_ref<<" -- "<<end_in_ref<<" :: "<<SZ(rev_index)<<endl;
        if (lidx >= ridx) counter = 0;
        else counter = (ridx - lidx + 1);
    }
    //FILE* ffp = fopen((char*)read_name.c_str(), "w+");
    cout << read_name << endl;
    kmer_found = -1, nxt_found = 0;
    int mxm = 1;
    for (int ii = 0 ; ii <= L - K; ++ii)
    {
        kkmer = red_seq.substr(ii, K);
        kkmer = insert_gap_read(kkmer);
        long long kmer = pat2num(kkmer);
        //printf("%s %lld\n", num2pat(kmer, ADK).c_str(), kmer);
        if (flag)
        {
            if (retf.find(kmer) != retf.end())
            {
                nxt_found++;
                if (kmer_found == -1)
                {
                    kmer_found = ii;
                }
                else
                {
                    mxm = max(mxm, (ii - kmer_found - 1));
                    kmer_found = ii;
                }
//                for (auto t : retf[kmer])
//                {
//                    printf("%d ", t.fs);
//                }
//                printf("\n");
            }
//            else printf("Not found in reference.\n");
        }
        else
        {
            if (retr.find(kmer) != retr.end())
            {
                nxt_found++;
                if (kmer_found == -1)
                {
                    kmer_found = ii;
                }
                else
                {
                    mxm = max(mxm, (ii - kmer_found - 1));
                    kmer_found = ii;
                }
//                for (auto t : retr[kmer])
//                {
//                    printf("%d ", t.fs);
//                }
//                printf("\n");
            }
//            else printf("Not found in reference.\n");
        }
    }
    cout << start_in_ref << " " << end_in_ref << " " << counter << "\n";
    cout << nxt_found << "\n";

}