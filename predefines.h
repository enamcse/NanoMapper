#include <sdsl/suffix_arrays.hpp>
#include <unordered_map>
#include <ctime>
#include <deque>
#include <cstdio>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <string>

using namespace sdsl;
using namespace std;

#include <zlib.h>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread);
const long long INF = 200000000000000000;

#define SZ(s)     ((int)s.size())
#define unq(vec)  stable_sort(vec.begin(),vec.end());\
                  vec.resize(distance(vec.begin(),unique(vec.begin(),vec.end())));
#define all(a)     a.begin(),a.end()

/////// These header files are included for Thread///////////////
// #include <chrono>                                              //
// #include "ThreadPool.h"                                        //
// #include <chrono>                                              //
/////////////////////////////////////////////////////////////////


int K = 14; ///minimizer length
int W = 24; ///minimizer windows size
int ADK = 10;
long long kmask;
int kmer_found, nxt_found;

///REV_REFF is the reverse complement of REFF -- reference. 20161005
string REFF, REV_REFF;

struct my_data
{
    long long minim;
    int idx;
    my_data(long long t, int p): minim(t), idx(p) {}
};

struct single_read{
	string seq;
	string name;
	int offs,strt,alignment;
	bool forward;
	single_read(kseq_t *seqs)
	{
		seq = string(seqs->seq.s);
        name = string(seqs->name.s);
        vector<string> strings;
        istringstream f(name);
        string s;
        while (getline(f, s, '_'))
        {
            strings.emplace_back(s);
        }
        offs = stoi(strings[5]);
        strt = stoi(strings[1]);
        alignment = (strings[2] != "unaligned");
        forward = (strings[4] == "F");
	}
};

unordered_map<long long, bool > retf, retr;
vector<int> for_index, rev_index;

#include "bio_func.h"
#include "minimizers.h"
#include "ModifiedFM.h"

void init()
{
    for (int i = 0; i < ADK; i++)
    {
        kmask <<= 2;
        kmask |= 3;
    }
}