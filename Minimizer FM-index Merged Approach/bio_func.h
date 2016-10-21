// modified in this version 20161005
long long pat2num(string s)
{
    long long ret = 0;
    int now = 0, kk = 0;
    for (int i = s.size() - 1; i >= 0; i--)
    {
        if (s[i] == '_') continue;
        kk = 0;
        if (s[i] == 'G' || s[i] == 'T') kk |= 2;
        if (s[i] == 'C' || s[i] == 'T') kk |= 1;
        ret |= (kk << now);
        now += 2;
    }
    return ret;
}

string num2pat(long long num, int k)
{
    string s = "";
    int tem = 0;
    for (int i = 0; i < k; i++)
    {
        tem = num & 3;
        if (tem == 0) s.push_back('A');
        else if (tem == 1) s.push_back('C');
        else if (tem == 2) s.push_back('G');
        else s.push_back('T');
        num >>= 2;
    }
    reverse(s.begin(), s.end());
    return s;
}

// added in this version 20161006
string num2pat2(long long num, int k)
{
    string s = "";
    int tem = 0;
    for (int i = 0; i < k; i++)
    {
        tem = num & 3;
        if (tem == 0) s.push_back('A');
        else if (tem == 1) s.push_back('C');
        else if (tem == 2) s.push_back('G');
        else s.push_back('T');
        if ((i % 3 == 1) && s.size() < k) s.push_back('_'), i++;
        num >>= 2;
    }
    reverse(s.begin(), s.end());
    return s;
}

int nt2num(char &x)
{
    if (x == 'A') return 0;
    else if (x == 'C') return 1;
    else if (x == 'G') return 2;
    return 3;
}

char num2nt(int &x)
{
    if (x & 1)
    {
        if (x & 2) return 'T';
        else return 'C';
    }
    else
    {
        if (x & 2) return 'G';
        else return 'A';
    }
}

//unordered_map<char, char> base;
long long rev_comp(long long x, int k)
{
    long long rev = 0;
    for (int i = 0; i < k; i++)
    {
        rev <<= 2;
        rev |= ( (x & 3) ^ 3);
        x >>= 2;
    }
    return rev;
}

// modified in this version 20161005
void insert_gap(string &s, const int init = 0) // parameter changed, optional parameter added to fix insert position
{
    int L = SZ(s);
    for (int i = ((2 + init) % 3); i < L; i += 3) // modified here
    {
        s[i] = '_'; // modified here
    }
    return ;
}

// modified from reverse_comp() in this version 20161005, added REV_REFF
void reverse_comp_of_ref()
{
    REV_REFF = REFF;
    reverse(REV_REFF.begin(), REV_REFF.end());

    int L = SZ(REV_REFF);
    for (int i = 0; i < L ; ++i)
    {
        //if (REV_REFF[i] == '_') REV_REFF[i] = '_';
        if (REV_REFF[i] == 'A')  REV_REFF[i] = 'T';
        else if (REV_REFF[i] == 'T')  REV_REFF[i] = 'A';
        else if (REV_REFF[i] == 'G')  REV_REFF[i] = 'C';
        else if (REV_REFF[i] == 'C')  REV_REFF[i] = 'G';

    }
}

string insert_gap_read(string str)
{
    int L = SZ(str);
    for (int i = 0; i < L; i++)
    {
        if ((i + 1) % 3 == 0 && i != 0) str[i] = '_';
    }
    return str;
}