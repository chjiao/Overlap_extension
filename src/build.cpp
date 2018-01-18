/*
Build index for fasta file and save as binary files
To do:
save as 3 files, (1) alphabet, bwt and rev_bwt (2) C, Occ and rev_Occ (3) index_array and rev_index_array
*/

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <dirent.h>
#include <errno.h>
#include <ctime>
#include <sstream>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include "radix.h"

using namespace std;
typedef unsigned int uint;

void upper_str(string &seq){
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

int getdir (const string& dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        //cout << "Error opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if (string(dirp->d_name)!=string(".") && string(dirp->d_name)!=string("..")){
            files.push_back(dir+string(dirp->d_name));
        }
    }
    closedir(dp);
    return 0;
}

int max_array(int *A,int len){
    int max_num=A[0];
    for (int i=0;i<len;i++){
        if (max_num<A[i]) max_num=A[i];
    }
    return max_num;
}

int min_array(int *A,int len){
    int min_num=A[0];
    for (int i=0;i<len;i++){
        if (min_num>A[i]) min_num=A[i];
    }
    return min_num;
}

vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

string get_bwt(uint* sa, uint seq_len, string& seq){
    string bwt;
    unordered_map<char, int> char_count;
    for(uint i=0; i<seq_len; i++){
        if(char_count.find(seq[i])==char_count.end()) char_count[seq[i]]=1;
        if(sa[i]==0) bwt+="$";
        else bwt+=seq[sa[i]-1];
    }
    return bwt;
}

tuple<string, uint*> get_C(string &seq){
    unordered_map<char, uint> char_count;
    for(uint i=0; i<seq.length(); i++){
        if(char_count.find(seq[i])==char_count.end()) char_count[seq[i]]=1;
        else char_count[seq[i]]+=1;
    }
    string alphabet;
    for(auto it=char_count.begin(); it!=char_count.end();it++) alphabet += it->first;
    sort(alphabet.begin(), alphabet.end());

    uint count=0;
    uint* C= new uint[alphabet.length()];
    for(uint i=0; i<alphabet.length(); i++){
        C[i] = count;
        count+= char_count[alphabet[i]];
    }
    return make_tuple(alphabet, C);
}

uint** get_Occ(string& bwt, string& alphabet, uint r=50){
    unordered_map<char, uint> char_map;
    unordered_map<char, uint> char_occ;
    uint len = bwt.length()/r;
    uint** Occ = new uint* [alphabet.length()];
    for(uint i=0; i<alphabet.length(); i++){
        Occ[i] = new uint [len];
        fill_n(Occ[i], len, 0);
        char_map[alphabet[i]]=i;
        char_occ[alphabet[i]]=0;
    }

    Occ[char_map[bwt[0]]][0] = 1;
    char_occ[bwt[0]] = 1;
    for(uint i=1; i<bwt.length(); i++){
        char_occ[bwt[i]] += 1;
        if(i%r==0){
            for(uint j=0; j<alphabet.length(); j++) Occ[j][i/r] = char_occ[alphabet[j]];
        }
    }
    return Occ;
}

unsigned int get_suffix_idx(unsigned int sa_val, unsigned int left, unsigned int right, vector<unsigned int> &seq_len_array){
    //seq_len_array: an array storing the sequence lengths
    uint mid = (left+right)/2;
    if(mid==left){
        if(sa_val>=seq_len_array[left]) return right;
        else return left;
    }
    if(sa_val<seq_len_array[mid]){
        right = mid;
        return get_suffix_idx(sa_val, left, right, seq_len_array);
    }
    else{
        left = mid;
        return get_suffix_idx(sa_val, left, right, seq_len_array);
    }
}

unsigned int* get_seq_index(unsigned int* sa, unsigned int sa_len, vector<unsigned int> &seq_len_array){
    unsigned int* result = new unsigned int [sa_len];
    for(unsigned int i=0; i<sa_len; i++){
        unsigned int index = get_suffix_idx(sa[i], 0, seq_len_array.size()-1, seq_len_array);
        result[i] = index;
    }
    return result;
}

char seq_alphabet[] = "ATUGCYRSWKMBDHVN";
char rev_alphabet[] = "TAACGRYSWMKVHDBN";
unordered_map<char, char> convert ({{'A','T'},{'T','A'},{'U','A'},{'C','G'},{'G','C'},{'Y','R'},{'R','Y'},{'S','S'}, \
{'W','W'},{'K','M'},{'M','K'},{'B','V'},{'D','H'},{'H','D'},{'V','B'},{'N','N'}});

char complement(char n)
{
    if(convert.find(n)!=convert.end()) return convert[n];
    else return ' ';
}

void rev_com(string &seq){
    transform(seq.begin(), seq.end(), seq.begin(), complement);
    reverse(seq.begin(), seq.end());
}

void save_index(string& seq_whole, string& rev_seq_whole, vector<uint>& seq_len_array, uint r, string& out_file_base, uint idx){
    ofstream ofile0, ofile1, ofile2, ofile3;
    string out_alphabet = out_file_base+".alphabet";
    string out_bwt = out_file_base+".bwt";
    string out_occ = out_file_base+".occ";
    string out_idx = out_file_base+".idx";

    ofile0.open(out_alphabet.c_str(), ios::out | ios::binary);
    ofile1.open(out_bwt.c_str(), ios::out | ios::binary);
    ofile2.open(out_occ.c_str(), ios::out | ios::binary);
    ofile3.open(out_idx.c_str(), ios::out | ios::binary);

    unsigned char *seq =new unsigned char[seq_whole.length()+1]; //'\0' as the string end
    unsigned char *rev_seq = new unsigned char[rev_seq_whole.length()+1];
    strcpy((char *)seq, seq_whole.c_str());
    strcpy((char *)rev_seq, rev_seq_whole.c_str());
    uint seq_len = seq_whole.length();  // total sequence length
    unsigned int *seq_sa = Radix(seq, seq_len).build(); // suffix array
    uint *rev_seq_sa = Radix(rev_seq, seq_len).build();
    unsigned int *seq_index_array = get_seq_index(seq_sa, seq_len, seq_len_array); // original sequence index

    // use consistent sequence ids for different partitions
    if(idx>0){
        for(uint i=0; i<seq_len; i++) seq_index_array[i] += idx;
    }
    uint *rev_seq_index_array = get_seq_index(rev_seq_sa, seq_len, seq_len_array);
    if(idx>0){
        for(uint i=0; i<seq_len; i++) rev_seq_index_array[i] += idx;
    }
    string seq_bwt = get_bwt(seq_sa, seq_len, seq_whole);
    delete [] seq;
    delete [] seq_sa;
    string rev_seq_bwt = get_bwt(rev_seq_sa, seq_len, rev_seq_whole);
    delete [] rev_seq;
    delete [] rev_seq_sa;
    tuple<string, uint*> tup = get_C(seq_whole);
    string alphabet = get<0>(tup);
    uint* C = get<1>(tup);
    uint** Occ = get_Occ(seq_bwt, alphabet, r);
    uint** rev_Occ = get_Occ(rev_seq_bwt, alphabet, r);

    // calculate the memory usage
    unsigned long long int memory_usage = 0; // bytes
    cout<<"The size of the sequence is:\t"<<sizeof(char)*seq_len<<endl;
    memory_usage += 2*sizeof(char)*seq_len; //+ bwt
    cout<<"The size of suffix array is:\t"<<sizeof(uint)*seq_len<<endl;
    memory_usage += 2*sizeof(uint)*seq_len; //+ index array
    cout<<alphabet<<endl;
    cout<<"The Occ array size is:\t"<<sizeof(uint)*alphabet.length()*seq_len<<endl;
    memory_usage += sizeof(uint)*alphabet.length()*seq_len;
    cout<<"The size of suffix array, bwt, index array, and Occ array is:\t"<<memory_usage<<endl;
    cout<<memory_usage/1024/1024<<"MB\n";

    cout<<"The total memory usage is:\t"<<2*memory_usage<<"bytes"<<"\t"<<2*memory_usage/1024/1024<<"MB"<<endl; // plus reverse

    // output the bwt and occ
    ofile0<<alphabet;
    ofile1<<seq_bwt<<rev_seq_bwt;
    for(uint i=0; i<alphabet.length(); i++) cout<<C[i]<<"\t";
    cout<<endl;
    ofile2.write((char*) C, alphabet.length()*sizeof(uint));
    for(uint i=0; i<alphabet.length(); i++) ofile2.write((char*) Occ[i], seq_len*sizeof(uint)/r);
    for(uint i=0; i<alphabet.length(); i++) ofile2.write((char*) rev_Occ[i], seq_len*sizeof(uint)/r);
    ofile3.write((char*) seq_index_array, seq_len*sizeof(uint));
    ofile3.write((char*) rev_seq_index_array, seq_len*sizeof(uint));
    ofile0.close();
    ofile1.close();
    ofile2.close();
    ofile3.close();

    // release the memory
    delete [] seq_index_array;
    delete [] rev_seq_index_array;
    delete [] C;
    for(uint i=0; i<alphabet.length(); i++){
        delete [] Occ[i];
        delete [] rev_Occ[i];
    }
    delete [] Occ;
    delete [] rev_Occ;
}

int main(int argc, char* argv[]){
    clock_t start_time=clock();
    char* fa_file;
    char* out_file;
    if (argc < 5){
        cout<<"Usage is -f <readfile> -o <outfile>\n";
        cin.get();
        exit(0);
    }
    else{
        for (int i=1;i<argc;i++){
            if (i+1 <argc){
                if (strcmp(argv[i],"-f")==0){
                    fa_file=argv[i+1];
                    i++;
                }
                else if (strcmp(argv[i],"-o")==0){
                    out_file=argv[i+1];
                    i++;
                }
                else{
                    cout<<argv[i]<<endl;
                    cout<<"Not enough or invalid arguments, please try again.\n";
                    return 0;
                }
            }
        }
    }

    uint r = 50;
    uint d = 5; // number of partitions of the fasta file
    string out_file_base(out_file);
    string str_fa_file1(fa_file);
    vector<string> reads_title; // save all the reads title in a vector
    vector<string> readsData;
    unordered_map<string, uint> reads_map;
    uint read_num = 0;
    string title;
    string seq;

    if (str_fa_file1.substr(str_fa_file1.length()-3,3)==string(".fa") || str_fa_file1.substr(str_fa_file1.length()-6,6)==string(".fasta")\
     || str_fa_file1.substr(str_fa_file1.length()-4)==string(".fna")){
        string line;
        string header=">";
        ifstream f (str_fa_file1.c_str());

        if (f.is_open()){
            while(getline (f,line))
            {
                if (line[0]=='>'){
                    if(!title.empty() and !seq.empty()){
                        upper_str(seq);   // convert to upper case
                        reads_map[title] = read_num;
                        readsData.push_back(seq);
                        //seq_whole += seq + "$";
                        //seq_idx += seq.length() + 1;
                        //reverse(seq.begin(), seq.end());
                        //rev_seq_whole += seq+"$";
                        //seq_len_array.push_back(seq_idx);
                        reads_title.push_back(title);
                        read_num++;
                    }
                    title=line;
                    seq.clear();
                }
                else{
                    seq+=line;
                }
            }
            // the last sequence
            upper_str(seq);
            reads_map[title] = read_num;
            readsData.push_back(seq);
            //seq_whole += seq+"$";
            //seq_idx += seq.length() + 1;
            //reverse(seq.begin(), seq.end());
            //rev_seq_whole += seq+"$";
            //seq_len_array.push_back(seq_idx);
            reads_title.push_back(title);
            read_num++;
            f.close();

            cout<<"The total number of reads number is:　"<<read_num<<endl;
        }
        else cout <<"Unable to open file\n";
    }
    else return 0;

    uint p=read_num/d;
    cout<<read_num<<endl;
    uint seq_idx_all=0;
    for(uint i=0; i<d; i++){
        vector<uint> seq_len_array;
        uint seq_idx = 0; //The total sequence length
        string seq_whole="";
        string rev_seq_whole = ""; //Concatenate the reverse complement sequences of the first fasta file
        string prefix = out_file_base+'_'+to_string(i);
        cout<<prefix<<endl;
        uint id_start=i*p;
        uint id_end = i==d-1?read_num:(i+1)*p;
        for(uint j=id_start; j<id_end; j++){
            seq = readsData[j];
            seq_whole += seq +'$';
            seq_idx += seq.length() + 1;
            //cout<<j<<':'<<seq_idx<<" ";
            seq_idx_all += seq.length()+1;
            reverse(seq.begin(), seq.end());
            rev_seq_whole += seq+"$";
            seq_len_array.push_back(seq_idx);
        }
        //cout<<endl;
        cout<<seq_idx<<endl;
        cout<<seq_whole.substr(0, 20)<<endl;
        uint idx = i==0?0:i*p;
        save_index(seq_whole, rev_seq_whole, seq_len_array, r, prefix, idx);
        seq_len_array.clear();
        seq_whole.clear();
        rev_seq_whole.clear();
    }
    cout<<"The total length of the fa file is:\t"<<seq_idx_all<<endl;

    clock_t now_time=clock();
    double elapsed_secs = double(now_time - start_time) / CLOCKS_PER_SEC;
    cout<<"Running time: "<<elapsed_secs<<endl;

    return 0;
}
