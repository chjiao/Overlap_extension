/*
Find overlap with saved bwt
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
    //unordered_map<char, uint> C;
    uint* C= new uint[alphabet.length()];
    for(uint i=0; i<alphabet.length(); i++){
        C[i] = count;
        //C[alphabet[i]] = count;
        count+= char_count[alphabet[i]];
    }
    return make_tuple(alphabet, C);
}

uint** get_Occ(string& bwt, string& alphabet){
    unordered_map<char, uint> char_map;
    uint** Occ = new uint* [alphabet.length()];
    for(uint i=0; i<alphabet.length(); i++){
        Occ[i] = new uint [bwt.length()];
        fill_n(Occ[i], bwt.length(), 0);
        char_map[alphabet[i]]=i;
    }

    Occ[char_map[bwt[0]]][0] = 1;
    for(uint i=1; i<bwt.length(); i++){
        for(uint j=0; j<alphabet.length(); j++) Occ[j][i] = Occ[j][i-1];
        Occ[char_map[bwt[i]]][i] = Occ[char_map[bwt[i]]][i-1] + 1;
    }
    return Occ;
}

void find_overlap(string& T, string& bwt, string& alphabet, uint cutoff, uint* C, uint** Occ, uint* seq_index_array, unordered_map<uint, uint>& saved_reads, vector<uint>& results, uint r){
    unordered_map<char, uint> char_map;
    for(uint i=0; i<alphabet.length(); i++) char_map[alphabet[i]]=i;

    uint i=T.length() -1;
    char c = T[i];
    char c1 = char_map[c]>=alphabet.length()-1?-1:alphabet[char_map[c]+1];
    uint b_start = C[char_map[c]];
    uint b_end = c1==-1?bwt.length()-1:C[char_map[c1]]-1;
    uint start1, end1;
    uint occ_start, occ1, occ_end, occ2;
    while(b_start<=b_end && i>=2){
        c=T[i-1];
        occ_start = (b_start -1)/r;
        occ1 = Occ[char_map[c]][occ_start];
        for(uint k=occ_start*r+1; k<=b_start-1; k++){
            if(bwt[k]==c) occ1++;
        }
        occ_end = b_end/r;
        occ2 = Occ[char_map[c]][occ_end];
        for(uint k=occ_end*r+1; k<=b_end; k++){
            if(bwt[k]==c) occ2++;
        }
        //cout<<"Updated occ: "<<occ1<<"\t"<<occ2<<endl;

        //b_start = C[char_map[c]]+Occ[char_map[c]][b_start-1];
        //b_end = C[char_map[c]]+Occ[char_map[c]][b_end]-1;
        b_start = C[char_map[c]]+occ1;
        b_end = C[char_map[c]]+occ2-1;
        //cout<<"Updated loc: "<<b_start<<"\t"<<b_end<<endl;
        i--;
        if(T.length()-i>=cutoff && b_start<=b_end){
            occ_start = (b_start -1)/r;
            occ1 = Occ[char_map['$']][occ_start];
            for(uint k=occ_start*r+1; k<=b_start-1; k++){
                if(bwt[k]=='$') occ1++;
            }
            occ_end = b_end/r;
            occ2 = Occ[char_map['$']][occ_end];
            for(uint k=occ_end*r+1; k<=b_end; k++){
                if(bwt[k]=='$') occ2++;
            }
            //cout<<"Updated occ ($): "<<occ1<<"\t"<<occ2<<endl;
            //start1 = C[char_map['$']]+Occ[char_map['$']][b_start-1];
            //end1 = C[char_map['$']]+Occ[char_map['$']][b_end]-1;
            start1 = C[char_map['$']]+occ1;
            end1 = C[char_map['$']]+occ2-1;
            if(start1<=end1){
                for(uint j=b_start; j<=b_end; j++){
                    if(bwt[j]=='$'){
                        if(saved_reads.find(seq_index_array[j])==saved_reads.end()){
                            saved_reads[seq_index_array[j]]=T.length()-i; // save the overlap size
                            results.push_back(seq_index_array[j]);
                        }
                    }
                }
            }
        }
    }
}

void rev_com(string&);
void find_all_overlap(unordered_map<uint, string>& seeds, string& bwt, string& rev_bwt, string& alphabet, uint cutoff, \
                      uint* C, uint** Occ, uint** rev_Occ, uint* seq_index_array, uint* rev_seq_index_array, unordered_map<uint, uint>& saved_reads, vector<uint>& results, uint r){
    for(auto it=seeds.begin(); it!=seeds.end(); it++){
        string temp = it->second;
        find_overlap(temp, bwt, alphabet, cutoff, C, Occ, seq_index_array, saved_reads, results, r);
        string rev_temp(temp);
        reverse(rev_temp.begin(), rev_temp.end());
        find_overlap(rev_temp, rev_bwt, alphabet, cutoff, C, rev_Occ, rev_seq_index_array, saved_reads, results, r);
        rev_com(temp);
        find_overlap(temp, bwt, alphabet, cutoff, C, Occ, seq_index_array, saved_reads, results, r);
        rev_temp = temp;
        reverse(rev_temp.begin(), rev_temp.end());
        find_overlap(rev_temp, rev_bwt, alphabet, cutoff, C, rev_Occ, rev_seq_index_array, saved_reads, results, r);
    }
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

unordered_map<uint, string> read_sam(char* sam_file, unordered_map<string, uint>& reads_map, vector<string>& readsData){
    unordered_map<uint, string> result;
    ifstream f(sam_file);
    string line;
    if(f.is_open()){
        while(getline(f, line)){
            if(line[0]=='@') continue;
            vector<string> lmap = split(line, '\t');
            if(stoi(lmap[1]) != 4){
                if(reads_map.find(lmap[0])!=reads_map.end()){
                    uint idx = reads_map[lmap[0]];
                    result[idx] = readsData[idx];
                }
                else cout<<"Read: "<<lmap[0]<<" not found."<<endl;
            }
        }
    }
    return result;
}

int main(int argc, char* argv[]){
    clock_t start_time=clock();
    char* sam_file;
    char* fa_file;
    char* index_file;
    uint cutoff = 0;
    char* out_file;
    if (argc < 11){
        cout<<"Usage is -S <samfile> -x <index> -f <readfile> -c <cutoff> -o <outfile>\n";
        cin.get();
        exit(0);
    }
    else{
        for (int i=1;i<argc;i++){
            if (i+1 <argc){
                if (strcmp(argv[i],"-S")==0){
                    sam_file = argv[i+1];
                    i++;
                }
                else if (strcmp(argv[i],"-x")==0){
                    index_file=argv[i+1];
                    i++;
                }
                else if (strcmp(argv[i],"-f")==0){
                    fa_file=argv[i+1];
                    i++;
                }
                else if(strcmp(argv[i], "-c")==0){
                    char* cut_tmp = argv[i+1];
                    cutoff = atoi(cut_tmp);
                    i++;
                }
                else if (strcmp(argv[i],"-o")==0){
                    out_file=argv[i+1];
                    i++;
                }
                else{
                    cout<<argv[i]<<endl;
                    cout<<"Not enough or invalid arguments, please try again.\n";
                    exit(0);
                }
            }
        }
    }

    string str_fa_file1(fa_file);
    //vector<uint> seq_len_array;
    //uint seq_idx = 0; //The total sequence length
    //string seq_whole="";
    //string rev_seq_whole = ""; //Concatenate the reverse complement sequences of the first fasta file
    vector<string> reads_title; // save all the reads title in a vector
    vector<string> readsData;
    unordered_map<string, uint> reads_map;
    uint read_num = 0;

    if (str_fa_file1.substr(str_fa_file1.length()-3,3)==string(".fa") || str_fa_file1.substr(str_fa_file1.length()-6,6)==string(".fasta")\
     || str_fa_file1.substr(str_fa_file1.length()-4)==string(".fna")){
        string line;
        string header=">";
        ifstream f (str_fa_file1.c_str());
        string title;
        string seq;

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
                    vector<string> lmap = split(line, ' ');
                    title=lmap[0].substr(1);
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

    // read built bwt
    uint r=50;
    uint d=5;
    vector<string> alphabet;
    //vector<uint*> C;
    uint** C = new uint* [d];
    vector<string> bwt;
    vector<string> rev_bwt;
    //vector<uint**> Occ;
    uint*** Occ = new uint** [d];
    //vector<uint**> rev_Occ;
    uint*** rev_Occ = new uint** [d];
    vector<uint*> seq_index_array;
    //uint** seq_index_array = new uint* [d];
    vector<uint*> rev_seq_index_array;
    //uint** rev_seq_index_array = new uint* [d];

    string index_base(index_file);
    for(uint i=0; i<d; i++){
        string in_alphabet = index_base+"_"+to_string(i)+".alphabet";
        string in_bwt = index_base+"_"+to_string(i)+".bwt";
        string in_occ = index_base+"_"+to_string(i)+".occ";
        string in_idx = index_base+"_"+to_string(i)+".idx";
        cout<<in_alphabet<<"\t"<<in_bwt<<"\t"<<in_occ<<"\t"<<in_idx<<endl;

        ifstream ofile0(in_alphabet.c_str());
        ifstream ofile1, ofile2, ofile3;

        ofile1.open(in_bwt.c_str());
        ofile2.open(in_occ.c_str(), ios::in | ios::binary);
        ofile3.open(in_idx.c_str(), ios::in | ios::binary);
        string p_alphabet;
        getline(ofile0, p_alphabet);
        alphabet.push_back(p_alphabet);

        string seq_bwt;
        getline(ofile1, seq_bwt);
        bwt.push_back(seq_bwt.substr(0, seq_bwt.length()/2));
        rev_bwt.push_back((seq_bwt.substr(seq_bwt.length()/2)));
        cout<<"bwt readed!"<<endl;
        uint seq_len = seq_bwt.length()/2;
        cout<<"seq_len: "<<seq_len<<endl;
        seq_bwt.clear();

        C[i] = new uint[p_alphabet.length()];
        cout<<C[i][1]<<endl;
        ofile2.read((char*) C[i], p_alphabet.length()*sizeof(uint));
        //C.push_back(p_C);
        Occ[i] = new uint* [p_alphabet.length()];
        rev_Occ[i] = new uint* [p_alphabet.length()];
        for(uint j=0; j<p_alphabet.length(); j++){
            Occ[i][j] = new uint [seq_len];
            ofile2.read((char*) Occ[i][j], seq_len*sizeof(uint)/r);
            cout<<Occ[i][j][1]<<endl;
            //Occ.push_back(p_Occ);
        }
        for(uint j=0; j<p_alphabet.length(); j++){
            rev_Occ[i][j] = new uint [seq_len];
            ofile2.read((char*) rev_Occ[i][j], seq_len*sizeof(uint)/r);
            cout<<rev_Occ[i][j][1]<<endl;
            //rev_Occ.push_back(p_rev_Occ);
        }
        cout<<"Occ readed!"<<endl;

        uint* index_array = new uint [seq_len];
        uint* rev_index_array = new uint [seq_len];
        ofile3.read((char*) index_array, seq_len*sizeof(uint));
        seq_index_array.push_back(index_array);
        ofile3.read((char*) rev_index_array, seq_len*sizeof(uint));
        rev_seq_index_array.push_back(rev_index_array);
        cout<<"Index readed!"<<endl;

        ofile0.close();
        ofile1.close();
        ofile2.close();
        ofile3.close();
    }
    cout<<"Index files readed!"<<endl;

    // overlap extension
    unordered_map<uint, string> seeds = read_sam(sam_file, reads_map, readsData);
    cout<<"The number of seed reads is: "<<seeds.size()<<endl;
    unordered_map<uint, uint> saved_reads;
    for(auto it=seeds.begin(); it!=seeds.end(); it++) saved_reads[it->first] = 1; // initialize with the seeds
    vector<uint> result;
    uint iter = 0;
    while(seeds.size()!=0){
        for(uint i=0; i<d; i++){
            find_all_overlap(seeds, bwt[i], rev_bwt[i], alphabet[i], cutoff, C[i], Occ[i], rev_Occ[i], seq_index_array[i], rev_seq_index_array[i], saved_reads, result, r);
        }
        cout<<"Iteration: "<<iter<<", recruited reads number: "<<result.size()<<endl;
        seeds.clear();
        for(uint i=0; i<result.size(); i++) seeds[result[i]] = readsData[result[i]]; // use the new recruited reads for next iteration
        result.clear();
        iter++;
    }
    cout<<"The total number of recruited reads (including seed reads) is: "<<saved_reads.size()<<endl;

    for(uint i=0; i<d; i++){
        delete [] seq_index_array[i];
        delete [] rev_seq_index_array[i];
        delete [] C[i];
        for(uint j=0; j<alphabet[i].length(); j++){
            delete [] Occ[i][j];
            delete [] rev_Occ[i][j];
        }
        delete [] Occ[i];
        delete [] rev_Occ[i];
    }
    delete [] C;
    delete [] Occ;
    delete rev_Occ;


    // output the results;
    ofstream ofile(out_file);
    for(auto it=saved_reads.begin(); it!=saved_reads.end(); it++){
        uint idx = it->first;
        ofile<<">"<<reads_title[idx]<<endl;
        ofile<<readsData[idx]<<endl;
    }
    ofile.close();

    clock_t now_time=clock();
    double elapsed_secs = double(now_time - start_time) / CLOCKS_PER_SEC;
    cout<<"Running time: "<<elapsed_secs<<endl;

    return 0;
}
