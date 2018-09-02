/*
Build index for fasta file and save as binary files
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
    uint len = (bwt.length()-1)/r +1;
    //cout<<"BWT length: "<<bwt.length()<<'\t'<<"r: "<<r<<'\t'<<len<<endl;
    uint** Occ = new uint* [alphabet.length()];
    for(uint i=0; i<alphabet.length(); i++){
        Occ[i] = new uint [len];
        fill_n(Occ[i], len, 0);
        char_map[alphabet[i]]=i;
        char_occ[alphabet[i]]=0;
    }

    Occ[char_map[bwt[0]]][0] = 1;
    char_occ[bwt[0]] = 1;
    cout<<"alphabet: "<<alphabet<<endl;
    for(uint i=1; i<bwt.length(); i++){
        char_occ[bwt[i]] += 1;
        if(i%r==0){
        	//cout<<i<<endl;
            for(uint j=0; j<alphabet.length(); j++) Occ[j][i/r] = char_occ[alphabet[j]];
        }
    }
    //cout<<"Occ finished!"<<endl;
    return Occ;
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

void get_read_idx(string& bwt, uint* seq_sa, unordered_map<uint, uint>& ID_map, uint* seq_index_array){
    uint seq_len = bwt.length();
    uint dol_idx, arr_idx=0;
    for(uint i=0; i<seq_len; i++){
        if(bwt[i]=='$'){
            if(seq_sa[i]==0) dol_idx=seq_len-1;
            else dol_idx = seq_sa[i]-1;
            if(ID_map.find(dol_idx)==ID_map.end()) cout<<"Dollar Index error!"<<endl;
            seq_index_array[arr_idx] = ID_map[dol_idx];
            //cout<<ID_map[dol_idx]<<endl;
            arr_idx++;
        }
    }
    if(arr_idx!=ID_map.size()) cout<<"Read Index array error!"<<endl;
}

string encode_DNA(string& seq){
    int L = seq.length();
    int N=L%2==1?(L+1):L;
    string res;
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    for(int i=0; i<N; i+=2){
        char c1 = seq[i];
        char c2;
        if(i+1>=L) c2 = '&'; //append $ at the end of the bwt string if the length is odd
        else c2 = seq[i+1];
        char c=0;
        if(c1=='A'){ c=c|(1<<4);}
        else if(c1=='C'){ c|= 1<<5;}
        else if(c1=='G'){ c|= 1<<6;}
        else if(c1=='T'){ c|= 1<<7;}
        else if(c1=='N'){ c|= 1<<7; c|= 1<<5;}
        else if(c1=='$') ;
        else{ cout<<"Unknow character: "<<c1<<endl;}

        if(c2=='A'){ c=c|(1);}
        else if(c2=='C'){ c|= 1<<1;}
        else if(c2=='G'){ c|= 1<<2;}
        else if(c2=='T'){ c|= 1<<3;}
        else if(c2=='N'){ c|= 1<<3; c|= 1;}
        else if(c2=='$') ;
        else if(c2=='&'){ c|= 1<<3; c|= 1<<2; c|= 1<<1; c|= 1;} //use 1111 to represent sequence of the odd length
        else{ cout<<"Unknown character: "<<c2<<endl;}
        res+=c;
    }
    return res;
}

string decode_DNA(string& bit_seq, uint L, unordered_map<char,char>alphabet){
    string seq;
    for(uint i=0; i<bit_seq.length(); i++){
        char c = bit_seq[i];
        char c1 = alphabet[(c&240)>>4];
        char c2 = alphabet[(c&15)];
        seq+=c1;
        if(i<L/2){ seq+=c2;}
    }
    return seq;
}

void save_index(string& seq_whole, string& rev_seq_whole, unordered_map<uint,uint>& ID_map, uint r, string& out_file_base, uint idx){
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

    uint read_num = ID_map.size(); // reads number
    uint* seq_index_array = new uint [ID_map.size()];
    uint* rev_seq_index_array = new uint[ID_map.size()];

    string seq_bwt = get_bwt(seq_sa, seq_len, seq_whole);
    get_read_idx(seq_bwt, seq_sa, ID_map, seq_index_array);
    delete [] seq;
    delete [] seq_sa;
    string rev_seq_bwt = get_bwt(rev_seq_sa, seq_len, rev_seq_whole);
    get_read_idx(rev_seq_bwt, rev_seq_sa, ID_map, rev_seq_index_array);
    ID_map.clear();

    //cout<<rev_seq_whole<<endl;
    //cout<<rev_seq_bwt<<endl;
    delete [] rev_seq;
    delete [] rev_seq_sa;
    tuple<string, uint*> tup = get_C(seq_whole);
    string alphabet = get<0>(tup);
    uint* C = get<1>(tup);
    cout<<"C calculated"<<endl;
    uint** Occ = get_Occ(seq_bwt, alphabet, r);
    cout<<"OCC calculated"<<endl;
    uint seq_len_ratio = (seq_len-1)/r +1;
    //cout<<"The OCC is: "<<endl;
    //for(uint i=0; i<seq_len_ratio; i++) cout<<Occ[0][i]<<endl;
    uint** rev_Occ = get_Occ(rev_seq_bwt, alphabet, r);
    cout<<alphabet<<endl;

    // calculate the memory usage
    unsigned long long int memory_usage = 0; // bytes
    cout<<"The size of the sequence is:\t"<<sizeof(char)*seq_len<<endl;
    memory_usage += 4*sizeof(char)*seq_len; //+ bwt + rev_seq + rev_bwt
    cout<<"The size of suffix array is:\t"<<sizeof(uint)*seq_len<<endl;
    memory_usage += sizeof(uint)*seq_len;
    cout<<"The read index array is:\t"<<sizeof(uint)*read_num<<endl;
    memory_usage += 2*sizeof(uint)*read_num; //+ index array
    cout<<"The Occ array size is:\t"<<sizeof(uint)*alphabet.length()*seq_len/r<<endl;
    memory_usage += 2*sizeof(uint)*alphabet.length()*seq_len/r; // reverse OCC array
    cout<<"The size of suffix array, bwt, index array, and Occ array is:\t"<<memory_usage<<endl;
    cout<<memory_usage/1024.0/1024.0<<"MB\n";

    //cout<<"The total memory usage is:\t"<<2*memory_usage<<"bytes"<<"\t"<<2*memory_usage/1024/1024<<"MB"<<endl; // plus reverse

    // output the bwt and occ
    ofile0<<alphabet;
    string bit_bwt = encode_DNA(seq_bwt);
    string rev_bit_bwt = encode_DNA(rev_seq_bwt);
    ofile1<<bit_bwt<<rev_bit_bwt;
    for(uint i=0; i<alphabet.length(); i++) cout<<C[i]<<"\t";
    cout<<endl;
    ofile2.write((char*) C, alphabet.length()*sizeof(uint));
    for(uint i=0; i<alphabet.length(); i++) ofile2.write((char*) Occ[i], seq_len_ratio*sizeof(uint));
    for(uint i=0; i<alphabet.length(); i++) ofile2.write((char*) rev_Occ[i], seq_len_ratio*sizeof(uint));
    ofile3.write((char*) seq_index_array, read_num*sizeof(uint));
    ofile3.write((char*) rev_seq_index_array, read_num*sizeof(uint));
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
    uint d = 5; // number of partitions of the fasta file
    if (argc < 5){
        cout<<"Usage is -f <readfile> [-k <partition number (5)] -o <outfile>\n";
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
                else if(strcmp(argv[i], "-k")==0){ 
                    char* k_tmp = argv[i+1];
                    d = atoi(k_tmp);
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

    //uint r = 50;
    uint r = 100;
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
    uint seq_idx_all=0; //The total sequence length
    for(uint i=0; i<d; i++){
        vector<uint> seq_len_array;
        uint seq_idx = 0; //The total sequence length for this partition
        string seq_whole="";
        string rev_seq_whole = ""; //Concatenate the reverse complement sequences of the first fasta file
        string prefix = out_file_base+'_'+to_string(i);
        cout<<prefix<<endl;
        uint id_start=i*p;
        uint id_end = i==d-1?read_num:(i+1)*p;
        unordered_map<uint,uint> ID_map;
        for(uint j=id_start; j<id_end; j++){
            seq = readsData[j];
            seq_whole += seq +'$';
            seq_idx += seq.length() + 1;
            ID_map[seq_idx-1] = j;
            //cout<<j<<':'<<seq_idx<<" ";
            seq_idx_all += seq.length()+1;
            reverse(seq.begin(), seq.end());
            rev_seq_whole += seq+"$";
            seq_len_array.push_back(seq_idx);
        }
        //cout<<endl;
        cout<<seq_idx<<endl;
        cout<<seq_whole.substr(0, 20)<<endl;
        uint idx = i*p;
        save_index(seq_whole, rev_seq_whole, ID_map, r, prefix, idx);
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
