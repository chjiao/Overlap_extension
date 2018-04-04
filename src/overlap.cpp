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

void find_overlap(string& T, string& bwt, string& alphabet, uint cutoff, uint* C, uint** Occ, uint* seq_index_array, \
unordered_map<uint, uint>& saved_reads, vector<uint>& results, uint r, uint min_read_idx, uint max_read_idx){
    unordered_map<char, uint> char_map;
    for(uint i=0; i<alphabet.length(); i++) char_map[alphabet[i]]=i;

    uint i=T.length() -1;
    char c = T[i];
    if(char_map.find(c)==char_map.end()) return;
    char c1 = char_map[c]>=alphabet.length()-1?-1:alphabet[char_map[c]+1]; // the next character
    uint b_start = C[char_map[c]];
    uint b_end = c1==-1?bwt.length()-1:C[char_map[c1]]-1;
    //cout<<"i: "<<i<<T[i]<<'\t'<<"Initial loc: "<<b_start<<"\t"<<b_end<<endl;
    //cout<<b_start<<'\t'<<b_end<<endl;
    uint start1, end1;
    uint occ_start, occ1, occ_end, occ2;
    while(b_start<=b_end && i>=2){
        c=T[i-1];
        if(char_map.find(c)==char_map.end()) return;
        occ_start = (b_start -1)/r;
        occ1 = Occ[char_map[c]][occ_start];
        occ_end = b_end/r;
        occ2 = Occ[char_map[c]][occ_end];
        //cout<<"Initial occ: "<<occ_start<<'\t'<<occ_end<<'\t'<<occ1<<"\t"<<occ2<<endl;
        for(uint k=occ_start*r+1; k<=b_start-1; k++){
            if(bwt[k]==c) occ1++;
        }
        for(uint k=occ_end*r+1; k<=b_end; k++){
            if(bwt[k]==c) occ2++;
        }
        //cout<<"Updated occ: "<<occ1<<"\t"<<occ2<<endl;

        b_start = C[char_map[c]]+occ1;
        b_end = C[char_map[c]]+occ2-1;
        //cout<<"i: "<<i-1<<T[i-1]<<'\t'<<"Updated loc: "<<b_start<<"\t"<<b_end<<endl;
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
            if(occ2<=occ1) continue;
            start1 = C[char_map['$']]+occ1;
            end1 = C[char_map['$']]+occ2-1;
            //cout<<"occ1: "<<occ1<<endl;
            //cout<<"occ2: "<<occ2<<endl;
            if(start1<=end1){
                //cout<<start1<<'\t'<<end1<<endl;
                for(uint j=start1; j<=end1; j++){
                    uint read_id = seq_index_array[j]==max_read_idx?min_read_idx:seq_index_array[j]+1;
                    //cout<<"Read ID is: "<<read_id<<endl;
                    if(saved_reads.find(read_id)==saved_reads.end()){
                        saved_reads[read_id]=T.length()-i; // save the overlap size
                        results.push_back(read_id);
                    }
                }
            }
        }
    }
}

void rev_com(string&);
void find_all_overlap(unordered_map<uint, string>& seeds, string& bwt, string& rev_bwt, string& alphabet, uint cutoff, \
                      uint* C, uint** Occ, uint** rev_Occ, uint* seq_index_array, uint* rev_seq_index_array, unordered_map<uint, uint>& saved_reads, vector<uint>& results, uint r){

    uint max_idx = 0, min_idx=0, rev_max_idx=0, rev_min_idx=0;
    uint read_num=C[1];
    for(uint i=0; i<read_num; i++){
        if(seq_index_array[i]>max_idx) max_idx=seq_index_array[i];
        if(i==0) min_idx=seq_index_array[0];
        else if (seq_index_array[i]<min_idx) min_idx=seq_index_array[i];

        if(rev_seq_index_array[i]>rev_max_idx) rev_max_idx=rev_seq_index_array[i];
        if(i==0) rev_min_idx = rev_seq_index_array[0];
        else if (rev_seq_index_array[i]<rev_min_idx) rev_min_idx=rev_seq_index_array[i];
    }

    for(auto it=seeds.begin(); it!=seeds.end(); it++){
        string temp = it->second;
        find_overlap(temp, bwt, alphabet, cutoff, C, Occ, seq_index_array, saved_reads, results, r, min_idx, max_idx);
        string rev_temp(temp);
        reverse(rev_temp.begin(), rev_temp.end());
        find_overlap(rev_temp, rev_bwt, alphabet, cutoff, C, rev_Occ, rev_seq_index_array, saved_reads, results, r, rev_min_idx, rev_max_idx);
        rev_com(temp);
        find_overlap(temp, bwt, alphabet, cutoff, C, Occ, seq_index_array, saved_reads, results, r, min_idx, max_idx);
        rev_temp = temp;
        reverse(rev_temp.begin(), rev_temp.end());
        find_overlap(rev_temp, rev_bwt, alphabet, cutoff, C, rev_Occ, rev_seq_index_array, saved_reads, results, r, rev_min_idx, rev_max_idx);
    }
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

void output_fasta(char* outfile, unordered_map<uint, string> seeds, vector<string> reads_title){
    ofstream ofile(outfile);
    for(auto it=seeds.begin(); it!=seeds.end(); it++){
        uint idx = it->first;
        ofile<<">"<<reads_title[idx]<<endl;
        ofile<<it->second<<endl;
    }
    ofile.close();
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
    uint** C = new uint* [d];
    vector<string> bwt;
    vector<string> rev_bwt;
    uint*** Occ = new uint** [d];
    uint*** rev_Occ = new uint** [d];
    vector<uint*> seq_index_array;
    vector<uint*> rev_seq_index_array;

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
        ofile2.read((char*) C[i], p_alphabet.length()*sizeof(uint));
        cout<<C[i][1]<<endl;
        uint p_read_num = C[i][1];
        //C.push_back(p_C);
        Occ[i] = new uint* [p_alphabet.length()];
        rev_Occ[i] = new uint* [p_alphabet.length()];
        uint seq_len_ratio = ((seq_len-1)/r+1);
        for(uint j=0; j<p_alphabet.length(); j++){
            Occ[i][j] = new uint [seq_len_ratio];
            ofile2.read((char*) Occ[i][j], seq_len_ratio*sizeof(uint));
            //cout<<Occ[i][j][1]<<endl;
            //Occ.push_back(p_Occ);
        }
        //cout<<"The readed OCC is: "<<endl;
        //for(uint m=0; m<seq_len_ratio; m++) cout<<Occ[i][0][m]<<endl;
        //cout<<"End"<<endl;
        for(uint j=0; j<p_alphabet.length(); j++){
            rev_Occ[i][j] = new uint [seq_len_ratio];
            ofile2.read((char*) rev_Occ[i][j], seq_len_ratio*sizeof(uint));
            //cout<<rev_Occ[i][j][1]<<endl;
            //rev_Occ.push_back(p_rev_Occ);
        }
        //for(uint m=0; m<seq_len_ratio; m++) cout<<rev_Occ[i][0][m]<<endl;
        //cout<<"End"<<endl;
        cout<<"Occ readed!"<<endl;

        uint* index_array = new uint [p_read_num];
        uint* rev_index_array = new uint [p_read_num];
        ofile3.read((char*) index_array, p_read_num*sizeof(uint));
        seq_index_array.push_back(index_array);
        ofile3.read((char*) rev_index_array, p_read_num*sizeof(uint));
        rev_seq_index_array.push_back(rev_index_array);
        //cout<<index_array[1]<<'\t'<<rev_index_array[1]<<endl;
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
    char seed_out[] = "seed_reads.fa";
    output_fasta(seed_out, seeds, reads_title);

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
        cout<<"Seeds number: "<<seeds.size()<<endl;
        //for(auto it=seeds.begin(); it!=seeds.end(); it++) cout<<it->first<<'\t'<<it->second<<endl;
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
    delete [] rev_Occ;


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
