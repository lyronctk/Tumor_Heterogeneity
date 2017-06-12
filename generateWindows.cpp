#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <assert.h>
#include <time.h>
#include <stack>
#include <unordered_map>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ generateWindows.cpp
// ./a.out

const int MAX_POSITION=3*1e8 ;
const int NUM_CHROMOSOMES=23;
const string convert_chr[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"};

struct Base{
  stack<string> sRead, eRead; //start read & end read
};
struct Chr{
  Chr(): bases(MAX_POSITION){}
  vector<Base> bases; //bases are 1-indexed
  int minPosition, maxPosition;
};

Chr sequence[NUM_CHROMOSOMES]; //stores 'events' when a read starts/ends 0-indexed
int windowLength;
ofstream f;
void initializeSequence(){
  int last_chr=-1;
  while(!cin.eof()){
    int start, chr, length;
    string C, L, read;
    cin >> C >> start >> L >> read;
    L.pop_back();
    length = stoi(L);
    chr = C.back()=='X' ? 22 : ((int)(C.back()-'0'))-1;

    if(length<windowLength){
      cout << "-----Warning: Window length is greater than the length of a read --window defaulted to length of read (" << L << ")" << endl;
      windowLength = length;
    }
    assert(start+length<MAX_POSITION && "A sequence ends outside of MAX_POSITION range, adjust the constant in the code.");
    if(chr!=last_chr){
      sequence[chr].minPosition = start;
      last_chr = chr;
      cout << "----Initializing chr" << convert_chr[chr] << endl;
    }

    string mutations="";
    for(int i=0; i<read.size(); i++){
      if(read[i]=='=') 
        continue;
      mutations += "chr";
      mutations += to_string(chr);
      mutations += "-";
      mutations += to_string(start+i);
      mutations += "-";
      mutations += read[i];
      mutations += " ";
    }

    sequence[chr].bases[start].sRead.push(mutations);
    sequence[chr].bases[start+length].eRead.push(mutations);
    sequence[chr].maxPosition = start+length;
  }
}


void generateDepth(){
  for(int i=0; i<NUM_CHROMOSOMES; i++){
    int depth=0;
    for(int j=sequence[i].minPosition; j+windowLength<=sequence[i].maxPosition; j++){ //slide window
      depth += sequence[i].bases[j].sRead.size();
      if(depth > 0)
        f << "chr" << convert_chr[i] << ":" << j << "-" << j+windowLength-1 << " " << depth << endl;
      depth -= sequence[i].bases[j+windowLength].eRead.size();
    }
  }
}


int main(){ // <mutatedrows> <windowlength>
  string fileName;
  cin >> fileName >> windowLength;

  freopen("testcase.txt", "r", stdin); //fileName.c_str()
  f.open("windows.txt");  

  clock_t t;
  t = clock();
  cout << "------Generating..." << endl;

  initializeSequence();
  generateDepth();

  cout << "----Done" << endl;
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << endl;

  f.close();
  return 0;
}