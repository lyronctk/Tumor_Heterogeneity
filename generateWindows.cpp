#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <assert.h>
#include <time.h>
using namespace std;
// clang++ -std=c++11 -stdlib=libc++ generateWindows.cpp
// ./a.out

const int MAX_POSITION=30; //*1e8 
const int NUM_CHROMOSOMES=23;
const string convert_chr[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"};

struct Base{
  vector<int> sRead, eRead; //start read & end read
};
struct Chr{
  Chr(): bases(MAX_POSITION){}
  vector<Base> bases; //bases are 1-indexed
  int minPosition, maxPosition;
};
struct FirstRead{
  int start, chr, length;
  string read;
};


Chr sequence; //stores 'events' when a read starts/ends 
int windowLength;
ofstream f;
FirstRead firstRead;
void initializeSequence(int chrNum){
  int readId=2;
  bool firstLoop = true;
  while(!cin.eof()){
    int start, chr, length;
    string C, L, read;
    cin >> C >> start >> L >> read;
    L.pop_back();
    length = stoi(L);
    chr = C.back()=='X' ? 22 : ((int)(C.back()-'0'))-1;

    if(length<windowLength){
      cout << "Warning: Window length is greater than the length of a read --window defaulted to length of read (" << L << ")" << endl;
      windowLength = length;
    }
    assert(start+length<MAX_POSITION && "A sequence ends outside of MAX_POSITION range, adjust the constant in the code.");
    if(chr != chrNum) {
      firstRead = {start, chr, length, read};
      break;
    }
    if(chr == 0){
      sequence.minPosition = start;
    }

    sequence.bases[start].sRead.push_back(readId);
    sequence.bases[start+length].eRead.push_back(readId);
    sequence.maxPosition = start+length;

    readId++;
  }
}


void generateDepth(int chrNum){
  f << "ON CHR" << convert_chr[chrNum] << "| minPosition: " << sequence.minPosition << " | maxPosition: " << sequence.maxPosition << endl;
  int depth=0;
  for(int i=sequence.minPosition; i+windowLength<=sequence.maxPosition; i++){ //slide window
    depth += sequence.bases[i].sRead.size();
    sequence.bases[i].sRead.clear(); //reset for next sequence

    if(depth > 0)
      f << "chr" << convert_chr[chrNum] << ":" << i << "-" << i+windowLength-1 << " " << depth << endl;

    depth -= sequence.bases[i+windowLength].eRead.size();
    sequence.bases[i+windowLength].eRead.clear();
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
  for(int i=0; i<NUM_CHROMOSOMES; i++){
    cout << "----" << "Processing chr" << convert_chr[i] << endl;
    if(i > 0){
      sequence.bases[firstRead.start].sRead.push_back(1);
      sequence.bases[firstRead.start+firstRead.length].eRead.push_back(1);
      sequence.minPosition = firstRead.start;
    }
    initializeSequence(i);
    generateDepth(i);
  }
  cout << "----Done" << endl;
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << endl;

  f.close();
  return 0;
}