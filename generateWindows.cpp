#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <assert.h>
#include <time.h>
#include <stack>
#include <map>
#include <vector>
using namespace std;
// git add generateWindows.cpp filterMutatedRows.sh clonal_tumor_wrapper.sh
// clang++ -std=c++11 -stdlib=libc++ generateWindows.cpp
// ./a.out


const int N_BASES=2*1e5; //max # bases per tile  //CONSTANT STILL NEEDS TO BE FIXED
const int N_TILES=1*1e3; //max # tiles
const int N_CHROMOSOMES=23;
const string convert_chr[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"};


struct Base{
  stack<string> sRead, eRead; //start read & end read
};
struct Chr{
  Chr(): bases(N_TILES, vector<Base>(N_BASES)){} //stores 'events' when a read starts/ends
  vector<vector<Base> > bases; //[tile#][base(1-indexed)]
};
struct Tile{
  int chr, sTile, eTile;
};  


ofstream fOut;
ifstream fSelector, fRows;

Chr sequence[N_CHROMOSOMES]; //0-indexed
int tiles[N_TILES], windowLength;


void initializeSequence(int tilePos){
  int start, chr, length;
    string C, L, read;
    fRows >> C >> start >> L >> read;
    L.pop_back();
    length = stoi(L);

    if(length<windowLength){
      cout << "-----Warning: Window length is greater than the length of a read --window defaulted to length of read (" << L << ")" << endl;
      windowLength = length;
    }
    assert(start+length<MAX_POSITION && "Total number of bases in reads exceeds N_BASES, change the constant in code to allocate more memory.");

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


void initializeTiles(){
  int tilePos=0, prefixSum=0;
  while(!fSelector.eof()){
    int selectorChr, sTile, eTile; string selectorChrStr;
    cin >> selectorChrStr >> sTile >> eTile;
    selectorChr = selectorChrStr.back()=='X' ? 22 : ((int)(selectorChrStr.back()-'0'))-1;
    tiles[tilePos] = {chr, sTile, eTile};

    initializeSequence(tilePos);
    tilePos++;
  }
}


void generateWindows(){
  map<string, int> mutations; //<mutation_set, count>
  for(int i=0; i<NUM_CHROMOSOMES; i++){
    for(int j=sequence[i].minPosition; j+windowLength<=sequence[i].maxPosition; j++){ //slide window
      //add read-start events to current map of mutations
      pair<map<string, int>::iterator, bool> ret;
      while(!sequence[i].bases[j].sRead.empty()){
        ret = mutations.emplace(sequence[i].bases[j].sRead.top(), 1);
        sequence[i].bases[j].sRead.pop();
        if(ret.second == false) //same set of mutations already occured, increases count
          (ret.first->second)++;
      }

  //  f << "chr" << convert_chr[i] << ":" << j << "-" << j+windowLength-1 << " " << depth << endl;
      for(map<string, int>::iterator it=mutations.begin(); it != mutations.end(); ++it)
        fOut << it->second << " " << it->first << "chr" << convert_chr[i] << ":" << j << "-" << j+windowLength-1 << endl;

      //delete read-end events from current map of mutations
      map<string, int>::iterator it;
      while(!sequence[i].bases[j+windowLength].eRead.empty()){
        it = mutations.find(sequence[i].bases[j].eRead.top());
        sequence[i].bases[j].eRead.pop();
        if(it->second > 1) //more than one occurence of the set of mutations, just decrease count instead of erasing
          (it->second)--;
        else
          mutations.erase(it);
      }
    }
  }
}


int main(){ // <mutatedrows> <selector> <windowlength> <output>
  string rows, selector, output;
  cin >> fileName >> selector >> windowLength >> output;

  fOut.open(output);  
  fSelector.open(selector);
  fRows.open(fileName);

  clock_t t;
  t = clock();
  cout << "------Generating..." << endl;

  initializeTiles();
  generateWindows();

  cout << "----Done" << endl;
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << endl;

  f.close();
  return 0;
}