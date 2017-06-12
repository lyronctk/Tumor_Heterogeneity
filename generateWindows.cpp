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

const int N_BASES=5*1e3; //max # bases per tile  //CONSTANT STILL NEEDS TO BE FIXED
const int N_TILES=1*1e3; //max # tiles
const string convert_chr[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"};


struct Base{
  stack<string> sRead, eRead; //start read & end read
};
struct Tile{
  int chr, start, end;
};  


ofstream fOut;
ifstream fSelector, fRows;

Base bases[N_TILES][N_BASES];
Tile tiles[N_TILES];
int windowLength;


void processSequence(int tileNum){
  while(!fRows.eof()){
    int start, chr, length; string C, L, read;
    fRows >> C >> start >> L >> read;
    L.pop_back();
    length = stoi(L);
    chr = C.back()=='X' ? 22 : ((int)(C.back()-'0'))-1;

    if(chr != tiles[tileNum].chr || start>tiles[tileNum].end) //this means first read of next tile was just processed
      break;
    int tilePosition = start-tiles[tileNum].start;
    if(tilePosition<0 || start+length-1>tiles[tileNum].end) //read not completely in tile
      continue;
    if(windowLength>length){
      cout << "-----Warning: Window length is greater than the length of a read --window defaulted to length of read (" << L << ")" << endl;
      windowLength = length;
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

    bases[tileNum][tilePosition].sRead.push(mutations);
    bases[tileNum][tilePosition+length].eRead.push(mutations);
  }
}


int tileNum=0;
void processTiles(){
  int lastChr=-1;
  while(!fSelector.eof()){
    int chr, sTile, eTile; string chrStr;
    cin >> chrStr >> sTile >> eTile;
    chr = chrStr.back()=='X' ? 22 : ((int)(chrStr.back()-'0'))-1;
    tiles[tileNum] = {chr, sTile, eTile};

    if(lastChr != chr){
      lastChr = chr;
      cout << "----Processing chr" << chr << "..." << endl;
    }

    assert(eTile-sTile<N_BASES && "Number of base pairs in a single tile exceeds N_BASES, change the constant in the code to allocate more memory.");
    processSequence(tileNum);
    tileNum++;
  }
}


void createWindows(){
  map<string, int> mutations; //<mutations, count>
  for(int i=0; i<tileNum; i++){
    Tile curTile = tiles[tileNum];
    for(int j=0; j<=curTile.end-curTile.start-windowLength+1; j++){
      //add read-start events to current map of mutations
      pair<map<string, int>::iterator, bool> ret;
      while(!bases[i][j].sRead.empty()){
        ret = mutations.emplace(bases[i][j].sRead.top(), 1);
        bases[i][j].sRead.pop();
        if(ret.second == false) //same set of mutations already occured, increases count
          (ret.first->second)++;
      }

      //f << "chr" << convert_chr[i] << ":" << j << "-" << j+windowLength-1 << " " << depth << endl;
      for(map<string, int>::iterator it=mutations.begin(); it != mutations.end(); ++it)
        fOut << it->second << " " << it->first << "chr" << convert_chr[curTile.chr] << ":" << j+curTile.start << "-" << j+curTile.start+windowLength-1 << endl;

      //delete read-end events from current map of mutations
      map<string, int>::iterator it;
      while(!bases[i][j+windowLength].eRead.empty()){
        it = mutations.find(bases[i][j].eRead.top());
        bases[i][j].eRead.pop();
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
  cin >> rows >> selector >> windowLength >> output;

  fOut.open(output);  
  fSelector.open(selector);
  fRows.open(rows);

  clock_t t;
  t = clock();
  cout << "------Generating..." << endl;

  processTiles();
  createWindows();

  cout << "----Done" << endl;
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << endl;

  fOut.close();
  fSelector.close();
  fRows.close();
  return 0;
}