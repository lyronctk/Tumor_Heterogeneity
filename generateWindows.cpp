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

const int N_BASES=3*1e3; //max # bases per tile  
const int N_TILES=3*1e3; //max # tiles


struct Base{
  stack<string> sRead, eRead; //start read & end read
};
struct Tile{
  string chr;
  int start, end;
};  


ofstream fOut;
ifstream fSelector, fRows;

Base bases[N_TILES][N_BASES];
Tile tiles[N_TILES];
int windowLength;


void processSequence(int tileNum){
  while(!fRows.eof()){
    int start, length; string chr, L, read;
    fRows >> chr >> start >> L >> read;
    if(read == "") continue; //blank line
    L.pop_back();
    length = stoi(L);

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
      mutations += chr;
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
  string lastChr="";
  while(!fSelector.eof()){
    int sTile, eTile; string chr;
    fSelector >> chr >> sTile >> eTile;
    if(chr == "") continue; //blank line
    tiles[tileNum] = {chr, sTile, eTile};

    if(lastChr != chr){
      lastChr = chr;
      cout << "----Processing " << chr << "..." << endl;
    }

    assert(tileNum<N_TILES-5 && "Number of tiles exceeds N_TILES, change the constant in the code to allocate more memory.");
    assert(eTile-sTile<N_BASES && "Number of base pairs in a single tile exceeds N_BASES, change the constant in the code to allocate more memory.");
    processSequence(tileNum);
    tileNum++;
  }
}


void createWindows(){
  map<string, int> mutations; //<mutations, count>
  string lastChr = "";
  for(int i=0; i<tileNum; i++){
    Tile curTile = tiles[i];
    if(curTile.chr != lastChr){
      lastChr = curTile.chr;
      cout << "----Generating windows in " << lastChr << "..." << endl;
    }
    for(int j=0; j<=curTile.end-curTile.start-windowLength+1; j++){
      //add read-start events to current map of mutations
      pair<map<string, int>::iterator, bool> ret;
      while(!bases[i][j].sRead.empty()){
        ret = mutations.emplace(bases[i][j].sRead.top(), 1);
        bases[i][j].sRead.pop();
        if(ret.second == false) //same set of mutations already occured, increases count
          (ret.first->second)++;
      }

      //fOut << chr << ":" << j << "-" << j+windowLength-1 << " " << depth << endl;
      for(map<string, int>::iterator it=mutations.begin(); it != mutations.end(); ++it)
        fOut << it->second << " " << it->first << curTile.chr << ":" << j+curTile.start << "-" << j+curTile.start+windowLength-1 << endl;

      //delete read-end events from current map of mutations
      map<string, int>::iterator it;
      while(!bases[i][j+windowLength].eRead.empty()){
        it = mutations.find(bases[i][j+windowLength].eRead.top());
        bases[i][j+windowLength].eRead.pop();
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
  cout << " *" << endl;
  createWindows();

  cout << "----Done" << endl;
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << endl;

  fOut.close();
  fSelector.close();
  fRows.close();
  return 0;
}