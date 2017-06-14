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
// echo "DLBCL021-Tumor.mutatedrows.txt selector.bed Sample_DLBCL021_Normal.singleindex-deduped.sorted.freq.paired.Q30.txt Sample_DLBCL021_Tumor.singleindex-deduped.sorted.freq.paired.Q30.txt DLBCL021-Tumor.depth.txt 50 windows.txt" | ./a.out

const int N_BASES=3*1e3; //max # bases per tile  
const int N_TILES=3*1e3; //max # tiles
const double MUTATION_DEPTH_THRESHOLD=0.01;


struct Base{
  stack<string> sRead, eRead; //start read & end read
  int sDepth=0, eDepth=0;
};
struct Tile{
  string chr;
  int start, end;
};  


ofstream fOut;
ifstream fSelector, fRows, fNormal, fTumor, fDepths;

map<string, string> normalsAndErrors; //mutations that should be ignored during processTiles()
Base bases[N_TILES][N_BASES];
Tile tiles[N_TILES];
int windowLength, readLength=-1;


bool aboveThreshold(int depth, int totalPairs){
  if((double)totalPairs/(double)depth > MUTATION_DEPTH_THRESHOLD)
    return true;
  return false;
}


int normalPairs[8];
void processNormalMutations(){
  string filler;
  fNormal >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler;

  while(!fNormal.eof()){
    string chr;
    int pos, depth;
    fNormal >> chr >> pos >> depth >> filler >> filler >> filler >> normalPairs[0] >> normalPairs[1] >> normalPairs[2] >> normalPairs[3] >> normalPairs[4] >> normalPairs[5] >> normalPairs[6] >> normalPairs[7];
    if(chr == "") continue; 

    string mutations="";
    if(aboveThreshold(depth, normalPairs[0]+normalPairs[1])) mutations += 'A';
    if(aboveThreshold(depth, normalPairs[2]+normalPairs[3])) mutations += 'C';
    if(aboveThreshold(depth, normalPairs[4]+normalPairs[5])) mutations += 'T';
    if(aboveThreshold(depth, normalPairs[6]+normalPairs[7])) mutations += 'G';

    string key = ""; key += chr; key += ":"; key += to_string(pos);
    if(mutations != "")
      normalsAndErrors[key] = mutations;
  }
}


void calculateDepths(int tileNum){
  while(!fDepths.eof()){
    int start, length; string chr, L;
    fDepths >> chr >> start >> L;
    if(L=="") continue;
    L.pop_back(); length = stoi(L);

    if(chr != tiles[tileNum].chr || start>tiles[tileNum].end) //this means first read of next tile was just processed
      break;
    int tilePosition = start+length-tiles[tileNum].start-windowLength+1;;
    if(tilePosition<0 || start+length-1>tiles[tileNum].end) //read not completely in tile
      continue;

    bases[tileNum][tilePosition].sDepth++;
    bases[tileNum][tilePosition+length].eDepth++;
  }
}


void processSequence(int tileNum){
  while(!fRows.eof()){
    int start, length; string chr, L, read;
    fRows >> chr >> start >> L >> read;
    if(read == "") continue; //blank line
    L.pop_back(); length = stoi(L);
    readLength = length;

    if(chr != tiles[tileNum].chr || start>tiles[tileNum].end) //this means first read of next tile was just processed
      break;
    int tilePosition = start+length-tiles[tileNum].start-windowLength+1;
    if(tilePosition<0 || start+length-1>tiles[tileNum].end) //read not completely in tile
      continue;
    if(windowLength>length){
      cout << "-----Warning: Window length is greater than the length of a read --window defaulted to length of read (" << L << ")" << endl;
      windowLength = length;
    }

    string mutations="", key;
    map<string, string>::iterator it;
    for(int i=0; i<read.size(); i++){
      if(read[i]=='=') 
        continue;  

      key = ""; key += chr; key += ":"; key += to_string(start+i); 
      it = normalsAndErrors.find(key);
      if(it != normalsAndErrors.end())
        if((it->second).find(read[i]) != string::npos)
          continue;

      mutations += chr;
      mutations += "-";
      mutations += to_string(start+i);
      mutations += "-";
      mutations += read[i];
      mutations += "\t";
    }

    if(mutations != ""){
      bases[tileNum][tilePosition].sRead.push(mutations);
      bases[tileNum][tilePosition+length].eRead.push(mutations);
    }
  }
}


int tileNum=0;
void processTiles(){
  while(!fSelector.eof()){
    int sTile, eTile; string chr;
    fSelector >> chr >> sTile >> eTile;
    if(chr == "") continue; //blank line
    tiles[tileNum] = {chr, sTile, eTile};

    assert(tileNum<N_TILES-5 && "Number of tiles exceeds N_TILES, change the constant in the code to allocate more memory.");
    assert(eTile-sTile<N_BASES && "Number of base pairs in a single tile exceeds N_BASES, change the constant in the code to allocate more memory.");
    processSequence(tileNum);
    calculateDepths(tileNum);
    tileNum++;
  }
}

int errorPairs[8];
void processErrors(){
  string filler;
  fTumor >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler >> filler;

  while(!fTumor.eof()){
    string chr;
    int pos, depth;
    fTumor >> chr >> pos >> depth >> filler >> filler >> filler >> errorPairs[0] >> errorPairs[1] >> errorPairs[2] >> errorPairs[3] >> errorPairs[4] >> errorPairs[5] >> errorPairs[6] >> errorPairs[7];
    if(chr == "") continue; 

    string mutations="";
    if(errorPairs[0]+errorPairs[1]==1) mutations += 'A';
    if(errorPairs[2]+errorPairs[3]==1) mutations += 'C';
    if(errorPairs[4]+errorPairs[5]==1) mutations += 'T';
    if(errorPairs[6]+errorPairs[7]==1) mutations += 'G';

    string key = ""; key += chr; key += ":"; key += to_string(pos);
    if(mutations != "")
      normalsAndErrors[key] = mutations;
  }
}


int convertMutationToPosition(string mutation){
  string pos = "";
  bool dashFound = false;
  for(char c : mutation){
    if(c == '-' && dashFound)
      break;
    if(dashFound)
      pos += c;
    if(c == '-')
      dashFound = true;
  }
  return stoi(pos);
}


void generateWindows(){
  map<string, int> mutations; //<mutations, count>
  int depthCtr=0;
  for(int i=0; i<tileNum; i++){
    Tile curTile = tiles[i];
    for(int j=0; j<=curTile.end-curTile.start-windowLength+1+readLength*2; j++){
      map<string, int> dedupedOutput;
      //add read-start events to current map of mutations
      pair<map<string, int>::iterator, bool> ret;
      while(!bases[i][j].sRead.empty()){
        ret = mutations.emplace(bases[i][j].sRead.top(), 1);
        bases[i][j].sRead.pop();
        if(ret.second == false) //same set of mutations already occured, increases count
          (ret.first->second)++;
      }
      depthCtr += bases[i][j].sDepth;

      int windowStart = j-readLength+curTile.start+windowLength-1;
      int windowEnd = windowStart+windowLength-1;

      if(windowStart >= curTile.start && windowEnd <= curTile.end){
        for(map<string, int>::iterator it=mutations.begin(); it != mutations.end(); ++it){
          string mutationsInRange = "", currentMutation;
          istringstream totalMutations(it->first);
          while(totalMutations >> currentMutation){
            int mutationPosition = convertMutationToPosition(currentMutation);
            if(mutationPosition >= windowStart && mutationPosition <= windowEnd){
              mutationsInRange += currentMutation;
              mutationsInRange += "\t";
            }
          }
          if(mutationsInRange != ""){
            ret = dedupedOutput.emplace(mutationsInRange, 1);
            if(ret.second == false)
              (ret.first->second)++;
          }
        }

        for(map<string, int>::iterator it=dedupedOutput.begin(); it != dedupedOutput.end(); ++it)
          fOut << "\t" << it->second << "\t" << it->first << curTile.chr << ":" << windowStart << "-" << windowEnd << "\t" << depthCtr << endl;
      }

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
      depthCtr -= bases[i][j+windowLength].eDepth;
      
    }
  } 
}


int main(){ // <mutatedrows> <selector> <normal> <tumor> <depth> <windowlength> <output>
  string rows, selector, normal, tumor, depth, output;
  cin >> rows >> selector >> normal >> tumor >> depth >> windowLength >> output;

  fOut.open(output);  
  fSelector.open(selector);
  fRows.open(rows);
  fNormal.open(normal);
  fTumor.open(tumor);
  fDepths.open(depth);

  clock_t t;
  t = clock();
  cout << "------Program start (generateWindows)..." << endl;

  cout << "----Processing normal mutations..." << endl;
  processNormalMutations();
  cout << "----Processing errors..." << endl;
  processErrors();
  cout << "----Processing and filtering reads..." << endl;
  processTiles();
  cout << "----Generating windows..." << endl;
  generateWindows();

  cout << "----Done" << endl;
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << endl;

  fOut.close();
  fSelector.close();
  fRows.close();
  fNormal.close();
  fTumor.close();
  fDepths.close();
  return 0;
}