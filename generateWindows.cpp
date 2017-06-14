#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <time.h>
#include <vector>
#include <string>
#include <stack>
#include <map>
using namespace std;
// git add generateWindows.cpp filterMutatedRows.sh clonal_tumor_wrapper.sh
// clang++ -std=c++11 -stdlib=libc++ generateWindows.cpp
// echo "DLBCL021-Tumor.mutatedrows.txt selector.bed Sample_DLBCL021_Normal.singleindex-deduped.sorted.freq.paired.Q30.txt Sample_DLBCL021_Tumor.singleindex-deduped.sorted.freq.paired.Q30.txt DLBCL021-Tumor.depth.txt DLBCL021_Tumor_summarized.cosmic.split.txt 50 2 windows.txt" | ./a.out

const int N_BASES=3*1e3; //max # bases per tile  
const int N_TILES=3*1e3; //max # tiles
const int N_COSMIC=3;
const double MUTATION_DEPTH_THRESHOLD=0.01;


struct Base{
  stack<string> sRead, eRead; //start and end of deduped tumor reads 
  int sDepth=0, eDepth=0; //start and end of non-deduped tumor reads
};
struct Tile{
  string chr;
  int start, end;
};  


ofstream fOut;
ifstream fSelector, fRows, fNormal, fTumor, fDepths, fCosmic;

map<string, string> normalsAndErrors; //mutations that should be ignored during processTiles()
Base bases[N_TILES][N_BASES];
Tile tiles[N_TILES];
int windowLength, readLength=-1;


pair<int, int> allelePercentages[N_COSMIC];
double AP=0;
void processCosmicSplit(){
  int numLines=0;
  while(!fCosmic.eof()){
    int depth; string percentStr;
    fCosmic >> depth >> percentStr;
    if(percentStr == "") continue;
    int percent = stoi(percentStr.substr(0, 4));

    if(depth < 200)
      continue;
    if(numLines < N_COSMIC){
      allelePercentages[numLines].first = depth; 
      allelePercentages[numLines].second = percent;
    }
    numLines++;
  }
  if(numLines < N_COSMIC)
    cout << "--WARNING: There are less than " << N_COSMIC << " usable depth/AF pairs in the cosmic split file. The program will use the average AF of the " << numLines << " pairs."; 

  int nTop = min(N_COSMIC, (int)(0.05*numLines+0.5));
  if(nTop == 0)
    return;
  if(allelePercentages[0].second - allelePercentages[nTop-1].second > 10)
    cout << "--WARNING: The difference between the highest allele percentage and the lowest that are being averaged is greater than 10%";

  for(int i=0; i<nTop; i++)
    AP += allelePercentages[i].second;
  AP /= nTop;
  AP /= 100;
}


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
    streampos temp = fDepths.tellg();
    int start, length; string chr, L;
    fDepths >> chr >> start >> L;
    if(L=="") continue;
    L.pop_back(); length = stoi(L);

    if(chr != tiles[tileNum].chr || start>tiles[tileNum].end){ //this means first read of next tile was just processed
      fDepths.seekg(temp);
      break;
    }
    int tilePosition = start+length-tiles[tileNum].start-windowLength+1;;
    if(tilePosition<0) //read will not completely cover any window inside the tile
      continue;

    bases[tileNum][tilePosition].sDepth++;
    bases[tileNum][tilePosition+length].eDepth++;
  }
}


void processSequence(int tileNum){
  while(!fRows.eof()){
    streampos temp = fRows.tellg();
    int start, length; string chr, L, read;
    fRows >> chr >> start >> L >> read;
    if(read == "") continue; //blank line
    L.pop_back(); length = stoi(L);
    readLength = length;

    if(chr != tiles[tileNum].chr || start>tiles[tileNum].end){ //this means first read of next tile was just processed
      fRows.seekg(temp);
      break;
    }
    int tilePosition = start+length-tiles[tileNum].start-windowLength+1;
    if(tilePosition<0) //read will not completely cover any window inside the tile
      continue;
    if(windowLength>length){
      cout << "--WARNING: Window length is greater than the length of a read --window defaulted to length of read (" << L << ")" << '\n';
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

    assert(tileNum<N_TILES-5 && "--Number of tiles exceeds N_TILES, change the constant in the code to allocate more memory.");
    assert(eTile-sTile<N_BASES && "--Number of base pairs in a single tile exceeds N_BASES, change the constant in the code to allocate more memory.");
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


void generateWindows(int slide){
  map<string, int> mutations; //<mutations, count>
  int depthCtr=0, slideCtr;
  for(int i=0; i<tileNum; i++){
    Tile curTile = tiles[i];
    slideCtr = -1;
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
        slideCtr++;
        if(slideCtr%slide == 0){ //how much to shift sliding window every print 
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
              ret = dedupedOutput.emplace(mutationsInRange, it->second);
              if(ret.second == false)
                ret.first->second += it->second;
            }
          }

          //added to different map because an 'event' that is cut off by the window could match another shorter 'event', thus creating two separate but identical 'events'
          for(map<string, int>::iterator it=dedupedOutput.begin(); it != dedupedOutput.end(); ++it)
            fOut << "\t" << it->second << "\t" << it->first << curTile.chr << ":" << windowStart << "-" << windowEnd << "\t" << depthCtr << '\n';
        }
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


int main(){ // <mutatedrows> <selector> <normal> <tumor> <depth> <cosmic> <windowLength> <windowSlide> <output>
  std::ios::sync_with_stdio(false); cin.tie(NULL);

  int slide;
  string rows, selector, normal, tumor, depth, output, cosmic;
  cin >> rows >> selector >> normal >> tumor >> depth >> cosmic >> windowLength >> slide >> output;

  fOut.open(output);  
  fSelector.open(selector);
  fRows.open(rows);
  fNormal.open(normal);
  fTumor.open(tumor);
  fDepths.open(depth);
  fCosmic.open(cosmic);

  clock_t t;
  t = clock();
  cout << "------Program start (generateWindows)..." << '\n';

  cout << "----Processing cosmic split file..." << "\n";
  processCosmicSplit();
  cout << "----Processing normal mutations..." << '\n';
  processNormalMutations();
  cout << "----Processing errors..." << '\n';
  processErrors();
  cout << "----Processing and filtering reads..." << '\n';
  processTiles();
  cout << "----Generating windows..." << '\n';
  generateWindows(slide);

  cout << "----Done" << '\n';
  cout << "------Executed in " << ((float)(clock()-t))/CLOCKS_PER_SEC << " seconds." << '\n';

  fOut.close();
  fSelector.close();
  fRows.close();
  fNormal.close();
  fTumor.close();
  fDepths.close();
  fCosmic.close();
  return 0;
}