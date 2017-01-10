#ifndef MACH_DOSAGES_H
#define MACH_DOSAGES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unordered_map>
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/algorithm/string.hpp>
//#include "libs/iostreams/src/zlib.cpp"
using namespace std;
using namespace boost::iostreams;
using namespace boost;

void process_mach(string childFile,string motherFile,string cInfo,string mInfo,string outFile,string pheno,string chr,bool transmitted,unsigned int nsnp,bool interactions);
void process_child_mach(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,int loopStart,int loopEnd,int cInds,unordered_map<string,int>& chash,vector<string>& mids,vector<string>& mstr,unordered_map<string, string>& pMap,vector<vector<float>>& mHapArray,unordered_map<string,int>& mline,string pHead,string outFile,string childFile);
void process_child_mach_interaction(filtering_ostream& outhet,int loopStart,int loopEnd,int cInds,unordered_map<string,int>& chash,vector<string>& mids,vector<string>& mstr,unordered_map<string, string>& pMap,vector<vector<float>>& mHapArray,unordered_map<string,int>& mline,string pHead,string outFile,string childFile);

#endif
