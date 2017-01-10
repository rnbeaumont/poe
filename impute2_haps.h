//impute2.h
#ifndef IMPUTE2_HAPS_H
#define IMPUTE2_HAPS_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
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

void process_impute2_haps(string childFile,string motherFile,string cSamp,string mSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample,bool interactions);
void process_child_impute2_haps(filtering_ostream& outpat, filtering_ostream& outmat, filtering_ostream& outhet,vector<string>& toprint,int loopnum,int cSnp,vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr);
void process_child_impute2_haps_interaction(filtering_ostream& outhet,vector<string>& toprint,int loopnum,int cSnp,vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr);

#endif
