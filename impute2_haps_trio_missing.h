//impute2_haps_trio_missing.h
#ifndef IMPUTE2_HAPS_TRIO_MISSING_H
#define IMPUTE2_HAPS_TRIO_MISSING_H
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

void process_impute2_haps_trio(string childFile,string motherFile,string fatherFile,string cSamp,string mSamp,string fSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample);
void process_child_impute2_haps_trio(filtering_ostream& outpat, filtering_ostream& outmat, filtering_ostream& outhet,vector<string>& toprint,int loopnum,int cSnp,vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,vector<vector<float>>& fHapArray,std::unordered_map<std::string,int>& mhash,std::unordered_map<std::string,int>& fhash,vector<string>& cids,string chr,bool transmitted);
void compare_trio_genotypes(vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,vector<vector<float>>& fHapArray,int j,int i,int mind,int find,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted,int loopnum);
void compare_pair_genotypes(vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,int j,int i,int mind,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted,int loopnum,bool mother);

#endif
