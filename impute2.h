//impute2.h
#ifndef IMPUTE2_H
#define IMPUTE2_H
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

void process_impute2_probs(string childFile,string motherFile,string cSamp,string mSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample,bool interaction);
void process_child_impute2_probs(filtering_ostream& outpat, filtering_ostream& outmat, filtering_ostream& outhet,string toprint,int cSnp,vector<float>& cArray,vector<float>& mArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr);
void process_impute2_probs_bgen(string childFile,string motherFile,string cSamp,string mSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample,bool interaction);
void process_child_impute2_probs_interaction(filtering_ostream& outhet,string toprint,int cSnp,vector<float>& cArray,vector<float>& mArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr);

#endif
