//impute2.h
#ifndef IMPUTE2_TRIO_H
#define IMPUTE2_TRIO_H
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

void process_impute2_trio_probs(string childFile,string motherFile,string fatherFile,string cSamp,string mSamp,string fSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample);
void process_child_impute2_trio_probs(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,string toprint,int cSnp,vector<float>& cArray,vector<float>& mArray,vector<float>& fArray,std::unordered_map<std::string,int>& mhash,std::unordered_map<std::string,int>& fhash,vector<string>& cids,string chr,bool transmitted);
void compare_trio_genotypes_unphased(vector<float>& cArray,vector<float>& mArray,vector<float>& fArray,int i,int mind,int find,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted);
void compare_pair_genotypes_unphased(vector<float>& cArray,vector<float>& mArray,int i,int mind,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted,bool mother);
void process_impute2_trio_probs_bgen(string childFile,string motherFile,string fatherFile,string cSamp,string mSamp,string fSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample);

#endif
