// general_functions.h
#ifndef GENERAL_FUNCTIONS_H
#define GENERAL_FUNCTIONS_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <unordered_map>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
//#include "libs/iostreams/src/zlib.cpp"
using namespace std;
using namespace boost::iostreams;
using namespace boost;

struct gen{
public:
	gen();			//default constructor
	gen(string file);	//constructor to open file when it's all written
	unsigned int get_num_lines() const {return nSnps;}	//get number of lines/snps
	unsigned int get_num_samples() const {return nSample;}	//get numbr of samples
	void set_num_samples(unsigned int nSamp);	//set number of samples
	vector<float> probs;	//probabilities for last read variant
	vector<vector<float>> haps;	//haplotypes array
	bool getline();		//get the next line
	vector<string> fields;	//line split
	void split_fields();	//convert fields to floats
	void split_haps(unsigned int z);	//convert fields to floats
	string outString();	//construct the start of the gen file line
	void read_num_lines();	//set nSnp by reading number of lines

private:
	string filename;	//filename
	unsigned int nSnps;	//number of lines/snps
	unsigned int nSample;	//number of samples
	boost::iostreams::filtering_istream openFile;	//boost file for reading
	std::ifstream iFile;	//file for reading
	unsigned int nSampLine1;	//number of Samples in the first SNP to check sample file corresponds to this file
	string line;		//last line read
};

struct bgen{
public:
	bgen();			//default constructor
	bgen(string file);	//constructor to open file when it's all written
	void read_variant_id();	//read the next variant ID
	void read_variant_probabilities();	//read the probabilities from the previously read variant
	void skip_variant_probabilities();	//skip the probabilities and go to the next variant data block
	unsigned int get_m() const {return M;}		// get M
	unsigned int get_n() const {return N;}		// get N
	string get_rsid() const {return rsidS;}		//get rsidS
	unsigned int get_pos() const {return pos;}	//get pos
	unsigned int cSnp() const {return currentSnp;}	//get currentSnp
	string get_aA() const {return aA;}		//get aA
	string get_aB() const {return aB;}		//get aB
	vector<float> get_probs() {return probs;}	// get probs
	vector<float> probs;	// probabilities for last read variant

private:
	unsigned int offset;	//offset from start of file where genotype blocks start
	unsigned int LH;	//length in bytes of the header
	unsigned int M;		//number of SNPs
	unsigned int N;		//number of individuals
	bool compressed;	//whether the genotype data is zlib compressed
	int layout;		//layout (bgen version v1.x)
	bool sampleIds;		//sampleIDs - whether sample IDs are stored in the file
	unsigned int nv;	//number of individuals in current data block (must equal N but for some reason stored in each block too)
	unsigned int K;		//number of alleles (if layout==2, else undefined)
	unsigned int currentSnp;		//counter to keep track of which SNP has just been read
	bool readVarHead;	//bool to say whether we're part way through a variant or not
	vector<unsigned char> vid,rsid;	//variant identifier and rsid of current variant
	vector<unsigned char> chr;	//chromosome identifier (not always chromosome, especially if encoded using qctool!!!!!!!!!!!!!)
	unsigned int pos;	//position of current variant
	unsigned int nbytes;	// length of compressed data block (if compressed, else uninitiated!!)
	vector<unsigned char> alleleB;	//alleles
	vector<unsigned char> alleleA;
	string rsidS;	//string version of rsid
	string aA,aB;	// string versions of the alleles
	ifstream input;	// ifstream containing opened bgen file
};

void read_pheno( string pFile, unordered_map<string, string>& phenoMap,string& pHead );
void read_sample( string infoFile, int numlines, std::unordered_map<std::string, int> hash, vector <string> cids );
void read_sample_2( string cSamp,string mSamp,std::unordered_map<std::string,int>& mhash,vector<string>& cids);
void read_sample_3( string cSamp,string mSamp,string fSamp,std::unordered_map<std::string,int>& mhash,std::unordered_map<std::string,int>& fhash,vector<string>& cids);
int get_num_lines( string lineFile );
void check_gen_length(string inFile,unsigned int nCol);
void check_dosage_length(string inFile,unsigned int nCol);
void read_info( string infoFile, int numlines, std::unordered_map<std::string, int> hash, vector <string> cids );
void read_info_2( string cInfo,string mInfo,std::unordered_map<std::string,int>& chash,vector<string>& mids,vector<string>& mstr,string chr);
void print_sample(string cSamp,string outFile);
void open_bgen(ifstream& input,bgen& opened_bgen);
bgen read_variant_id_bgen(bgen opened_bgen,ifstream& input);
bgen read_variant_probabilities_bgen(bgen opened_bgen,ifstream& input);
bgen skip_variant_probabilities_bgen(bgen opened_bgen,ifstream& input);

#endif
