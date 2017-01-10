#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include "general_functions.h"
#include "impute2.h"
#include "impute2_haps.h"
#include "impute2_trio.h"
#include "impute2_haps_trio.h"
#include "mach.h"
//#include "libs/iostreams/src/zlib.cpp"
// #include "boost/iostreams/filtering_streambuf.hpp"
// #include "boost/iostreams/filtering_stream.hpp"
// #include "boost/iostreams/copy.hpp"
// #include "boost/iostreams/filter/gzip.hpp"
// #include "boost/algorithm/string.hpp"
#define STRLENPAR 500
using namespace std;
struct Parameters {
	string childFile;
	string motherFile;
	string fatherFile;
	string cSamp;
	string mSamp;
	string fSamp;
	string cInfo;
	string mInfo;
	string fInfo;
	string outFile;
	string pheno;
	unsigned int nsnp;
	string chr;
	bool unordered;
	bool phased;
	bool transmitted;
	bool sample;
	bool minimac;
	bool interaction;
	int childN;
	int childBgenN;
	int cInfN;
	int motherN;
	int mInfN;
	int fatherN;
	int fInfN;
	int outN;
	int chrN;
	int cSampN;
	int mSampN;
	int fSampN;
	int phenoN;
};
Parameters read_options(int &argc,char *argv[]);
void print_usage();
void print_options(Parameters result);
void process_dosages();
void print_header();
void print_required(Parameters result);

int main( int argc, char *argv[] ){
	using namespace std;
	print_header();
	Parameters inUse=read_options(argc,argv);
/* CHECK IF ALL OPTIONS SPECIFIED AND PRINT USAGE IF THERE AREN'T ENOUGH ARGUMENTS!!!! IF THERE ARE THEN PRINT THE OPTIONS */
	if(inUse.childN==0 && inUse.childBgenN==0 && inUse.cSampN==0 && inUse.motherN==0 && inUse.mSampN==0 && inUse.outN==0 && inUse.chrN==0){
		print_usage();
	}
	if(!((inUse.childN==1 || inUse.childBgenN==1) &&inUse.cSampN==1 && inUse.motherN==1 && inUse.mSampN==1 && inUse.outN==1 && inUse.chrN==1)){
		print_options(inUse);
		print_required(inUse);
		print_usage();
	}
	print_options(inUse);
//	Decide which process to run based on supplied options
	if(inUse.cSampN==1 && inUse.mSampN==1 && inUse.phased==0 && inUse.fatherN==0 && inUse.fSampN==0 && inUse.childBgenN==0 && inUse.minimac==false){	//Impute2 dyad unphased
		if(inUse.childFile.compare(inUse.childFile.length()-4,4,"bgen")==0){
			if(inUse.transmitted){
				process_impute2_probs_bgen(inUse.motherFile,inUse.childFile,inUse.mSamp,inUse.cSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample,inUse.interaction);
			}else{
				process_impute2_probs_bgen(inUse.childFile,inUse.motherFile,inUse.cSamp,inUse.mSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample,inUse.interaction);
			}
		}else{
			if(inUse.transmitted){
				process_impute2_probs(inUse.motherFile,inUse.childFile,inUse.mSamp,inUse.cSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample,inUse.interaction);
			}else{
				process_impute2_probs(inUse.childFile,inUse.motherFile,inUse.cSamp,inUse.mSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample,inUse.interaction);
			}
		}
	}else if(inUse.cSampN==1&& inUse.mSampN==1 && inUse.phased==1 && inUse.fatherN==0 && inUse.fSampN==0 && inUse.childBgenN==0 && inUse.minimac==false){ //Impute2 dyad phased
		if(inUse.transmitted){
			process_impute2_haps(inUse.motherFile,inUse.childFile,inUse.mSamp,inUse.cSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample,inUse.interaction);
		}else{
			process_impute2_haps(inUse.childFile,inUse.motherFile,inUse.cSamp,inUse.mSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample,inUse.interaction);
		}
	}else if(inUse.cSampN==1 && inUse.mSampN==1 && inUse.phenoN==1 && inUse.fatherN==0 && inUse.fSampN==0 && inUse.minimac==true){ //Mach dyad unphased
		if(inUse.transmitted){
			process_mach(inUse.motherFile,inUse.childFile,inUse.mSamp,inUse.cSamp,inUse.outFile,inUse.pheno,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.interaction);
		}else{
			process_mach(inUse.childFile,inUse.motherFile,inUse.cSamp,inUse.mSamp,inUse.outFile,inUse.pheno,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.interaction);
		}
	}else if(inUse.cSampN==1 && inUse.mSampN==1 && inUse.phased==0 && inUse.fatherN==1 && inUse.fSampN==1 && inUse.minimac==false && inUse.interaction==false){ //Impute2 trio unphased
		if(inUse.childFile.compare(inUse.childFile.length()-4,4,"bgen")==0){
			process_impute2_trio_probs_bgen(inUse.childFile,inUse.motherFile,inUse.fatherFile,inUse.cSamp,inUse.mSamp,inUse.fSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample);
		}else{
			process_impute2_trio_probs(inUse.childFile,inUse.motherFile,inUse.fatherFile,inUse.cSamp,inUse.mSamp,inUse.fSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample);
		}
	}else if(inUse.cSampN==1 && inUse.mSampN==1 && inUse.phased==1 && inUse.fatherN==1 && inUse.fSampN==1 && inUse.minimac==false && inUse.interaction==false){ //Impute2 trio phased
		process_impute2_haps_trio(inUse.childFile,inUse.motherFile,inUse.fatherFile,inUse.cSamp,inUse.mSamp,inUse.fSamp,inUse.outFile,inUse.chr,inUse.transmitted,inUse.nsnp,inUse.sample);
	}else{
		cout<<"Invaid/incompatible options supplied"<<endl;
		return 1;
	}
	puts( "DONE!\n\n" );
	return 0;
}

Parameters read_options(int &argc,char *argv[]){
	int m,n,	/* Loop counters */
	    l,		/* String length */
	    x;		/* Exit code. */
        char argcase[10];	/* List buffer */
	Parameters result;
//initialise result
	result.nsnp=100000;
	result.unordered=false;
	result.phased=false;
	result.transmitted=false;
	result.sample=true;
	result.minimac=false;
	result.interaction=false;
	result.childN=0;
	result.childBgenN=0;
	result.cInfN=0;
	result.motherN=0;
	result.mInfN=0;
	result.outN=0;
	result.chrN=0;
	result.cSampN=0;
	result.mSampN=0;
	result.phenoN=0;

	for( n=1; n < argc; n++ )	/* Scan through args */
	{
		switch( (int)argv[n][0] )	/* Check for option character */
		{
/*		case "--":	*/
		case '-': x = 0;	/* Bail out if 1. */
			  l = strlen( argv[n] );
			for( m = 1; m < l; ++m )	/* scan through options */
			{
/*				argcase = &argv[n][1];*/
				if( strcmp(&argv[n][1],"child")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -child\n" );
						exit( 1 );
					}else{
						result.childFile=&argv[n+1][0];
						result.childN=1;
/*						printf( "String = %s\n", s );*/
					}
					x = 1;
				}else if( strcmp(&argv[n][1],"child_bgen")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -child\n" );
						exit( 1 );
					}else{
						result.childFile=&argv[n+1][0];
						result.childBgenN=1;
/*						printf( "String = %s\n", s );*/
					}
					x = 1;
				}else if( strcmp(&argv[n][1],"child_sample")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -child_sample\n" );
						exit( 1 );
					}else{
						result.cSamp=&argv[n+1][0];
						result.cSampN=1;
					}
				}else if( strcmp(&argv[n][1],"child_info")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -child_info\n" );
						exit( 1 );
					}else{
						result.cInfo=&argv[n+1][0];
						result.cInfN=1;
					}
				}else if( strcmp(&argv[n][1],"mother")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -mother\n" );
						exit( 1 );
					}else{
						result.motherFile=&argv[n+1][0];
						result.motherN=1;
					}
					x = 1;
				}else if( strcmp(&argv[n][1],"mother_sample")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -mother_sample\n" );
						exit( 1 );
					}else{
						result.mSamp=&argv[n+1][0];
						result.mSampN=1;
					}
				}else if( strcmp(&argv[n][1],"mother_info")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -mother_info\n" );
						exit( 1 );
					}else{
						result.mInfo=&argv[n+1][0];
						result.mInfN=1;
					}
				}else if( strcmp(&argv[n][1],"father")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -father\n" );
						exit( 1 );
					}else{
						result.fatherFile=&argv[n+1][0];
						result.fatherN=1;
					}
					x = 1;
				}else if( strcmp(&argv[n][1],"father_sample")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -father_sample\n" );
						exit( 1 );
					}else{
						result.fSamp=&argv[n+1][0];
						result.fSampN=1;
					}
				}else if( strcmp(&argv[n][1],"father_info")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -father_info\n" );
						exit( 1 );
					}else{
						result.fInfo=&argv[n+1][0];
						result.fInfN=1;
					}
				}else if( strcmp(&argv[n][1],"pheno")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -pheno\n" );
						exit( 1 );
					}else{
						result.pheno=&argv[n+1][0];
						result.phenoN=1;
					}
				}else if( strcmp(&argv[n][1],"out")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string! following option -out\n" );
						exit( 1 );
					}else{
						result.outFile=&argv[n+1][0];
						result.outN=1;
					}
				}else if( strcmp(&argv[n][1],"chr")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -chr\n" );
						exit( 1 );
					}else{
						result.chr=&argv[n+1][0];
						//sscanf(&argv[n+1][0], "%d", &chr);	/* should really check that sscanf returns either 1 or 2 and that chr is between 1 and 22 */
						result.chrN=1;
					}
				}else if( strcmp(&argv[n][1],"n_snp")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -n_snp\n" );
						exit( 1 );
					}else{
						sscanf(&argv[n+1][0], "%d", &result.nsnp);	/* do some error checking!!! */
					}
				}else if( strcmp(&argv[n][1],"unordered")==0 ){
					result.unordered=true;	/* do some error checking!!! */
				}else if( strcmp(&argv[n][1],"phased")==0 ){
					result.phased=true; /* do some error checking!!! */
				}else if( strcmp(&argv[n][1],"transmitted")==0 ){
					result.transmitted=true;	/* do some error checking!!! */
				}else if( strcmp(&argv[n][1],"no_print_sample")==0 ){
					result.sample=false;	/* do some error checking!!! */
				}else if( strcmp(&argv[n][1],"minimac")==0 ){
					result.minimac=true;
				}else if( strcmp(&argv[n][1],"interaction")==0 ){
					result.interaction=true;
				}else{
					strcpy( argcase, &argv[n][1] );
					printf( "Unrecognised option : %s\n", &argv[n][0] ); /* THIS NEEDS FIXING TO OUTPUT THE WHOLE ARGUMENT!!!!!!! */
					print_usage() ;
				}
				if( x == 1 ) {
					break;
				}
			}
			break;
		default:  /*printf( "Unrecognised option : %s\n", &argv[n][0] );	Not option -- text. */
			break;
		}
	}
	return result;
}

void print_usage( ){
	printf( "\nUsage:\n" );
	printf( "\t-child           child's dosage file\n" );
	printf( "\t-child_sample    childs sample file\n" );
	printf( "\t-mother          mothers dosage file\n" );
	printf( "\t-mother_sample   mothers sample file\n" );
	printf( "\t-father          fathers dosage file\n" );
	printf( "\t-father_sample   fathers sample file\n" );
	printf( "\t-out             output_file_prefix\n" );
	printf( "\t-chr             chromosome\n" );
//	printf( "\t-pheno           phenotype file\n" );
	printf( "\t-n_snp           (optional) number of snps to evaluate at once Default: 100000.\n" );
	printf( "\t                 Increasing this number dramatically increases memory usage\n" );
	printf( "\t-phased          (optional) phased haplotypes in impute2 format provided\n" );
	printf( "\t-transmitted     (optional) generate mother's transmitted and untransmitted alleles\n" );
	printf( "\t                 instead of maternally and paternally inherited alleles of fetus\n" );
	printf( "\n\tFor a full list of options see accompanying documentation\n\n\n");
	exit(0);
	return;
}

void print_options(Parameters result){
        cout<<"\nOptions in use:"<<endl;
        if(result.childN){
                cout<<"-child            "<<result.childFile<<endl;
        }
        if(result.cSampN){
                cout<<"-child_sample     "<<result.cSamp<<endl;
        }
        if(result.cInfN){
                cout<<"-child_info       "<<result.cInfo<<endl;
        }
        if(result.motherN){
                cout<<"-mother           "<<result.motherFile<<endl;
        }
        if(result.mSampN){
                cout<<"-mother_sample    "<<result.mSamp<<endl;
        }
        if(result.mInfN){
                cout<<"-mother_info      "<<result.mInfo<<endl;
        }
        if(result.fatherN){
                cout<<"-father           "<<result.fatherFile<<endl;
        }
        if(result.fSampN){
                cout<<"-father_sample    "<<result.fSamp<<endl;
        }
        if(result.fInfN){
                cout<<"-father_info      "<<result.fInfo<<endl;
        }
	if(result.phenoN){
		cout<<"-pheno            "<<result.pheno<<endl;
	}
	if(result.minimac || result.phased){
		cout<<"-n_snp            "<<result.nsnp<<endl;
	}
	if(result.chrN){
		cout<<"-chr              "<<result.chr<<endl;
	}
        if(result.outN){
                cout<<"-out              "<<result.outFile<<endl;
        }
	if(result.unordered){
		cout<<"-unordered        true"<<endl;
	}
	if(result.phased){
		cout<<"-phased           true"<<endl;
	}
	if(result.transmitted){
		cout<<"-transmitted      true"<<endl;
	}
	if(result.minimac){
		cout<<"-minimac          true"<<endl;
	}
	if(result.interaction){
		cout<<"-interaction      true"<<endl;
	}
        cout<<endl;
	return;
}

void print_header(){
	printf("\n\nCopyright 2016 Robin Beaumont\nr.beaumont@exeter.ac.uk\n");
	return;
}

void print_required(Parameters result){
        cout<<"\nPlease supply at least the following options:"<<endl;
        if(result.childN==0){
                cout<<"-child            "<<result.childFile<<endl;
        }
        if(result.cSampN==0){
                cout<<"-child_sample     "<<result.cSamp<<endl;
        }
        if(result.motherN==0){
                cout<<"-mother           "<<result.motherFile<<endl;
        }
        if(result.mSampN==0){
                cout<<"-mother_sample    "<<result.mSamp<<endl;
        }
	if(result.chrN==0){
		cout<<"-chr              "<<result.chr<<endl;
	}
        if(result.outN==0){
                cout<<"-out              "<<result.outFile<<endl;
        }
        cout<<endl;
	return;
}
