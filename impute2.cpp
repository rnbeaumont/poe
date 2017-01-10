#include "general_functions.h"
#include "impute2.h"

void process_impute2_probs(string childFile,string motherFile,string cSamp,string mSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample,bool interaction){
	string str;
	vector<string> mfields, cfields;
	unsigned int linec=0,linem=0,cInds,mInds;
	gen cGen(childFile);
	gen mGen(motherFile);
	std::unordered_map<std::string,int> mhash;
	vector<string> cids;
	// get number of SNPs in both mother and child files - need to decrement due to header line in info files. This is still number of SNPs not last index!!!
	cInds=get_num_lines(cSamp);
	cInds-=2;
	cGen.set_num_samples(cInds);
	mInds=get_num_lines(mSamp);
	mInds-=2;
	mGen.set_num_samples(mInds);
	cGen.probs.resize(3*cInds,0);
	mGen.probs.resize(3*mInds,0);
	printf( "processing sample files\n" );
	// check gen file corresponds to the sample file by comparing number of probabiliries in line 1 with number of variants in sample file TODO
	// read sample files into hashes and arrays (I think I need both for both but need to check)
	read_sample_2(cSamp,mSamp,mhash,cids);
	if(sample){
		print_sample(cSamp,outFile);
	}
	//loop through files:
	//open output files
	printf( "opening output files\n" );
	ofstream outpat;
	filtering_ostream outpatboost;
	ofstream outmat;
	filtering_ostream outmatboost;
	ofstream outhet;
	if(!interaction){
		if(transmitted){
			outpat.open(outFile+".untransmitted.gz",std::ios_base::out | std::ios_base::binary);
		}else{
			outpat.open(outFile+".father.gz",std::ios_base::out | std::ios_base::binary);
		}
		outpatboost.push(gzip_compressor());
		outpatboost.push(outpat);
		if(transmitted){
			outmat.open(outFile+".transmitted.gz",std::ios_base::out | std::ios_base::binary);
		}else{
			outmat.open(outFile+".mother.gz",std::ios_base::out | std::ios_base::binary);
		}
		outmatboost.push(gzip_compressor());
		outmatboost.push(outmat);
		outhet.open(outFile+".hets.gz",std::ios_base::out | std::ios_base::binary);
	}else{
		outhet.open(outFile+".interaction.gz",std::ios_base::out | std::ios_base::binary);
	}
	filtering_ostream outhetboost;
	outhetboost.push(gzip_compressor());
	outhetboost.push(outhet);
	//check if length of mfile is less then nsnp and if so replace nsnp by length of mfile
	//loop through mother's file
	//open mother file
	try{
		while(mGen.getline()){
			linem++;
			top:
			if(!cGen.getline()){
				break;
			}
			linec++;
			test:
			if(!mGen.fields[1].compare(cGen.fields[1]) && !mGen.fields[2].compare(cGen.fields[2])){
				mGen.split_fields();
				cGen.split_fields();
				if(interaction){
					process_child_impute2_probs_interaction(outhetboost,cGen.outString(),cInds,cGen.probs,mGen.probs,mhash,cids,chr);
				}else{
					process_child_impute2_probs(outpatboost,outmatboost,outhetboost,cGen.outString(),cInds,cGen.probs,mGen.probs,mhash,cids,chr);
				}
			}else if(stoi(mGen.fields[2])>stoi(cGen.fields[2])){
				if(!cGen.getline()){
					break;
				}
				linec++;
				goto test;
			}else if(stoi(mGen.fields[2])<stoi(cGen.fields[2])){
				if(!mGen.getline()){
					break;
				}
				linem++;
				goto test;
			}else if(!mGen.fields[2].compare(cGen.fields[2])){
				if(!mGen.getline()){
					break;
				}
				linem++;
				goto top;
			}else{
				cerr<<"\n\ninvalid position on line "<<linec<<" of child file or line "<<linem<<" of mother file\n\n"<<endl;
				exit(1);
			}
		}
	}catch(const gzip_error& e){
		cerr<<e.what()<<endl;
	}

	//close output files here
	outhetboost.pop();
	outmatboost.pop();
	outpatboost.pop();
	outhet.close();
	outmat.close();
	outpat.close();
	return;
}

void process_child_impute2_probs(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,string toprint,int cInds,vector<float>& cArray,vector<float>& mArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr){
	string str;
	vector<string> phap;
	vector<string> mhap;
	vector<string> hhap;
	vector<string> fields;
	ofstream pfile;
	int mind;

	//print start of gen file line! TODO
	outpat<<chr<<" "<<toprint;
	outmat<<chr<<" "<<toprint;
	outhet<<chr<<" "<<toprint;
	for(int i=0;i<cInds;i++){
		//check child snp is in mother
		if(!(mhash.find(cids[i])==mhash.end())){
			//find mother index
			mind=mhash[cids[i]];
			//find genotypes
			if(mArray[3*mind]>0.5){ //AA
				if(cArray[3*i]>0.5){ // AA
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else if(cArray[3*i+1]>0.5){ // Aa
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}
			}else if(mArray[3*mind+1]>0.5){ //Aa
				if(cArray[3*i]>0.5){ // AA
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else if(cArray[3*i+2]>0.5){ // aa
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}
			}else if(mArray[3*mind+2]>0.5){ //aa
				if(cArray[3*i+1]>0.5){ // Aa
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else if(cArray[3*i+2]>0.5){ // aa
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else{
			outpat<<" 0 0 0";
			outmat<<" 0 0 0";
			outhet<<" 0 0 0";
		}
	}
	outpat<<"\n";
	outmat<<"\n";
	outhet<<"\n";
	return;
}

void process_impute2_probs_bgen(string childFile,string motherFile,string cSamp,string mSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample,bool interactions){
	string str;
	vector<string> mfields, cfields;
	unsigned int cInds,mInds,cSnp,mSnp;
	std::unordered_map<std::string,int> mhash;
	vector<string> cids;

	bgen cBgen(childFile),mBgen(motherFile);
	cSnp=cBgen.get_m();
	mSnp=mBgen.get_m();
	/* get number of SNPs in both mother and child files - need to decrement due to header line in info files. This is still number of SNPs not last index!!! */
	printf( "processing sample files\n" );
	cInds=get_num_lines(cSamp);
	cInds-=2;
	mInds=get_num_lines(mSamp);
	mInds-=2;
	//check number of individuals in sample file matches number of individuals in bgen file
	if(cInds!=cBgen.get_n()){
		cerr<<"ERROR: "<<childFile<<" contains differing number of individuals to sample file!"<<endl;
		exit(1);
	}
	if(mInds!=mBgen.get_n()){
		cerr<<"ERROR: "<<motherFile<<" contains differing number of individuals to sample file!"<<endl;
		exit(1);
	}
	/* read sample files into hashes and arrays (I think I need both for both but need to check) */
	read_sample_2(cSamp,mSamp,mhash,cids);
	if(sample){
		print_sample(cSamp,outFile);
	}
	//loop through files:
	//open output files
	printf( "opening output files\n" );
	ofstream outpat;
	filtering_ostream outpatboost;
	ofstream outmat;
	filtering_ostream outmatboost;
	ofstream outhet;
	if(!interactions){
		if(transmitted){
			outpat.open(outFile+".untransmitted.gz",std::ios_base::out | std::ios_base::binary);
		}else{
			outpat.open(outFile+".father.gz",std::ios_base::out | std::ios_base::binary);
		}
		outpatboost.push(gzip_compressor());
		outpatboost.push(outpat);
		if(transmitted){
			outmat.open(outFile+".transmitted.gz",std::ios_base::out | std::ios_base::binary);
		}else{
			outmat.open(outFile+".mother.gz",std::ios_base::out | std::ios_base::binary);
		}
		outmatboost.push(gzip_compressor());
		outmatboost.push(outmat);
		outhet.open(outFile+".hets.gz",std::ios_base::out | std::ios_base::binary);
	}else{
		outhet.open(outFile+".interaction.gz",std::ios_base::out | std::ios_base::binary);
	}
	filtering_ostream outhetboost;
	outhetboost.push(gzip_compressor());
	outhetboost.push(outhet);
	//check if length of mfile is less then nsnp and if so replace nsnp by length of mfile
	if(mSnp<nsnp){nsnp=mSnp;}
	//loop through mother's file
	try{
		while(mBgen.cSnp()<mSnp){
			//read the next SNP from mother file
			mBgen.read_variant_id();
			top:
			if(cBgen.cSnp()<cSnp-1){
				cBgen.read_variant_id();
			}else{
				break;
			}
			test:
			if(!mBgen.get_rsid().compare(cBgen.get_rsid()) && mBgen.get_pos()==cBgen.get_pos()){
				mBgen.read_variant_probabilities();
				cBgen.read_variant_probabilities();
				if(interactions){
					process_child_impute2_probs_interaction(outhetboost,cBgen.get_rsid()+" "+(static_cast<ostringstream&>(ostringstream().seekp(0)<<cBgen.get_pos()).str())+" "+cBgen.get_aA()+" "+cBgen.get_aB(),cInds,cBgen.probs,mBgen.probs,mhash,cids,chr);
				}else{
					process_child_impute2_probs(outpatboost,outmatboost,outhetboost,cBgen.get_rsid()+" "+(static_cast<ostringstream&>(ostringstream().seekp(0)<<cBgen.get_pos()).str())+" "+cBgen.get_aA()+" "+cBgen.get_aB(),cInds,cBgen.probs,mBgen.probs,mhash,cids,chr);
				}
			}else if(mBgen.get_pos()>cBgen.get_pos()){
				if(cBgen.cSnp()<cSnp-1){
					//skip rest of line in cfile
					cBgen.skip_variant_probabilities();
					//read next snp header
					cBgen.read_variant_id();
				}else{
					break;
				}
				goto test;
			}else if(mBgen.get_pos()<cBgen.get_pos()){
				if(mBgen.cSnp()<mSnp-1){
					//skip rest of line in mfile
					mBgen.skip_variant_probabilities();
					//read next snp header
					mBgen.read_variant_id();
				}else{
					break;
				}
				goto test;
			}else if(mBgen.get_pos()==cBgen.get_pos()){
				if(mBgen.cSnp()<mSnp-1){
					//skip rest of line in mfile
					mBgen.skip_variant_probabilities();
					//read next snp header
					mBgen.read_variant_id();
				}else{
					break;
				}
				goto top;
			}else{
				cerr<<"\nError reading SNP "<<mBgen.get_rsid()<<". Exiting"<<endl;
				exit(1);
			}
		}
	}catch(const gzip_error& e){
		cerr<<e.what()<<endl;
	}

	//close output files here
	outhetboost.pop();
	if(!interactions){
		outmatboost.pop();
		outpatboost.pop();
	}
	outhet.close();
	if(!interactions){
		outmat.close();
		outpat.close();
	}
	return;
}

void process_child_impute2_probs_interaction(filtering_ostream& outhet,string toprint,int cInds,vector<float>& cArray,vector<float>& mArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr){
	string str;
	vector<string> phap;
	vector<string> mhap;
	vector<string> hhap;
	vector<string> fields;
	ofstream pfile;
	int mind;

	//print start of gen file line! TODO
	outhet<<chr<<" "<<toprint;
	for(int i=0;i<cInds;i++){
		//check child snp is in mother
		if(!(mhash.find(cids[i])==mhash.end())){
			//find mother index
			mind=mhash[cids[i]];
			//find genotypes
			if(mArray[3*mind]>mArray[3*mind+1] && mArray[3*mind]>mArray[3*mind+2]){ //AA
				if(cArray[3*i]>cArray[3*i+1] && cArray[3*i]>cArray[3*i+2]){ // AA
					outhet<<" 1 0 0";
				}else{
					outhet<<" 0 1 0";
				}
			}else if(mArray[3*mind+1]>mArray[3*mind] && mArray[3*mind+1]>mArray[3*mind+2]){ //Aa
				if(cArray[3*i+1]>cArray[3*i] && cArray[3*i+1]>cArray[3*i+2]){ // Aa
					outhet<<" 1 0 0";
				}else{
					outhet<<" 0 1 0";
				}
			}else if(mArray[3*mind+2]>mArray[3*mind] && mArray[3*mind+2]>mArray[3*mind+1]){ //aa
				if(cArray[3*i+2]>cArray[3*i] && cArray[3*i+2]>cArray[3*i+1]){ // aa
					outhet<<" 1 0 0";
				}else{
					outhet<<" 0 1 0";
				}
			}else{
				outhet<<" 0 0 0";
			}
		}else{
			outhet<<" 0 0 0";
		//	cout<<"error2"<<endl;
		}
	}
	outhet<<"\n";
	return;
}
