#include "general_functions.h"
#include "impute2_trio.h"

void process_impute2_trio_probs(string childFile,string motherFile,string fatherFile,string cSamp,string mSamp,string fSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample){
	string str;
	vector<string> mfields, cfields, ffields;
	unsigned int linec=0,linem=0,linef=0,cInds,mInds,fInds;
	gen cGen(childFile);
	gen mGen(motherFile);
	gen fGen(fatherFile);
	std::unordered_map<std::string,int> mhash;
	std::unordered_map<std::string,int> fhash;
	vector<string> cids;
	/* get number of SNPs in both mother and child files - need to decrement due to header line in info files. This is still number of SNPs not last index!!! */
	cInds=get_num_lines(cSamp);
	cInds-=2;
	cGen.set_num_samples(cInds);
	mInds=get_num_lines(mSamp);
	mInds-=2;
	mGen.set_num_samples(mInds);
	fInds=get_num_lines(fSamp);
	fInds-=2;
	fGen.set_num_samples(fInds);
	cGen.probs.resize(3*cInds,0);
	mGen.probs.resize(3*mInds,0);
	fGen.probs.resize(3*fInds,0);
	printf( "processing sample files\n" );
	// check gen file corresponds to the sample file by comparing number of probabilities in line 1 with number of variants in sample file TODO
	/* read sample files into hashes and arrays (I think I need both for both but need to check) */
	read_sample_3(cSamp,mSamp,fSamp,mhash,fhash,cids);
	if(sample){
		print_sample(cSamp,outFile);
	}
	//loop through files:
	//open output files
	printf( "opening output files\n" );
	ofstream outpat;
	if(transmitted){
		outpat.open(outFile+".untransmitted.gz",std::ios_base::out | std::ios_base::binary);
	}else{
		outpat.open(outFile+".father.gz",std::ios_base::out | std::ios_base::binary);
	}
	filtering_ostream outpatboost;
	outpatboost.push(gzip_compressor());
	outpatboost.push(outpat);
	ofstream outmat;
	if(transmitted){
		outmat.open(outFile+".transmitted.gz",std::ios_base::out | std::ios_base::binary);
	}else{
		outmat.open(outFile+".mother.gz",std::ios_base::out | std::ios_base::binary);
	}
	filtering_ostream outmatboost;
	outmatboost.push(gzip_compressor());
	outmatboost.push(outmat);
	ofstream outhet(outFile+".hets.gz",std::ios_base::out | std::ios_base::binary);
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
			if(!fGen.getline()){
				break;
			}
			linef++;
			test:
			if(!mGen.fields[1].compare(cGen.fields[1]) && !mGen.fields[2].compare(cGen.fields[2]) && !mGen.fields[1].compare(fGen.fields[1]) && !mGen.fields[2].compare(fGen.fields[2])){
				mGen.split_fields();
				cGen.split_fields();
				fGen.split_fields();
				process_child_impute2_trio_probs(outpatboost,outmatboost,outhetboost,cGen.outString(),cInds,cGen.probs,mGen.probs,fGen.probs,mhash,fhash,cids,chr,transmitted);
			}else if(stoi(mGen.fields[2])>stoi(cGen.fields[2])){
				if(!cGen.getline()){
					break;
				}
				linec++;
				if(stoi(mGen.fields[2])>stoi(fGen.fields[2])){
					if(!fGen.getline()){
						break;
					}
					linef++;
				}
				goto test;
			}else if(stoi(mGen.fields[2])<stoi(cGen.fields[2])){	//add father from here
				if(!mGen.getline()){
					break;
				}
				linem++;
				if(stoi(fGen.fields[2])<stoi(cGen.fields[2])){
					if(!fGen.getline()){
						break;
					}
					linef++;
				}
				goto test;
			}else if(stoi(fGen.fields[2])>stoi(cGen.fields[2])){
				if(!cGen.getline()){
					break;
				}
				linec++;
				if(stoi(fGen.fields[2])>stoi(mGen.fields[2])){
					if(!mGen.getline()){
						break;
					}
					linem++;
				}
				goto test;
			}else if(stoi(mGen.fields[2])>stoi(fGen.fields[2])){
				if(!fGen.getline()){
					break;
				}
				linef++;
				goto test;
			}else if(!mGen.fields[2].compare(cGen.fields[2])){
				if(!mGen.getline()){
					break;
				}
				linem++;
				goto top;
			}else{
				cerr<<"\n\ninvalid position on line "<<linec<<" of child file or line "<<linem<<" of mother file or line "<<linef<<" of father file\n\n"<<endl;
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

void process_child_impute2_trio_probs(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,string toprint,int cInds,vector<float>& cArray,vector<float>& mArray,vector<float>& fArray,std::unordered_map<std::string,int>& mhash,std::unordered_map<std::string,int>& fhash,vector<string>& cids,string chr,bool transmitted){
	int mind;
	int find;
	//print start of gen file line!
	outpat<<chr<<" "<<toprint;
	outmat<<chr<<" "<<toprint;
	outhet<<chr<<" "<<toprint;
	for(int i=0;i<cInds;i++){
		//check child snp is in mother
		if(!(mhash.find(cids[i])==mhash.end()) && !(fhash.find(cids[i])==fhash.end())){
			//find mother index
			mind=mhash[cids[i]];
			find=fhash[cids[i]];
			//find genotypes
			compare_trio_genotypes_unphased(cArray,mArray,fArray,i,mind,find,outpat,outmat,outhet,transmitted);
		}else if(!(mhash.find(cids[i])==mhash.end())){
			mind=mhash[cids[i]];
			compare_pair_genotypes_unphased(cArray,mArray,i,mind,outpat,outmat,outhet,transmitted,true);
		}else if(!(fhash.find(cids[i])==fhash.end())){
			find=fhash[cids[i]];
			compare_pair_genotypes_unphased(cArray,fArray,i,find,outpat,outmat,outhet,transmitted,false);
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

void compare_trio_genotypes_unphased(vector<float>& cArray,vector<float>& mArray,vector<float>& fArray,int i,int mind,int find,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted){
	if(mArray[3*mind]>0.5){ //AA
		if(fArray[3*find]>0.5){ //AA
			if(cArray[3*i]>0.5){ // AA      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if(fArray[3*find+1]>0.5){ //Aa
			if(cArray[3*i]>0.5){ // AA      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
			}else if(cArray[3*i+1]>0.5){ // Aa
				if(transmitted){
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if(fArray[3*find+2]>0.5){ //aa
			if(cArray[3*i+1]>0.5){ // Aa
				if(transmitted){
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
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
	}else if(mArray[3*mind+1]>0.5){ //Aa
		if(fArray[3*find]>0.5){ //AA
			if(cArray[3*i]>0.5){ // AA
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else if(cArray[3*i+1]>0.5){ // Aa	//transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 1 0 0";
				outmat<<" 0 1 0";
				outhet<<" 0 1 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if(fArray[3*find+1]>0.5){ //Aa
			if(cArray[3*i]>0.5){ // AA
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else if(cArray[3*i+2]>0.5){ // aa
				if(transmitted){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if(fArray[3*find+2]>0.5){ //aa
			if(cArray[3*i+1]>0.5){ // Aa      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 0 1 0";
				outmat<<" 1 0 0";
				outhet<<" 1 0 0";
			}else if(cArray[3*i+2]>0.5){ // aa
				if(transmitted){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
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
	}else if(mArray[3*mind+2]>0.5){ //aa
		if(fArray[3*find]>0.5){ //AA
			if(cArray[3*i+1]>0.5){ // Aa
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if(fArray[3*find+1]>0.5){ //Aa
			if(cArray[3*i+1]>0.5){ // Aa
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}
			}else if(cArray[3*i+1]>0.5){ // aa      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 0 1 0";
				outmat<<" 0 1 0";
				outhet<<" 0 0 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if(fArray[3*find+2]>0.5){ //aa
			if(cArray[3*i+1]>0.5){ // aa      //transmitted and inherited here are the same so don't need to test which one we want
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

void compare_pair_genotypes_unphased(vector<float>& cArray,vector<float>& mArray,int i,int mind,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted,bool mother){
	if(mArray[3*mind]>0.5){ //AA
		if(cArray[3*i]>0.5){ //AA
			if(mother || !transmitted){
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
			}
		}else if(cArray[3*i+1]>0.5){ // Aa
			if(mother){
				if(transmitted){
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}
			}else{
				if(transmitted){
					outpat<<" 0 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}
			}
		}else{
			outpat<<" 0 0 0";
			outmat<<" 0 0 0";
			outhet<<" 0 0 0";
		}
	}else if(mArray[3*mind+1]>0.5){ //Aa
		if(cArray[3*i]>0.5){ // AA
			if(mother){
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else{
				if(transmitted){
					outpat<<" 0 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}
		}else if(cArray[3*i+2]>0.5){ // aa
			if(mother){
				if(transmitted){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}
			}else{
				if(transmitted){
					outpat<<" 0 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}
			}
		}else{
			outpat<<" 0 0 0";
			outmat<<" 0 0 0";
			outhet<<" 0 0 0";
		}
	}else if(mArray[3*mind+2]>0.5){ //aa
		if(cArray[3*i+1]>0.5){ // Aa
			if(mother){
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}
			}else{
				if(transmitted){
					outpat<<" 0 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}
			}
		}else if(cArray[3*i+1]>0.5){ // aa
			if(mother || !transmitted){
				outpat<<" 0 1 0";
				outmat<<" 0 1 0";
				outhet<<" 0 0 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 1 0";
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

void process_impute2_trio_probs_bgen(string childFile,string motherFile,string fatherFile,string cSamp,string mSamp,string fSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample){
	string str;
	vector<string> mfields, cfields, ffields;
	unsigned int linec=0,linem=0,linef=0,cInds,mInds,fInds,cSnp,mSnp,fSnp;
	bgen cBgen(childFile),mBgen(motherFile),fBgen(fatherFile);
	cSnp=cBgen.get_m();
	mSnp=mBgen.get_m();
	fSnp=fBgen.get_m();
	std::unordered_map<std::string,int> mhash;
	std::unordered_map<std::string,int> fhash;
	vector<string> cids;
	/* get number of SNPs in both mother and child files - need to decrement due to header line in info files. This is still number of SNPs not last index!!! */
	printf( "processing sample files\n" );
	cInds=get_num_lines(cSamp);
	cInds-=2;
	mInds=get_num_lines(mSamp);
	mInds-=2;
	fInds=get_num_lines(fSamp);
	fInds-=2;

	//check number of individuals in sample file matches number of individuals in bgen file
	if(cInds!=cBgen.get_n()){
		cerr<<"ERROR: "<<childFile<<" contains differing number of individuals to sample file!"<<endl;
		exit(1);
	}
	if(mInds!=mBgen.get_n()){
		cerr<<"ERROR: "<<motherFile<<" contains differing number of individuals to sample file!"<<endl;
		exit(1);
	}
	if(fInds!=fBgen.get_n()){
		cerr<<"ERROR: "<<fatherFile<<" contains differing number of individuals to sample file!"<<endl;
		exit(1);
	}
	/* read sample files into hashes and arrays (I think I need both for both but need to check) */
	read_sample_3(cSamp,mSamp,fSamp,mhash,fhash,cids);
	if(sample){
		print_sample(cSamp,outFile);
	}
	//loop through files:
	//open output files
	printf( "opening output files\n" );
	ofstream outpat;
	if(transmitted){
		outpat.open(outFile+".untransmitted.gz",std::ios_base::out | std::ios_base::binary);
	}else{
		outpat.open(outFile+".father.gz",std::ios_base::out | std::ios_base::binary);
	}
	filtering_ostream outpatboost;
	outpatboost.push(gzip_compressor());
	outpatboost.push(outpat);
	ofstream outmat;
	if(transmitted){
		outmat.open(outFile+".transmitted.gz",std::ios_base::out | std::ios_base::binary);
	}else{
		outmat.open(outFile+".mother.gz",std::ios_base::out | std::ios_base::binary);
	}
	filtering_ostream outmatboost;
	outmatboost.push(gzip_compressor());
	outmatboost.push(outmat);
	ofstream outhet(outFile+".hets.gz",std::ios_base::out | std::ios_base::binary);
	filtering_ostream outhetboost;
	outhetboost.push(gzip_compressor());
	outhetboost.push(outhet);
	//check if length of mfile is less then nsnp and if so replace nsnp by length of mfile
	//loop through mother's file
	//open mother file
	try{
		while(mBgen.cSnp()<mSnp){
			mBgen.read_variant_id();
			top:
			if(cBgen.cSnp()<cSnp-1){
				cBgen.read_variant_id();
			}else{
				break;
			}
			if(fBgen.cSnp()<fSnp-1){
				fBgen.read_variant_id();
			}else{
				break;
			}
			test:
			if(!mBgen.get_rsid().compare(cBgen.get_rsid()) && mBgen.get_pos()==cBgen.get_pos() && !mBgen.get_rsid().compare(fBgen.get_rsid()) && mBgen.get_pos()==fBgen.get_pos()){
				mBgen.read_variant_probabilities();
				cBgen.read_variant_probabilities();
				fBgen.read_variant_probabilities();
				process_child_impute2_trio_probs(outpatboost,outmatboost,outhetboost,cBgen.get_rsid()+" "+(static_cast<ostringstream&>(ostringstream().seekp(0)<<cBgen.get_pos()).str())+" "+cBgen.get_aA()+" "+cBgen.get_aB(),cInds,cBgen.probs,mBgen.probs,fBgen.probs,mhash,fhash,cids,chr,transmitted);
			}else if(mBgen.get_pos()>cBgen.get_pos()){
				if(cBgen.cSnp()<cSnp-1){
					//skip rest of line in cfile
					cBgen.skip_variant_probabilities();
					//read next snp header
					cBgen.read_variant_id();
				}else{
					break;
				}
				if(mBgen.get_pos()>fBgen.get_pos()){
					if(fBgen.cSnp()<fSnp-1){
						//skip rest of line in ffile
						fBgen.skip_variant_probabilities();
						//read next snp header
						fBgen.read_variant_id();
					}else{
						break;
					}
				}
				goto test;
			}else if(mBgen.get_pos()<cBgen.get_pos()){	//add father from here
				if(mBgen.cSnp()<mSnp-1){
					//skip rest of line in cfile
					mBgen.skip_variant_probabilities();
					//read next snp header
					mBgen.read_variant_id();
				}else{
					break;
				}
				if(fBgen.get_pos()<cBgen.get_pos()){
					if(fBgen.cSnp()<fSnp-1){
						//skip rest of line in ffile
						fBgen.skip_variant_probabilities();
						//read next snp header
						fBgen.read_variant_id();
					}else{
						break;
					}
				}
				goto test;
			}else if(fBgen.get_pos()>cBgen.get_pos()){
				if(cBgen.cSnp()<cSnp-1){
					//skip rest of line in cfile
					cBgen.skip_variant_probabilities();
					//read next snp header
					cBgen.read_variant_id();
				}else{
					break;
				}
				if(fBgen.get_pos()>mBgen.get_pos()){
					if(mBgen.cSnp()<mSnp-1){
						//skip rest of line in cfile
						mBgen.skip_variant_probabilities();
						//read next snp header
						mBgen.read_variant_id();
					}else{
						break;
					}
				}
				goto test;
			}else if(mBgen.get_pos()>fBgen.get_pos()){
				if(fBgen.cSnp()<fSnp-1){
					//skip rest of line in ffile
					fBgen.skip_variant_probabilities();
					//read next snp header
					fBgen.read_variant_id();
				}else{
					break;
				}
				goto test;
			}else if(mBgen.get_pos()==cBgen.get_pos()){
				if(mBgen.cSnp()<mSnp-1){
					//skip rest of line in cfile
					mBgen.skip_variant_probabilities();
					//read next snp header
					mBgen.read_variant_id();
				}else{
					break;
				}
				linem++;
				goto top;
			}else{
				cerr<<"\n\ninvalid position on line "<<linec<<" of child file or line "<<linem<<" of mother file or line "<<linef<<" of father file\n\n"<<endl;
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
