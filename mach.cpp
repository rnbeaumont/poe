#include "general_functions.h"
#include "mach.h"

void process_mach(string childFile,string motherFile,string cInfo,string mInfo,string outFile,string pheno,string chr,bool transmitted,unsigned int nsnp,bool interactions){
	string str;
	vector<string> fields;
	int loopE;
	unsigned int cInds,mInds,cSnp,mSnp;
	unordered_map<string,int> chash;
	vector<string> mids, mstr;
	unordered_map<string, string> pMap;
	vector<vector<float>> mHapArray;
	unordered_map<string,int> mline;
	string pHead;
	/* get number of individuals in both mother and child files */
	cInds=get_num_lines(childFile);
	mInds=get_num_lines(motherFile);
	/* get number of SNPs in both mother and child files - need to decrement due to header line in info files. This is still number of SNPs not last index!!! */
	cSnp=get_num_lines(cInfo);
	cSnp--;
	mSnp=get_num_lines(mInfo);
	mSnp--;
	cout<<"processing info files"<<endl;
	// check dosage file corresponds to the info file by comparing number of dosages in line 1 with number of variants in info file
	check_dosage_length(childFile,cSnp+2);
	check_dosage_length(motherFile,mSnp+2);
	/* read info files into hashes and arrays (I think I need both for both but need to check) */
//	cout<<"processing child's info file\n";
//	read_info(cInfo,cSnp,chash,cids);
//	cout<<"processing mother's info file\n";
//	read_info(mInfo,mSnp,mhash,mids);
	read_info_2(cInfo,mInfo,chash,mids,mstr,chr);
	//read in phenotype file
	cout<<"processing phenotype file"<<endl;
	read_pheno(pheno,pMap,pHead);
	//loop through files:
	//open output files
	cout<<"opening output files"<<endl;
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
	//allocate array to hold genotypes/dosages (deallocate if already allocated!)
	mHapArray.resize(mInds,vector<float>(nsnp,0));
	//loop through mother's file
	for(unsigned int i=0;i<mSnp;i=i+nsnp){
		cout<<"processing mother's genotypes"<<endl;
		//open mother's file
		ifstream mFile(motherFile,ios::in | ios_base::binary);
		try{
			filtering_istream boost_in;
			boost_in.push(gzip_decompressor());
			boost_in.push(mFile);
			//read dosages into array and create unordered_map holding location of each mother
			int rCount=0;
			while(getline(boost_in,str)){
				split(fields,str,is_any_of(" \t"));
				mline[fields[0]]=rCount;
				if(i+nsnp<mSnp){
					loopE=nsnp;
				}else{
					loopE=mSnp-i;
				}
				for(int j=0;j<loopE;j++){
					string buf2=fields[i+j+2];
					float buf=stof(buf2);
					mHapArray[rCount][j]=buf;
				}
				rCount++;
			}
		}catch(const gzip_error& e){
			cout<<e.what()<<endl;
		}
		mFile.close();
		//process child's file
		if((i+nsnp)<mSnp){
			cout<<"processing child's genotypes"<<endl;
			if(interactions){
				process_child_mach_interaction(outhetboost,i,nsnp,cInds,chash,mids,mstr,pMap,mHapArray,mline,pHead,outFile,childFile);
			}else{
				process_child_mach(outpatboost,outmatboost,outhetboost,i,nsnp,cInds,chash,mids,mstr,pMap,mHapArray,mline,pHead,outFile,childFile);
			}
		}else{
			cout<<"processing child's genotypes for the last time"<<endl;
			if(interactions){
				process_child_mach_interaction(outhetboost,i,nsnp,cInds,chash,mids,mstr,pMap,mHapArray,mline,pHead,outFile,childFile);
			}else{
				process_child_mach(outpatboost,outmatboost,outhetboost,i,mSnp-i-1,cInds,chash,mids,mstr,pMap,mHapArray,mline,pHead,outFile,childFile);
			}
		}
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

void process_child_mach(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,int loopStart,int loopEnd,int cInds,unordered_map<string,int>& chash,vector<string>& mids,vector<string>& mstr,unordered_map<string, string>& pMap,vector<vector<float>>& mHapArray,unordered_map<string,int>& mline,string pHead,string outFile,string childFile){
	string str;
	vector<string> phap;
	vector<string> mhap;
	vector<string> hhap;
	vector<string> fields;
	vector<string> cnames;
	int** pharr;
	int** mharr;
	int** hharr;
	int cCount=0;
	ofstream pfile;
	pharr=new int*[cInds];
	mharr=new int*[cInds];
	hharr=new int*[cInds];
	for(int k=0;k<cInds;++k){
		pharr[k]=new int[loopEnd];
		mharr[k]=new int[loopEnd];
		hharr[k]=new int[loopEnd];
	}
	//allocate array for the three 
	//print phenotype file header
	if(loopStart==0){
		pfile.open(outFile+".pheno");
		pfile << "ID_1 ID_2 missing "<<pHead<<"\n0 0 0 P\n";
	}
	//loop through genotypes being tested here (first nsnp mother genotypes)- if they're in both files then initialise strings to hold them
	for(int i=0;i<loopEnd;i++){
		phap.push_back("");
		mhap.push_back("");
		hhap.push_back("");
		//if(chash.count(mids[i+loopStart])){
		if(!(chash.find(mids[i+loopStart]) == chash.end())){
			phap[i].append(mstr[i+loopStart]);
			mhap[i].append(mstr[i+loopStart]);
			hhap[i].append(mstr[i+loopStart]);
		}
	}
	//open child file
	ifstream file(childFile,ios_base::in | ios_base::binary);
	try{
		filtering_istream boost_in;
		boost_in.push(gzip_decompressor());
		boost_in.push(file);
		//loop through individuals in the file
		while(getline(boost_in, str)){
			//read line into buffer
			split(fields,str,is_any_of(" \t"));
			//check if mother exists
			cnames.push_back(fields[0]);
			if(!(mline.find(fields[0]) == mline.end())){
				//if this is the first loop then print phenotype
				if(loopStart==0){
					if(!(pMap.find(fields[0])==pMap.end())){
						pfile <<fields[0]<<" "<<fields[0]<<" 0 "<<pMap[fields[0]]<<"\n";
					}else{
						pfile <<fields[0]<<" "<<fields[0]<<" 0 NA\n";
					}
				}
				int mind=mline[fields[0]];
				//loop through 0 to nsnp
				for(int i=0;i<loopEnd;i++){
					//check if mother SNP is in child file
					if(!(chash.find(mids[i+loopStart])==chash.end())){
						//if so do stuff
						int pos=chash[mids[i+loopStart]]+2;
						float cprob=stof(fields[pos]);
						if( mHapArray[mind][i] <= 0.5 ){ //AA
							if( cprob <= 0.5 ){	//	AA
								pharr[cCount][i]=0;
								mharr[cCount][i]=0;
								hharr[cCount][i]=-1;
							}else if(cprob>0.5 && cprob<1.5){	//	Aa
								pharr[cCount][i]=1;
								mharr[cCount][i]=0;
								hharr[cCount][i]=0;
							}else{
								pharr[cCount][i]=-1;
								mharr[cCount][i]=-1;
								hharr[cCount][i]=-1;
							}
						}else if(mHapArray[mind][i]>0.5 && mHapArray[mind][i]<1.5){ //Aa
							if( cprob <= 0.5 ){	//	AA
								pharr[cCount][i]=0;
								mharr[cCount][i]=0;
								hharr[cCount][i]=-1;
							}else if(cprob>=1.5){	//	aa
								pharr[cCount][i]=1;
								mharr[cCount][i]=1;
								hharr[cCount][i]=-1;
							}else{
								pharr[cCount][i]=-1;
								mharr[cCount][i]=-1;
								hharr[cCount][i]=-1;
							}
						}else if(mHapArray[mind][i]>=1.5){ //aa
							if(cprob>0.5 && cprob<1.5){	//	Aa
								pharr[cCount][i]=0;
								mharr[cCount][i]=1;
								hharr[cCount][i]=1;
							}else if(cprob>=1.5){	//	aa
								pharr[cCount][i]=1;
								mharr[cCount][i]=1;
								hharr[cCount][i]=-1;
							}else{
								pharr[cCount][i]=-1;
								mharr[cCount][i]=-1;
								hharr[cCount][i]=-1;
							}
						}else{
							pharr[cCount][i]=-1;
							mharr[cCount][i]=-1;
							hharr[cCount][i]=-1;
						}
					}
				}
			}
			cCount++;
		}
		//print the haplotypes to the output files
		for(int i=0;i<loopEnd;i++){
			if(!(chash.find(mids[i+loopStart])==chash.end())){
				outpat<<phap[i];
				outmat<<mhap[i];
				outhet<<hhap[i];
				for(int k=0;k<cCount;k++){
					if(!(mline.find(cnames[k]) == mline.end())){
						if(pharr[k][i]==1){
							outpat<<" 0 1 0";
						}else if(pharr[k][i]==0){
							outpat<<" 1 0 0";
						}else{
							outpat<<" 0 0 0";
						}
						if(mharr[k][i]==1){
							outmat<<" 0 1 0";
						}else if(mharr[k][i]==0){
							outmat<<" 1 0 0";
						}else{
							outmat<<" 0 0 0";
						}
						if(hharr[k][i]==1){
							outhet<<" 0 1 0";
						}else if(hharr[k][i]==0){
							outhet<<" 1 0 0";
						}else{
							outhet<<" 0 0 0";
						}
					}
				}
				outpat<<"\n";
				outmat<<"\n";
				outhet<<"\n";
			}
		}
	}catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	//close the child file!
	file.close();
	if(loopStart==0){
		pfile.close();
	}
	return;
}

void process_child_mach_interaction(filtering_ostream& outhet,int loopStart,int loopEnd,int cInds,unordered_map<string,int>& chash,vector<string>& mids,vector<string>& mstr,unordered_map<string, string>& pMap,vector<vector<float>>& mHapArray,unordered_map<string,int>& mline,string pHead,string outFile,string childFile){
	string str;
	vector<string> hhap;
	vector<string> fields;
	vector<string> cnames;
	int** hharr;
	int cCount=0;
	ofstream pfile;
	hharr=new int*[cInds];
	for(int k=0;k<cInds;++k){
		hharr[k]=new int[loopEnd];
	}
	//allocate array for the three 
	//print phenotype file header
	if(loopStart==0){
		pfile.open(outFile+".pheno");
		pfile << "ID_1 ID_2 missing "<<pHead<<"\n0 0 0 P\n";
	}
	//loop through genotypes being tested here (first nsnp mother genotypes)- if they're in both files then initialise strings to hold them
	for(int i=0;i<loopEnd;i++){
		hhap.push_back("");
		//if(chash.count(mids[i+loopStart])){
		if(!(chash.find(mids[i+loopStart]) == chash.end())){
			hhap[i].append(mstr[i+loopStart]);
		}
	}
	//open child file
	ifstream file(childFile,ios_base::in | ios_base::binary);
	try{
		filtering_istream boost_in;
		boost_in.push(gzip_decompressor());
		boost_in.push(file);
		//loop through individuals in the file
		while(getline(boost_in, str)){
			//read line into buffer
			split(fields,str,is_any_of(" \t"));
			//check if mother exists
			cnames.push_back(fields[0]);
			if(!(mline.find(fields[0]) == mline.end())){
				//if this is the first loop then print phenotype
				if(loopStart==0){
					if(!(pMap.find(fields[0])==pMap.end())){
						pfile <<fields[0]<<" "<<fields[0]<<" 0 "<<pMap[fields[0]]<<"\n";
					}else{
						pfile <<fields[0]<<" "<<fields[0]<<" 0 NA\n";
					}
				}
				int mind=mline[fields[0]];
				//loop through 0 to nsnp
				for(int i=0;i<loopEnd;i++){
					//check if mother SNP is in child file
					if(!(chash.find(mids[i+loopStart])==chash.end())){
						//if so do stuff
						int pos=chash[mids[i+loopStart]]+2;
						float cprob=stof(fields[pos]);
						if( mHapArray[mind][i] <= 0.5 ){ //AA
							if( cprob <= 0.5 ){	//	AA
								hharr[cCount][i]=0;
							}else if(cprob>0.5 && cprob<1.5){	//	Aa
								hharr[cCount][i]=1;
							}else{
								hharr[cCount][i]=-1;
							}
						}else if(mHapArray[mind][i]>0.5 && mHapArray[mind][i]<1.5){ //Aa
							if( cprob <= 0.5 ){	//	AA
								hharr[cCount][i]=1;
							}else if(cprob>0.5 && cprob<1.5){	//	Aa
								hharr[cCount][i]=0;
							}else if(cprob>=1.5){	//	aa
								hharr[cCount][i]=1;
							}else{
								hharr[cCount][i]=-1;
							}
						}else if(mHapArray[mind][i]>=1.5){ //aa
							if(cprob>0.5 && cprob<1.5){	//	Aa
								hharr[cCount][i]=1;
							}else if(cprob>=1.5){	//	aa
								hharr[cCount][i]=0;
							}else{
								hharr[cCount][i]=-1;
							}
						}else{
							hharr[cCount][i]=-1;
						}
					}
				}
			}
			cCount++;
		}
		//print the haplotypes to the output files
		for(int i=0;i<loopEnd;i++){
			if(!(chash.find(mids[i+loopStart])==chash.end())){
				outhet<<hhap[i];
				for(int k=0;k<cCount;k++){
					if(!(mline.find(cnames[k]) == mline.end())){
						if(hharr[k][i]==1){
							outhet<<" 0 1 0";
						}else if(hharr[k][i]==0){
							outhet<<" 1 0 0";
						}else{
							outhet<<" 0 0 0";
						}
					}
				}
				outhet<<"\n";
			}
		}
	}catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	//close the child file!
	file.close();
	if(loopStart==0){
		pfile.close();
	}
	return;
}

