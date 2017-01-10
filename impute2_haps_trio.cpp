#include "general_functions.h"
#include "impute2_haps_trio.h"

void process_impute2_haps_trio(string childFile,string motherFile,string fatherFile,string cSamp,string mSamp,string fSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample){
	string str;
	vector<string> mfields,ffields, cfields;
	unsigned int linec=0,linem=0,linef=0;
	int nel=0;
	bool last=0;
	unsigned int mSnp,cInds,mInds,fInds;
	gen cGen(childFile);
	gen mGen(motherFile);
	gen fGen(fatherFile);
	vector<string> titles;
	unordered_map<string,int> mhash;
	unordered_map<string,int> fhash;
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
	cGen.haps.resize(nsnp,vector<float>(cInds*2,0));
	mGen.haps.resize(nsnp,vector<float>(mInds*2,0));
	fGen.haps.resize(nsnp,vector<float>(fInds*2,0));
	cout<<"Processing Sample files"<<endl;
	/* read sample files into hashes and arrays (I think I need both for both but need to check) */
	read_sample_3(cSamp,mSamp,fSamp,mhash,fhash,cids);
        if(sample){
                print_sample(cSamp,outFile);
        }
	//loop through files:
	//open output files
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
	mGen.read_num_lines();
	mSnp=mGen.get_num_lines();
	if(mSnp<nsnp){nsnp=mSnp;}
	//loop through mother's file
	//open mother file
	cout<<"Processing Genotypes"<<endl;
	try{
		// load in the first n_snp snps from both files, but need to align them 
		loopStart:
		nel=0;
		titles.clear();
		for(unsigned int z=0;z<nsnp;++z){		//TODO loop index needs to be updated for the last loop!!! (maybe, check it's OK!!!!)
			top1:
			if(mGen.getline()){
				linem++;
			}else{
				last=1;
				goto endc;
			}
			if(cGen.getline()){
				linec++;
			}else{
				last=1;
				goto endc;
			}
			if(fGen.getline()){
				linef++;
			}else{
				last=1;
				goto endc;
			}
			first:
			if(!mGen.fields[1].compare(cGen.fields[1]) && !mGen.fields[2].compare(cGen.fields[2]) && !mGen.fields[1].compare(fGen.fields[1]) && !mGen.fields[2].compare(fGen.fields[2])){
				nel++;
				mGen.split_haps(z);
				fGen.split_haps(z);
				cGen.split_haps(z);
				titles.push_back(cGen.fields[1]+" "+cGen.fields[2]+" "+cGen.fields[3]+" "+cGen.fields[4]);
			}else if(stoi(mGen.fields[2])>stoi(cGen.fields[2])){
				if(cGen.getline()){
					linec++;
				}else{
					last=1;
					goto endc;
				}
				if(stoi(mGen.fields[2])>stoi(fGen.fields[2])){
					if(fGen.getline()){
						linef++;
					}else{
						last=1;
						goto endc;
					}
				}
				goto first;
			}else if(stoi(mGen.fields[2])<stoi(cGen.fields[2])){
				if(!mGen.getline()){
					last=1;
					goto endc;
				}
				linem++;
				if(stoi(fGen.fields[2])<stoi(cGen.fields[2])){
					if(!fGen.getline()){
						last=1;
						goto endc;
					}
					linef++;
				}
				goto first;
			}else if(stoi(fGen.fields[2])>stoi(cGen.fields[2])){
				if(!cGen.getline()){
					last=1;
					goto endc;
				}
				linec++;
				if(stoi(fGen.fields[2])>stoi(mGen.fields[2])){
					if(!mGen.getline()){
						last=1;
						goto endc;
					}
					linem++;
				}
				goto first;
			}else if(stoi(mGen.fields[2])>stoi(fGen.fields[2])){
				if(!fGen.getline()){
					last=1;
					goto endc;
				}
				linef++;
				goto first;
			}else if(!mGen.fields[2].compare(cGen.fields[2])){
				if(!mGen.getline()){
					last=1;
					goto endc;
				}
				linem++;
				goto top1;
			}else{
				cerr<<"\n\ninvalid position on line "<<linec<<" of child file or line "<<linem<<" of mother file\n\n"<<endl;
				exit(1);
			}
		}
		endc:
		process_child_impute2_haps_trio(outpatboost,outmatboost,outhetboost,titles,nel,cInds,cGen.haps,mGen.haps,fGen.haps,mhash,fhash,cids,chr,transmitted);	//nel will tell process_child how many elements are in the array, should be the same as titles.size();
		if(!last){
			goto loopStart;
		}
//		process_child_impute2_haps_trio(outpatboost,outmatboost,outhetboost,titles,nel,cInds,cHapArray,mHapArray,fHapArray,mhash,fhash,cids,chr,transmitted);	//nel will tell process_child how many elements are in the array, should be the same as titles.size();
	}catch(const gzip_error& e){
		cerr<<e.what()<<endl;
	}
	//
	//close output files here
	outhetboost.pop();
	outmatboost.pop();
	outpatboost.pop();
	outhet.close();
	outmat.close();
	outpat.close();
	return;
}

void process_child_impute2_haps_trio(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,vector<string>& toprint,int loopnum,int cInds,vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,vector<vector<float>>& fHapArray,std::unordered_map<std::string,int>& mhash,std::unordered_map<std::string,int>& fhash,vector<string>& cids,string chr,bool transmitted){
	vector<string> mhap;
	vector<string> fhap;
	int mind;
	int find;
	//check that loopnum=size of toprint	TODO
	for(int j=0;j<loopnum;++j){
		//print start of gen file line!
		outpat<<chr<<" "<<toprint[j];
		outmat<<chr<<" "<<toprint[j];
		outhet<<chr<<" "<<toprint[j];
		for(int i=0;i<cInds;i++){
			//check child snp is in mother
			if(!(mhash.find(cids[i])==mhash.end()) && !(fhash.find(cids[i])==fhash.end())){
				//find mother index
				mind=mhash[cids[i]];
				find=fhash[cids[i]];
				//find genotypes
				compare_trio_genotypes(cHapArray,mHapArray,fHapArray,j,i,mind,find,outpat,outmat,outhet,transmitted,loopnum);
			}else if(!(mhash.find(cids[i])==mhash.end())){
				//find mother index
				mind=mhash[cids[i]];
				compare_pair_genotypes(cHapArray,mHapArray,j,i,mind,outpat,outmat,outhet,transmitted,loopnum,true);
			}else if(!(fhash.find(cids[i])==fhash.end())){
				//find father index
				find=fhash[cids[i]];
				compare_pair_genotypes(cHapArray,fHapArray,j,i,find,outmat,outpat,outhet,transmitted,loopnum,false);
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}
		outpat<<"\n";
		outmat<<"\n";
		outhet<<"\n";
	}
	return;
}

void compare_trio_genotypes(vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,vector<vector<float>>& fHapArray,int j,int i,int mind,int find,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted,int loopnum){
	if(mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]<0.5){ //AA
		if(fHapArray[j][2*find]<0.5 && fHapArray[j][2*find+1]<0.5){ //AA
			if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if((fHapArray[j][2*find]<0.5 && fHapArray[j][2*find+1]>0.5) || (fHapArray[j][2*find]>0.5 && fHapArray[j][2*find+1]<0.5)){ //Aa
			if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
			}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
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
		}else if(fHapArray[j][2*find]>0.5 && fHapArray[j][2*find+1]>0.5){ //aa
			if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
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
	}else if((mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]>0.5) || (mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]<0.5)){ //Aa
		if(fHapArray[j][2*find]<0.5 && fHapArray[j][2*find+1]<0.5){ //AA
			if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 1 0 0";
				outmat<<" 0 1 0";
				outhet<<" 0 1 0";
			}else{
				outpat<<" 0 0 0";
				outmat<<" 0 0 0";
				outhet<<" 0 0 0";
			}
		}else if((fHapArray[j][2*find]<0.5 && fHapArray[j][2*find+1]>0.5) || (fHapArray[j][2*find]>0.5 && fHapArray[j][2*find+1]<0.5)){ //Aa
			if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
				//So we have a triple het!
				int chap=round(cHapArray[j][2*i]);
				int mhap=round(mHapArray[j][2*mind]);
				int fhap=round(fHapArray[j][2*find]);
				//find haps
				int identical=0;
				std::array<int,3> crefhap,calthap,mrefhap,malthap,frefhap,falthap;
				if(j==0){
					crefhap[0]=round(cHapArray[j][2*i+chap]);
					calthap[0]=round(cHapArray[j][2*i+1-chap]);
					mrefhap[0]=round(mHapArray[j][2*mind+mhap]);
					malthap[0]=round(mHapArray[j][2*mind+1-mhap]);
					frefhap[0]=round(fHapArray[j][2*find+fhap]);
					falthap[0]=round(fHapArray[j][2*find+1-fhap]);	// we know these are all the same, and we need at least one to be different...
					crefhap[1]=round(cHapArray[j+1][2*i+chap]);
					calthap[1]=round(cHapArray[j+1][2*i+1-chap]);
					mrefhap[1]=round(mHapArray[j+1][2*mind+mhap]);
					malthap[1]=round(mHapArray[j+1][2*mind+1-mhap]);
					frefhap[1]=round(fHapArray[j+1][2*find+fhap]);
					falthap[1]=round(fHapArray[j+1][2*find+1-fhap]);
					if(crefhap[1]==mrefhap[1] && calthap[1]==malthap[1] && crefhap[1]==frefhap[1] && calthap[1]==falthap[1]){	// so if the other one is the same we need to search until 
						int k=j+2;												// we find a different one...
						while(identical==0 && k<loopnum){
							crefhap[2]=round(cHapArray[k][2*i+chap]);
							calthap[2]=round(cHapArray[k][2*i+1-chap]);
							mrefhap[2]=round(mHapArray[k][2*mind+mhap]);
							malthap[2]=round(mHapArray[k][2*mind+1-mhap]);
							frefhap[2]=round(fHapArray[k][2*find+fhap]);
							falthap[2]=round(fHapArray[k][2*find+1-fhap]);
							k++;
							identical=(!(crefhap[2]==mrefhap[2] && calthap[2]==malthap[2] && crefhap[2]==frefhap[2] && calthap[2]==falthap[2]));
						}
					}else{							// else we just take the next one for completeness
						crefhap[2]=round(cHapArray[j+2][2*i+chap]);
						calthap[2]=round(cHapArray[j+2][2*i+1-chap]);
						mrefhap[2]=round(mHapArray[j+2][2*mind+mhap]);
						malthap[2]=round(mHapArray[j+2][2*mind+1-mhap]);
						frefhap[2]=round(fHapArray[j+2][2*find+fhap]);
						falthap[2]=round(fHapArray[j+2][2*find+1-fhap]);
					}
				}else if(j==loopnum-1){
					crefhap[2]=round(cHapArray[j][2*i+chap]);
					calthap[2]=round(cHapArray[j][2*i+1-chap]);
					mrefhap[2]=round(mHapArray[j][2*mind+mhap]);
					malthap[2]=round(mHapArray[j][2*mind+1-mhap]);
					frefhap[2]=round(fHapArray[j][2*find+fhap]);
					falthap[2]=round(fHapArray[j][2*find+1-fhap]);
					crefhap[1]=round(cHapArray[j-1][2*i+chap]);
					calthap[1]=round(cHapArray[j-1][2*i+1-chap]);
					mrefhap[1]=round(mHapArray[j-1][2*mind+mhap]);
					malthap[1]=round(mHapArray[j-1][2*mind+1-mhap]);
					frefhap[1]=round(fHapArray[j-1][2*find+fhap]);
					falthap[1]=round(fHapArray[j-1][2*find+1-fhap]);
					if(crefhap[1]==mrefhap[1] && calthap[1]==malthap[1] && crefhap[1]==frefhap[1] && calthap[1]==falthap[1]){
						int k=j-2;
						while(identical==0 && k>=0){
							crefhap[0]=round(cHapArray[k][2*i+chap]);
							calthap[0]=round(cHapArray[k][2*i+1-chap]);
							mrefhap[0]=round(mHapArray[k][2*mind+mhap]);
							malthap[0]=round(mHapArray[k][2*mind+1-mhap]);
							frefhap[0]=round(fHapArray[k][2*find+fhap]);
							falthap[0]=round(fHapArray[k][2*find+1-fhap]);
							k--;
							identical=(!(crefhap[0]==mrefhap[0] && calthap[0]==malthap[0] && crefhap[0]==frefhap[0] && calthap[0]==falthap[0]));
						}
					}else{
						crefhap[0]=round(cHapArray[j-2][2*i+chap]);
						calthap[0]=round(cHapArray[j-2][2*i+1-chap]);
						mrefhap[0]=round(mHapArray[j-2][2*mind+mhap]);
						malthap[0]=round(mHapArray[j-2][2*mind+1-mhap]);
						frefhap[0]=round(fHapArray[j-2][2*find+fhap]);
						falthap[0]=round(fHapArray[j-2][2*find+1-fhap]);
					}
				}else{
					crefhap[0]=round(cHapArray[j-1][2*i+chap]);
					calthap[0]=round(cHapArray[j-1][2*i+1-chap]);
					mrefhap[0]=round(mHapArray[j-1][2*mind+mhap]);
					malthap[0]=round(mHapArray[j-1][2*mind+1-mhap]);
					frefhap[0]=round(fHapArray[j-1][2*find+fhap]);
					falthap[0]=round(fHapArray[j-1][2*find+1-fhap]);
					crefhap[1]=round(cHapArray[j][2*i+chap]);
					calthap[1]=round(cHapArray[j][2*i+1-chap]);
					mrefhap[1]=round(mHapArray[j][2*mind+mhap]);
					malthap[1]=round(mHapArray[j][2*mind+1-mhap]);
					frefhap[1]=round(fHapArray[j][2*find+fhap]);
					falthap[1]=round(fHapArray[j][2*find+1-fhap]);
					if(crefhap[0]==mrefhap[0] && calthap[0]==malthap[0] && crefhap[0]==frefhap[0] && calthap[0]==falthap[0]){
						int k=j+1;
						while(identical==0 && k<loopnum){
							crefhap[2]=round(cHapArray[k][2*i+chap]);
							calthap[2]=round(cHapArray[k][2*i+1-chap]);
							mrefhap[2]=round(mHapArray[k][2*mind+mhap]);
							malthap[2]=round(mHapArray[k][2*mind+1-mhap]);
							frefhap[2]=round(fHapArray[k][2*find+fhap]);
							falthap[2]=round(fHapArray[k][2*find+1-fhap]);
							k++;
							identical=(!(crefhap[2]==mrefhap[2] && calthap[2]==malthap[2] && crefhap[2]==frefhap[2] && calthap[2]==falthap[2]));
						}
						k=j-2;
						while(identical==0 && k>=0){
							crefhap[2]=round(cHapArray[k][2*i+chap]);
							calthap[2]=round(cHapArray[k][2*i+1-chap]);
							mrefhap[2]=round(mHapArray[k][2*mind+mhap]);
							malthap[2]=round(mHapArray[k][2*mind+1-mhap]);
							frefhap[2]=round(fHapArray[k][2*find+fhap]);
							falthap[2]=round(fHapArray[k][2*find+1-fhap]);
							k--;
							identical=(!(crefhap[2]==mrefhap[2] && calthap[2]==malthap[2] && crefhap[2]==frefhap[2] && calthap[2]==falthap[2]));
						}
					}else{
						crefhap[2]=round(cHapArray[j+1][2*i+chap]);
						calthap[2]=round(cHapArray[j+1][2*i+1-chap]);
						mrefhap[2]=round(mHapArray[j+1][2*mind+mhap]);
						malthap[2]=round(mHapArray[j+1][2*mind+1-mhap]);
						frefhap[2]=round(fHapArray[j+1][2*find+fhap]);
						falthap[2]=round(fHapArray[j+1][2*find+1-fhap]);
					}
				}
				if(crefhap==mrefhap && calthap==malthap && crefhap==frefhap && calthap==falthap){
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}else if(crefhap==mrefhap && calthap==falthap){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else if(calthap==malthap && crefhap==frefhap){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}
			}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
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
		}else if(fHapArray[j][2*find]>0.5 && fHapArray[j][2*find+1]>0.5){ //aa
			if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa      //transmitted and inherited here are the same so don't need to test which one we want
				outpat<<" 0 1 0";
				outmat<<" 1 0 0";
				outhet<<" 1 0 0";
			}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
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
	}else if(mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]>0.5){ //aa
		if(fHapArray[j][2*find]<0.5 && fHapArray[j][2*find+1]<0.5){ //AA
			if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
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
		}else if((fHapArray[j][2*find]<0.5 && fHapArray[j][2*find+1]>0.5) || (fHapArray[j][2*find]>0.5 && fHapArray[j][2*find+1]<0.5)){ //Aa
			if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}
			}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
				if(transmitted){
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
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
		}else if(fHapArray[j][2*find]>0.5 && fHapArray[j][2*find+1]>0.5){ //aa
			if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa      //transmitted and inherited here are the same so don't need to test which one we want
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

void compare_pair_genotypes(vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,int j,int i,int mind,filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,bool transmitted,int loopnum,bool mother){
	if(mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]<0.5){ //AA
		if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
			if(mother){
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
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
		}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
			if(transmitted){
				if(mother){
					outpat<<" 1 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}
			}else{
				if(mother){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
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
	}else if((mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]>0.5) || (mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]<0.5)){ //Aa
		if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
			if(transmitted){
				if(mother){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else{
				outpat<<" 1 0 0";
				outmat<<" 1 0 0";
				outhet<<" 0 0 0";
				}
		}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
			//So we have a double het!
			int chap=round(cHapArray[j][2*i]);
			int mhap=round(mHapArray[j][2*mind]);
			//find haps
			int identical=0;
			std::array<int,3> crefhap,calthap,mrefhap,malthap;
			if(j==0){
				crefhap[0]=round(cHapArray[j][2*i+chap]);
				calthap[0]=round(cHapArray[j][2*i+1-chap]);
				mrefhap[0]=round(mHapArray[j][2*mind+mhap]);
				malthap[0]=round(mHapArray[j][2*mind+1-mhap]);
				crefhap[1]=round(cHapArray[j+1][2*i+chap]);
				calthap[1]=round(cHapArray[j+1][2*i+1-chap]);
				mrefhap[1]=round(mHapArray[j+1][2*mind+mhap]);
				malthap[1]=round(mHapArray[j+1][2*mind+1-mhap]);
				if(crefhap[1]==mrefhap[1] && calthap[1]==malthap[1]){
					int k=j+2;
					while(identical==0 && k<loopnum){
						crefhap[2]=round(cHapArray[k][2*i+chap]);
						calthap[2]=round(cHapArray[k][2*i+1-chap]);
						mrefhap[2]=round(mHapArray[k][2*mind+mhap]);
						malthap[2]=round(mHapArray[k][2*mind+1-mhap]);
						k++;
						identical=(!(crefhap[2]==mrefhap[2] && calthap[2]==malthap[2]));
					}
				}else{
					crefhap[2]=round(cHapArray[j+2][2*i+chap]);
					calthap[2]=round(cHapArray[j+2][2*i+1-chap]);
					mrefhap[2]=round(mHapArray[j+2][2*mind+mhap]);
					malthap[2]=round(mHapArray[j+2][2*mind+1-mhap]);
				}
			}else if(j==loopnum-1){
				crefhap[2]=round(cHapArray[j][2*i+chap]);
				calthap[2]=round(cHapArray[j][2*i+1-chap]);
				mrefhap[2]=round(mHapArray[j][2*mind+mhap]);
				malthap[2]=round(mHapArray[j][2*mind+1-mhap]);
				crefhap[1]=round(cHapArray[j-1][2*i+chap]);
				calthap[1]=round(cHapArray[j-1][2*i+1-chap]);
				mrefhap[1]=round(mHapArray[j-1][2*mind+mhap]);
				malthap[1]=round(mHapArray[j-1][2*mind+1-mhap]);
				if(crefhap[1]==mrefhap[1] && calthap[1]==malthap[1]){
					int k=j-2;
					while(identical==0 && k>=0){
						crefhap[0]=round(cHapArray[k][2*i+chap]);
						calthap[0]=round(cHapArray[k][2*i+1-chap]);
						mrefhap[0]=round(mHapArray[k][2*mind+mhap]);
						malthap[0]=round(mHapArray[k][2*mind+1-mhap]);
						k--;
						identical=(!(crefhap[0]==mrefhap[0] && calthap[0]==malthap[0]));
					}
				}else{
					crefhap[0]=round(cHapArray[j-2][2*i+chap]);
					calthap[0]=round(cHapArray[j-2][2*i+1-chap]);
					mrefhap[0]=round(mHapArray[j-2][2*mind+mhap]);
					malthap[0]=round(mHapArray[j-2][2*mind+1-mhap]);
				}
			}else{
				crefhap[0]=round(cHapArray[j-1][2*i+chap]);
				calthap[0]=round(cHapArray[j-1][2*i+1-chap]);
				mrefhap[0]=round(mHapArray[j-1][2*mind+mhap]);
				malthap[0]=round(mHapArray[j-1][2*mind+1-mhap]);
				crefhap[1]=round(cHapArray[j][2*i+chap]);
				calthap[1]=round(cHapArray[j][2*i+1-chap]);
				mrefhap[1]=round(mHapArray[j][2*mind+mhap]);
				malthap[1]=round(mHapArray[j][2*mind+1-mhap]);
				if(crefhap[0]==mrefhap[0] && calthap[0]==malthap[0]){
					int k=j+1;
					while(identical==0 && k<loopnum){
						crefhap[2]=round(cHapArray[k][2*i+chap]);
						calthap[2]=round(cHapArray[k][2*i+1-chap]);
						mrefhap[2]=round(mHapArray[k][2*mind+mhap]);
						malthap[2]=round(mHapArray[k][2*mind+1-mhap]);
						k++;
						identical=(!(crefhap[2]==mrefhap[2] && calthap[2]==malthap[2]));
					}
					k=j-2;
					while(identical==0 && k>=0){
						crefhap[2]=round(cHapArray[k][2*i+chap]);
						calthap[2]=round(cHapArray[k][2*i+1-chap]);
						mrefhap[2]=round(mHapArray[k][2*mind+mhap]);
						malthap[2]=round(mHapArray[k][2*mind+1-mhap]);
						k--;
						identical=(!(crefhap[2]==mrefhap[2] && calthap[2]==malthap[2]));
					}
				}else{
					crefhap[2]=round(cHapArray[j+1][2*i+chap]);
					calthap[2]=round(cHapArray[j+1][2*i+1-chap]);
					mrefhap[2]=round(mHapArray[j+1][2*mind+mhap]);
					malthap[2]=round(mHapArray[j+1][2*mind+1-mhap]);
				}
			}
			if(mother){
				if(crefhap==mrefhap && calthap==malthap){
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}else if(crefhap==mrefhap){
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}else if(calthap==malthap){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}
			}else{
				if(crefhap==mrefhap && calthap==malthap){
					outpat<<" 0 0 0";
					outmat<<" 0 0 0";
					outhet<<" 0 0 0";
				}else if(crefhap==mrefhap){
					if(transmitted){
						outpat<<" 0 0 0";
						outmat<<" 0 1 0";
						outhet<<" 0 0 0";
					}else{
						outpat<<" 1 0 0";
						outmat<<" 0 1 0";
						outhet<<" 0 1 0";
					}
				}else if(calthap==malthap){
					if(transmitted){
						outpat<<" 0 0 0";
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
			}
		}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
			if(transmitted){
				if(mother){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}
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
	}else if(mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]>0.5){ //aa
		if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
			if(transmitted){
				if(mother){
					outpat<<" 0 1 0";
					outmat<<" 0 1 0";
					outhet<<" 0 0 0";
				}else{
					outpat<<" 0 0 0";
					outmat<<" 1 0 0";
					outhet<<" 0 0 0";
				}
			}else{
				if(mother){
					outpat<<" 1 0 0";
					outmat<<" 0 1 0";
					outhet<<" 0 1 0";
				}else{
					outpat<<" 0 1 0";
					outmat<<" 1 0 0";
					outhet<<" 1 0 0";
				}
			}
		}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
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
