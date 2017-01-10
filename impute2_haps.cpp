#include "general_functions.h"
#include "impute2_haps.h"

void process_impute2_haps(string childFile,string motherFile,string cSamp,string mSamp,string outFile,string chr,bool transmitted,unsigned int nsnp,bool sample,bool interactions){
	string str;
	vector<string> mfields, cfields;
	unsigned int linec=0,linem=0;
	int nel=0;
	unsigned int mSnp,cInds,mInds;
	bool last=0;
	gen cGen(childFile);
	gen mGen(motherFile);
	vector<string> titles;
	unordered_map<string,int> mhash;
	vector<string> cids;
	/* get number of SNPs in both mother and child files - need to decrement due to header line in info files. This is still number of SNPs not last index!!! */
	cInds=get_num_lines(cSamp);
	cInds-=2;
	cGen.set_num_samples(cInds);
	mInds=get_num_lines(mSamp);
	mInds-=2;
	mGen.set_num_samples(mInds);
	cGen.haps.resize(nsnp,vector<float>(cInds*2,0));
	mGen.haps.resize(nsnp,vector<float>(mInds*2,0));
	cout<<"Processing Sample files"<<endl;
	/* read sample files into hashes and arrays (I think I need both for both but need to check) */
	read_sample_2(cSamp,mSamp,mhash,cids);
        if(sample){
                print_sample(cSamp,outFile);
        }
	//loop through files:
	//open output files
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
		for(unsigned int z=0;z<nsnp;++z){
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
			first:
			if(!mGen.fields[1].compare(cGen.fields[1]) && !mGen.fields[2].compare(cGen.fields[2])){
				nel++;
				mGen.split_haps(z);
				cGen.split_haps(z);
				titles.push_back(cGen.fields[1]+" "+cGen.fields[2]+" "+cGen.fields[3]+" "+cGen.fields[4]);
			}else if(stoi(mGen.fields[2])>stoi(cGen.fields[2])){
				if(!cGen.getline()){
					last=1;
					goto endc;
				}
				linec++;
				goto first;
			}else if(stoi(mGen.fields[2])<stoi(cGen.fields[2])){
				if(!mGen.getline()){
					last=1;
					goto endc;
				}
				linem++;
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
		if(interactions){
			process_child_impute2_haps_interaction(outhetboost,titles,nel,cInds,cGen.haps,mGen.haps,mhash,cids,chr);	//nel will tell process_child how many elements are in the array, should be the same as titles.size();
		}else{
			process_child_impute2_haps(outpatboost,outmatboost,outhetboost,titles,nel,cInds,cGen.haps,mGen.haps,mhash,cids,chr);	//nel will tell process_child how many elements are in the array, should be the same as titles.size();
		}
		if(!last){
			goto loopStart;
		}
		//process_child_impute2_haps(outpatboost,outmatboost,outhetboost,titles,nel,cInds,cHapArray,mHapArray,mhash,cids,chr);	//nel will tell process_child how many elements are in the array, should be the same as titles.size();
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
	cout<<"done"<<endl;
	return;
}

void process_child_impute2_haps(filtering_ostream& outpat,filtering_ostream& outmat,filtering_ostream& outhet,vector<string>& toprint,int loopnum,int cInds,vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr){
	string str;
	vector<string> phap;
	vector<string> mhap;
	vector<string> hhap;
	vector<string> fields;
	ofstream pfile;
	int mind;

	//check that loopnum=size of toprint	TODO
	for(int j=0;j<loopnum;++j){		// Snps
		//print start of gen file line!
		outpat<<chr<<" "<<toprint[j];
		outmat<<chr<<" "<<toprint[j];
		outhet<<chr<<" "<<toprint[j];
		for(int i=0;i<cInds;i++){	// Individuals
			//check child snp is in mother
			if(!(mhash.find(cids[i])==mhash.end())){
				//find mother index
				mind=mhash[cids[i]];
				//find genotypes
				if(mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]<0.5){ //AA
					if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
						outpat<<" 1 0 0";
						outmat<<" 1 0 0";
						outhet<<" 0 0 0";
					}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
						outpat<<" 0 1 0";
						outmat<<" 1 0 0";
						outhet<<" 1 0 0";
					}else{
						outpat<<" 0 0 0";
						outmat<<" 0 0 0";
						outhet<<" 0 0 0";
					}
				}else if((mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]>0.5) || (mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]<0.5)){ //Aa
					if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
						outpat<<" 1 0 0";
						outmat<<" 1 0 0";
						outhet<<" 0 0 0";
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
					}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
						outpat<<" 0 1 0";
						outmat<<" 0 1 0";
						outhet<<" 0 0 0";
					}else{
						outpat<<" 0 0 0";
						outmat<<" 0 0 0";
						outhet<<" 0 0 0";
					}
				}else if(mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]>0.5){ //aa
					if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
						outpat<<" 1 0 0";
						outmat<<" 0 1 0";
						outhet<<" 0 1 0";
					}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
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
	}
	return;
}


void process_child_impute2_haps_interaction(filtering_ostream& outhet,vector<string>& toprint,int loopnum,int cInds,vector<vector<float>>& cHapArray,vector<vector<float>>& mHapArray,std::unordered_map<std::string,int>& mhash,vector<string>& cids,string chr){
	string str;
	vector<string> hhap;
	vector<string> fields;
	ofstream pfile;
	int mind;

	//check that loopnum=size of toprint	TODO
	for(int j=0;j<loopnum;++j){		// Snps
		//print start of gen file line!
		outhet<<chr<<" "<<toprint[j];
		for(int i=0;i<cInds;i++){	// Individuals
			//check child snp is in mother
			if(!(mhash.find(cids[i])==mhash.end())){
				//find mother index
				mind=mhash[cids[i]];
				//find genotypes
				if(mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]<0.5){ //AA
					if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
						outhet<<" 1 0 0";
					}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
						outhet<<" 0 1 0";
					}else{
						outhet<<" 0 0 0";
					}
				}else if((mHapArray[j][2*mind]<0.5 && mHapArray[j][2*mind+1]>0.5) || (mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]<0.5)){ //Aa
					if(cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]<0.5){ // AA
						outhet<<" 0 1 0";
					}else if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
						outhet<<" 1 0 0";
					}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
						outhet<<" 0 1 0";
					}else{
						outhet<<" 0 0 0";
					}
				}else if(mHapArray[j][2*mind]>0.5 && mHapArray[j][2*mind+1]>0.5){ //aa
					if((cHapArray[j][2*i]<0.5 && cHapArray[j][2*i+1]>0.5) || (cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]<0.5)){ // Aa
						outhet<<" 0 1 0";
					}else if(cHapArray[j][2*i]>0.5 && cHapArray[j][2*i+1]>0.5){ // aa
						outhet<<" 1 0 0";
					}else{
						outhet<<" 0 0 0";
					}
				}else{
					outhet<<" 0 0 0";
				}
			}else{
				outhet<<" 0 0 0";
			}
		}
		outhet<<"\n";
	}
	return;
}
