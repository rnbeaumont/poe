#include "general_functions.h"
#include "libs/iostreams/src/zlib.cpp"

extern string chr;

int get_num_lines( string lineFile ){
	using namespace std;
	using namespace boost::iostreams;

	int nLines=0;
	ifstream file(lineFile,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<lineFile<<" does not exists"<<endl;
		exit(-20);
	}
        try{
                filtering_istream boost_in;
                // CHECK THAT IT'S GZIPPED AND DON'T GZIP_DECOMPRESSOR IF NOT TODO
		if(lineFile.compare(lineFile.length()-2,2,"gz")==0){
	                boost_in.push(gzip_decompressor());
		}
                boost_in.push(file);
                for(string str; getline(boost_in, str); )
                {
                        /* increment counter as it's read a line */
                        nLines++;
                }
        }
        catch(const gzip_error& e){
                cerr << e.what() << endl;
        }
	file.close();
	return nLines;
}

void check_gen_length(string inFile,unsigned int nCol){
	vector<string> fields;
	ifstream file(inFile,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<inFile<<" does not exists"<<endl;
		exit(-20);
	}
	try{
		filtering_istream boost_in;
		if(inFile.compare(inFile.length()-2,2,"gz")==0){
			boost_in.push(gzip_decompressor());
		}
		boost_in.push(file);
		string str;
		getline(boost_in,str);
		split(fields,str,is_any_of(" \t"));
		if(fields.size()!=(nCol)){
			cerr<<"ERROR: "<<inFile<<" contains differing number of individuals to sample file!"<<endl;
			exit(1);
		}
	}catch(const gzip_error& e){
		cerr<<e.what()<<endl;
	}
	file.close();
	return;
}

void read_sample( string sampleFile, int numlines, unordered_map<string, int> hash, vector<string> cids ){
	using namespace std;
	using namespace boost::iostreams;
	using namespace boost;
	string str;
	vector <string> fields;
	int index=0;
	ifstream file(sampleFile,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<sampleFile<<" does not exists"<<endl;
		exit(-20);
	}
	try{
	        filtering_istream boost_in;
		if(sampleFile.compare(sampleFile.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(file);
		getline(boost_in,str);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			hash[fields[0]]=index;
			cids.push_back(fields[0]);
			index++;
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	file.close();
	return;
}

void read_sample_2(string cSamp,string mSamp,std::unordered_map<std::string,int>& mhash,vector<string>& cids){
	string str1;
	vector <string> fields;
	int index=0;
	ifstream file(cSamp,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<cSamp<<" does not exists"<<endl;
		exit(-20);
	}
	try{
	        filtering_istream boost_in;
		if(cSamp.compare(cSamp.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(file);
		//skip header lines
		getline(boost_in,str1);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			cids.push_back(fields[0]);
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	file.close();

	ifstream mfile(mSamp,ios_base::in | ios_base::binary);
	if(!mfile.is_open()){
		cerr<<"ERROR: file "<<mSamp<<" does not exists"<<endl;
		exit(-20);
	}
	index=0;
	try{
	        filtering_istream boost_in;
		if(mSamp.compare(mSamp.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(mfile);
		//skip header lines
		getline(boost_in,str1);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			mhash[fields[0]]=index;
			index++;
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	mfile.close();
	return;
}

void read_sample_3(string cSamp,string mSamp,string fSamp,std::unordered_map<std::string,int>& mhash,std::unordered_map<std::string,int>& fhash,vector<string>& cids){
	string str1;
	vector <string> fields;
	int index=0;
	ifstream file(cSamp,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<cSamp<<" does not exists"<<endl;
		exit(-20);
	}
	try{
	        filtering_istream boost_in;
		if(cSamp.compare(cSamp.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(file);
		//skip header lines
		getline(boost_in,str1);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			cids.push_back(fields[0]);
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	file.close();

	ifstream mfile(mSamp,ios_base::in | ios_base::binary);
	if(!mfile.is_open()){
		cerr<<"ERROR: file "<<mSamp<<" does not exists"<<endl;
		exit(-20);
	}
	index=0;
	try{
	        filtering_istream boost_in;
		if(mSamp.compare(mSamp.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(mfile);
		//skip header lines
		getline(boost_in,str1);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			mhash[fields[0]]=index;
			index++;
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	mfile.close();

	ifstream ffile(fSamp,ios_base::in | ios_base::binary);
	if(!ffile.is_open()){
		cerr<<"ERROR: file "<<fSamp<<" does not exists"<<endl;
		exit(-20);
	}
	index=0;
	try{
	        filtering_istream boost_in;
		if(fSamp.compare(fSamp.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(ffile);
		//skip header lines
		getline(boost_in,str1);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			fhash[fields[0]]=index;
			index++;
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	ffile.close();
	return;
}

void read_pheno( string pFile, unordered_map<string, string>& phenoMap,string& pHead ){
	using namespace boost;
	string str;
	vector<string> fields;
	ifstream file(pFile);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<pFile<<" does not exists"<<endl;
		exit(-20);
	}
	getline(file, str);
	split(fields,str,is_any_of("\t"),token_compress_on);
	pHead=fields[1];
	phenoMap["HEADER"]=fields[1];
	while(getline(file, str)){
		split(fields,str,is_any_of(" \t"),token_compress_on);
		phenoMap[fields[0]]=fields[1];
	}
	file.close();
	return;
}

void check_dosage_length(string inFile,unsigned int nCol){
	vector<string> fields;
	ifstream file(inFile,ios_base::in | ios_base::binary);
	try{
		filtering_istream boost_in;
		if(inFile.compare(inFile.length()-2,2,"gz")==0){
			boost_in.push(gzip_decompressor());
		}
		boost_in.push(file);
		string str;
		getline(boost_in,str);
		split(fields,str,is_any_of(" \t"));
		if(fields.size()!=(nCol)){
			cerr<<"ERROR: "<<inFile<<" contains differing number of SNPs to info file!\n"<<endl;
			exit(0);
		}
	}catch(const gzip_error& e){
		cerr<<e.what()<<endl;
	}
	file.close();
	return;
}

void read_info( string infoFile, int numlines, unordered_map<string, int> hash, vector<string> cids ){
	using namespace std;
	using namespace boost::iostreams;
	using namespace boost;
	string str;
	vector <string> fields;
	int index=0;
	ifstream file(infoFile,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<infoFile<<" does not exists"<<endl;
		exit(-20);
	}
	try{
	        filtering_istream boost_in;
		if(infoFile.compare(infoFile.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(file);
		getline(boost_in,str);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			hash[fields[0]]=index;
			cids.push_back(fields[0]);
			index++;
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	file.close();
	return;
}

void read_info_2(string cInfo,string mInfo,std::unordered_map<std::string,int>& chash,vector<string>& mids,vector<string>& mstr,string chr){
	string str1;
	vector <string> fields;
	int index=0;
	ifstream file(cInfo,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<cInfo<<" does not exists"<<endl;
		exit(-20);
	}
	try{
	        filtering_istream boost_in;
		if(cInfo.compare(cInfo.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(file);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			chash[fields[0]]=index;
			index++;
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	file.close();

	ifstream mfile(mInfo,ios_base::in | ios_base::binary);
	if(!mfile.is_open()){
		cerr<<"ERROR: file "<<mInfo<<" does not exists"<<endl;
		exit(-20);
	}
	index=0;
	try{
	        filtering_istream boost_in;
		if(mInfo.compare(mInfo.length()-2,2,"gz")==0){
		        boost_in.push(gzip_decompressor());
		}
	        boost_in.push(mfile);
		getline(boost_in,str1);
	        for(string str; getline(boost_in, str); )
	        {
	        	// it's read in a line - split it then store the ids in the unordered_map and the array
	        	split(fields,str,is_any_of(" \t"));
			mids.push_back(fields[0]);
			mstr.push_back(chr+" "+fields[0]+" 0 "+fields[1]+" "+fields[2]);
	        }
	}
	catch(const gzip_error& e){
	        cerr << e.what() <<endl;
	}
	mfile.close();
	return;
}

void print_sample(string sampleFile,string outFile){
        using namespace std;
        using namespace boost::iostreams;
        using namespace boost;
        vector <string> fields;
        ifstream file(sampleFile,ios_base::in | ios_base::binary);
	if(!file.is_open()){
		cerr<<"ERROR: file "<<sampleFile<<" does not exists"<<endl;
		exit(-20);
	}
	ofstream outsamp;
	outsamp.open(outFile+".sample",std::ios_base::out | std::ios_base::binary);
        try{
                filtering_istream boost_in;
                if(sampleFile.compare(sampleFile.length()-2,2,"gz")==0){
                        boost_in.push(gzip_decompressor());
                }
                boost_in.push(file);
                for(string str;getline(boost_in,str);){
	                outsamp<<str<<endl;
		}
        }
        catch(const gzip_error& e){
                cerr << e.what() <<endl;
        }
	outsamp.close();
        file.close();
        return;
}

bgen::bgen(){	//default constructor (needs open command TODO)
	aA="";		//initialise string variables
	aB="";		//initialise string variables
	rsidS="";	//initialise string variables
	currentSnp=0;	//initialise currentSnp to 0 to be able to increment
	readVarHead=false;	//whether the head of thevariant has been read but not the probabilities
}

bgen::bgen(string file){	//constructor to open file
	using namespace std;
	unsigned char buf[4];	//buffer
	aA="";		//initialise string variables
	aB="";		//initialise string variables
	rsidS="";	//initialise string variables
	currentSnp=0;	//initialise currentSnp to 0 to be able to increment
	readVarHead=false;	//whether the head of thevariant has been read but not the probabilities
	input.open(file,ios::binary);	//open file
	if(!input.is_open()){
		cerr<<"ERROR: file "<<file<<" does not exists"<<endl;
		exit(-20);
	}
	input.read((char*)(&offset),4);	//read offset
	input.read((char*)(&LH),4);	//read LH
	input.read((char*)(&M),4);	//read M
	input.read((char*)(&N),4);	//read N
	input.read((char*)(&buf[0]),4);	//read magic number
	////////// TODO CHECK THAT MAGIC NUMBER IS CORRECT!!!!!!!!////////////////////
	// check whether there's anything in the free data area, and if so skip it
	if((LH-20)>0){
		input.seekg((LH-20),ios_base::cur);
	}
	//read the flags
	//Check which format the data is in
	input.read((char*)(&buf[0]),4);
	compressed=((buf[0] >> 0) & 1);	//CompressedSNPBlocks
	layout=((buf[0] >> 2) & 15);	//layout
	sampleIds=((buf[0]>>31) & 1);	//sampleIds		//TODO if this is 1 then I need to read the sample block!!!
	//Header block read so skip to the start of the genotypes
	input.seekg((offset+4),ios_base::beg);
}

void bgen::read_variant_id(){
	using namespace std;
	unsigned int lrsid;
	unsigned char buf[4];
	unsigned long int uncompressedLength;
	//if v1.1 read number of individuals
	if(layout==1){
		input.read((char*)(&nv),4);
	}
	//read variant name
	//vid
	//set all bytes to zero as we'll be reading less than 4 bytes next!
	buf[0]=*(unsigned char*)"\0";
	buf[1]=*(unsigned char*)"\0";
	buf[2]=*(unsigned char*)"\0";
	buf[3]=*(unsigned char*)"\0";
	input.read((char*)(&buf),2);
	lrsid=(buf[1]<<8) | buf[0];
	vid.clear();
	vid.resize(lrsid);
	input.read((char*)(&vid[0]),lrsid);
//	input.read((char*)(&nv),4);//rsid
	input.read((char*)(&buf),2);
	lrsid=(buf[1]<<8) | buf[0];
	rsid.clear();
	rsid.resize(lrsid);
	input.read((char*)(&rsid[0]),lrsid);
	rsidS.clear();
	rsidS=string(rsid.begin(),rsid.end());
	//chr
	input.read((char*)(&buf),2);
	lrsid=(buf[1]<<8) | buf[0];
	chr.clear();
	chr.resize(lrsid);
	input.read((char*)(&chr[0]),lrsid);
	//position
	input.read((char*)(&pos),4);
	//if it's version 1.2 check how many alleles there are
	if(layout==2){
		input.read((char*)(&buf),2);
		K=(buf[1]<<8) | buf[0];
	}
	//allele A
	input.read((char*)(&lrsid),4);
	alleleA.clear();
	alleleA.resize(lrsid);
	input.read((char*)(&alleleA[0]),lrsid);
	aA.clear();
	aA=string(alleleA.begin(),alleleA.end());
	//allele B
	input.read((char*)(&lrsid),4);
	alleleB.clear();
	alleleB.resize(lrsid);
	input.read((char*)(&alleleB[0]),lrsid);
	aB.clear();
	aB=string(alleleB.begin(),alleleB.end());
	//if v1.2 read in the other alleles
	if(layout==2){
		for(unsigned int i=2;i<K;++i){
			input.read((char*)(&lrsid),4);
			alleleA.clear();
			alleleA.resize(lrsid);
			input.read((char*)(&alleleA[0]),lrsid);		//read these in properly and store TODO
		}
	}
	//if v1.2 read in length of rest of data block
	if(layout==2){
		input.read((char*)(&nbytes),4);
		if(compressed){
			input.read((char*)(&uncompressedLength),4);
			nbytes-=4;
		}
	}else if(layout==1){	//v1.1
		//length of compressed data block
		if(compressed){
			input.read((char*)(&nbytes),4);
		}else{
			nbytes=6*N;
		}
	}else{
		cerr<<"ERROR: only bgen version 1.1 and 1.2 supported"<<endl;
		exit(-7);
	}
	//check number of individuals is the same as the number of individuals in the header block, else exit
	if(nv!=N){
		cerr<<"ERROR: Snp ";
		for(vector<unsigned char>::const_iterator i=rsid.begin();i!=rsid.end();++i){
			cerr<<*i;
		}
		cerr<<" has fewer individuals than the bgen header block indicates. I don't know which individuals these are so exiting"<<endl;
		exit(10);
	}
	currentSnp++;
	probs.clear();
	readVarHead=true;
}

void bgen::read_variant_probabilities(){
	using namespace std;
	unsigned long int uncompressedLength;
	unsigned char *uncompressedBuf;
//	unsigned char probBuf[4];
	float probBuf;
	unsigned char *compressedData;
	int uncompressedValue;
//	unsigned char buf[4];
	unsigned int nalleles,minp,maxp,samplePloidy,phased,nbit;

	if(!readVarHead){
		read_variant_id();
	}
	probs.clear();
	if(layout==1){
		uncompressedLength=6*N;
		uncompressedBuf=new unsigned char[uncompressedLength];
		if(compressed){
			compressedData=new unsigned char[nbytes];
			input.read((char*)(&compressedData[0]),nbytes);
			uncompressedValue=uncompress(uncompressedBuf,&uncompressedLength,(const Bytef*)compressedData,nbytes);
			delete [] compressedData;
			if(uncompressedValue!=0){
				cout<<"Error reading variant, exiting!"<<endl;
				exit(10);
			}
		}else{
			input.read((char*)(&uncompressedBuf[0]),uncompressedLength);
		}
		for(unsigned int i=0;i<uncompressedLength;i+=2){
			probBuf=(uncompressedBuf[i+1]<<8) | uncompressedBuf[i];
			probs.push_back(probBuf/32768.0);
		}
		delete [] uncompressedBuf;
	}else if(layout==2){
		uncompressedBuf=new unsigned char[uncompressedLength];
		if(compressed){
			compressedData=new unsigned char[nbytes];
			input.read((char*)(&compressedData[0]),nbytes);
			uncompressedValue=uncompress(uncompressedBuf,&uncompressedLength,(const Bytef*)compressedData,nbytes);
			delete [] compressedData;
			if(uncompressedValue!=0){
				cout<<"Error reading variant, exiting!"<<endl;
				exit(10);
			}
		}else{
			input.read((char*)(&uncompressedBuf[0]),uncompressedLength);
		}
		//number of individuals in data block (must equal N)
		nv=(uncompressedBuf[3]<<24) | (uncompressedBuf[2]<<16) | (uncompressedBuf[1]<<8) | uncompressedBuf[0];	//TODO check this equals header (N)
		//number of alleles (must equal 
		nalleles=(uncompressedBuf[5]<<8) | uncompressedBuf[4];
		//min ploidy
		minp=uncompressedBuf[6];
		//max ploidy
		maxp=uncompressedBuf[7];
		//is the data phased
		phased=uncompressedBuf[8+nv];
		//how many bits are used to store each probability
		nbit=uncompressedBuf[9+nv];
		for(unsigned int i=0;i<nv;++i){
			//read ploidy of sample
			samplePloidy=uncompressedBuf[8+i];	//this isn't right!!!
		}

		delete [] uncompressedBuf;
	}
	readVarHead=false;
}

void bgen::skip_variant_probabilities(){
	using namespace std;
	if(!readVarHead){
		read_variant_id();
	}
	input.seekg(nbytes,ios_base::cur);
	readVarHead=false;
}

gen::gen(){
	
}

gen::gen(string file){
	using namespace std;
	using namespace boost::iostreams;
	vector<string> fields;
//	nSnps=0;
	filename=file;	//set filename
	ifstream fileOpen(file,ios_base::in | ios_base::binary);
	if(!fileOpen.is_open()){
		cerr<<"ERROR: file "<<file<<" does not exists"<<endl;
		exit(-20);
	}
	filtering_istream boost_in;
	try{
		if(file.compare(filename.length()-2,2,"gz")==0){
			boost_in.push(gzip_decompressor());
		}
		boost_in.push(fileOpen);
		string temp;
		std::getline(boost_in,temp);
		nSnps++;
		split(fields,temp,is_any_of(" \t"));
		nSampLine1=(fields.size()-5)/3;
//		for(string str;std::getline(boost_in,str);){
//			nSnps++;
//		}
	}catch(const gzip_error& e){
		cerr << e.what() << endl;
	}
	fileOpen.close();
	//read the number of SNPs, so now open the file for general reading
	try{
		iFile.open(file,ios_base::in | ios_base::binary);
		if(file.compare(filename.length()-2,2,"gz")==0){
			openFile.push(gzip_decompressor());
		}
		openFile.push(iFile);
	}catch(const gzip_error& e){
		cerr << e.what() << endl;
	}
}

void gen::read_num_lines(){
	ifstream fileOpen(filename,ios_base::in | ios_base::binary);
	filtering_istream boost_in;
	try{
		if(filename.compare(filename.length()-2,2,"gz")==0){
			boost_in.push(gzip_decompressor());
		}
		boost_in.push(fileOpen);
		string temp;
		std::getline(boost_in,temp);
		for(string str;std::getline(boost_in,str);){
			nSnps++;
		}
	}catch(const gzip_error& e){
		cerr << e.what() << endl;
	}
	fileOpen.close();
}

void gen::set_num_samples(unsigned int nSamp){
	nSample=nSamp;
	//we've set the number of samples, now check that it's the same as the number of samples in the gen file
	//or we'll get a segfault which will confuse the user
//	if(nSample!=nSampLine1){
//		cerr<<"ERROR: Number of individuals in sample file does not equal number of individuals in gen file!\n\n"<<endl;
//		exit(-2);
//	}
}

bool gen::getline(){
	bool read;
	if(std::getline(openFile,line)){
		boost::split(fields,line,boost::is_any_of(" \t"));
		read=1;
	}else{
		read=0;
	}
	return read;
}

void gen::split_fields(){
	for(unsigned int i=0;i<3*nSample;i++){
		if(fields[i+5].compare("NA") == 0){
			probs[i]=0;
		}else{
			float buf=stof(fields[i+5]);
			probs[i]=buf;
		}
	}
}

void gen::split_haps(unsigned int z){
	for(unsigned int i=0;i<2*nSample;i++){
		if(fields[i+5].compare("NA") == 0){
			haps[z][i]=0;
		}else{
			float buf=stof(fields[i+5]);
			haps[z][i]=buf;
		}
	}
}

string gen::outString(){
	return fields[1]+" "+fields[2]+" "+fields[3]+" "+fields[4];
}
