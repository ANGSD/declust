//gpl thorfinn@binf.ku.dk
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <ctime>
#include <string.h>
#include <math.h>

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

size_t *histogram = NULL;
size_t histogram_l = 4096;


int globalX(double x, int xlen, short int swath){
	return ((swath-1)*(xlen))+(x);
}
//int globalY(double y, int ylen, short int nTiles, short int tile){
int globalY(double y, int ylen, unsigned short int nTiles, unsigned short int tile){
	return ((nTiles-tile)*(ylen))+(y); 
}

//int getTileInfo(char *tileId, short int *surf, short int *swath, short int *tile) {
int getTileInfo(char *tileId, unsigned short int *surf, unsigned short int *swath, unsigned short int *tile) {
	std::string str(tileId);
	*surf = std::stoi(str.substr(0,1));
	*swath = std::stoi(str.substr(1,1));
	*tile = std::stoi(str.substr(2,2));
}

void tsktsk(){
	fprintf(stderr,"\t-> Incrementing histogram of duplicates from:%lu to  %lu\n",histogram_l,2*histogram_l);
	size_t *tmptmp = new size_t[2*histogram_l];
	for(int i=0;i<2*histogram_l;i++)
		tmptmp[i] = 0;
	for(int i=0;i<histogram_l;i++)
		tmptmp[i] = histogram[i];
	histogram_l = 2*histogram_l;
	histogram=tmptmp;
}


typedef struct{
	bam1_t **d;//data
	unsigned l;//at pos
	unsigned m;//maxpos;
	//int lid;//libid
}queue_t;

double pxdist=12000;
size_t totaldups=0;
size_t pcrdups=0;
size_t clustdups=0;
size_t nproc=0;
size_t purecount;

int nreads_per_pos=4;//assuming this is the number reads per pos. If larger the program will reallocate

char out_mode[5]="wb";

//below is just the average/mean calculated fancy
double CMA =0; //cumulative moving average

queue_t *init_queue_t(int l){
	//  fprintf(stderr,"initializing queue with: %d elements\n",l);
	queue_t *ret =(queue_t *) malloc(sizeof(queue_t));
	ret->d =(bam1_t **) malloc(l*sizeof(bam1_t*));
	ret->l=0;
	ret->m=l;
	//ret->lid=1;
	for(int i=0;i<ret->m;i++){
		//    fprintf(stderr,"queue[%d] init\n",i);
		ret->d[i] = bam_init1();
	}
	return ret;
}



void realloc_queue(queue_t *q){
	//  fprintf(stderr,"reallcing from: q:%d to q:%d\n",q->m,2*q->m);

	for(int i=0;0&&i<q->l;i++)
		fprintf(stderr,"inqueu[%d].pos:%d\n",i,(int)q->d[i]->core.pos);

	bam1_t **d2 = (bam1_t **) malloc(2*q->m*sizeof(bam1_t*));

	for(int i=0;i<q->l;i++)
		d2[i] = q->d[i];
	for(int i=q->l;i<2*q->m;i++){
		d2[i] = bam_init1();
		d2[i]->core.pos=-1;
	}

	free(q->d);
	q->d=d2;
	q->m=2*q->m;

	for(int i=0;0&&i<q->m;i++)
		fprintf(stderr,"onqueu[%d].pos:%d\n",i,(int)q->d[i]->core.pos);

}

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

// FIXME: we should also check the LB tag associated with each alignment
//unconstanted
char *bam_get_library(bam_hdr_t *h, const bam1_t *b)
{
	// Slow and inefficient.  Rewrite once we get a proper header API.
	const char *rg;
	char *cp = h->text;
	rg = (char *)bam_aux_get(b, "RG");

	if (!rg)
		return NULL;
	else
		rg++;

	// Header is guaranteed to be nul terminated, so this is valid.
	while (*cp) {
		char *ID, *LB;
		char last = '\t';

		// Find a @RG line
		if (strncmp(cp, "@RG", 3) != 0) {
			while (*cp && *cp != '\n') cp++; // skip line
			if (*cp) cp++;
			continue;
		}

		// Find ID: and LB: keys
		cp += 4;
		ID = LB = NULL;
		while (*cp && *cp != '\n') {
			if (last == '\t') {
				if (strncmp(cp, "LB:", 3) == 0)
					LB = cp+3;
				else if (strncmp(cp, "ID:", 3) == 0)
					ID = cp+3;
			}
			last = *cp++;
		}

		if (!ID || !LB)
			continue;

		// Check it's the correct ID
		if (strncmp(rg, ID, strlen(rg)) != 0 || ID[strlen(rg)] != '\t')
			continue;

		// Valid until next qual
		static char LB_text[1024];
		for (cp = LB; *cp && *cp != '\t' && *cp != '\n'; cp++)
			;
		strncpy(LB_text, LB, MIN(cp-LB, 1023));
		LB_text[MIN(cp-LB, 1023)] = 0;

		// Return it; valid until the next qual.
		return LB_text;
	}

	return NULL;
}


struct ltstr
{
	bool operator()(const char* s1, const char* s2) const
	{
		return strcmp(s1, s2) < 0;
	}
};

typedef std::map<char*,int,ltstr> aMap;
//doublemap

aMap char2int;

typedef struct{
	bam1_t *d;
	//unsigned long int xs;
	//unsigned long int ys;
	long long int xs;
	long long int ys;
	int seqlen;
}reldata;//<-releavant data

//typedef std::map<size_t,std::map<size_t,std::vector<reldata>> > dMap;


char *mystr =new char[2048];

double euc_dist(reldata &a,reldata &b){

	double dx=a.xs-b.xs;
	double dy=a.ys-b.ys;
	double dist = sqrt(dx*dx+dy*dy);

	return dist;
}

void print_clusters(std::vector<std::vector<int> > &clusters){
	fprintf(stderr,"[print clusters] cluster.size():%lu\n",clusters.size());
	for(int i=0;i<clusters.size();i++){
		fprintf(stderr,"[print clusters] cluster:%i \n",i);
		std::vector<int> &tmp = clusters[i];
		for(int j=0;j<tmp.size();j++)
			fprintf(stderr," %d ",tmp[j]);
		fprintf(stderr,"\n");
	}
}

//void plugin(std::map<size_t,std::vector<reldata> > &mymap,bam1_t *b,bam_hdr_t *hdr){
//void plugin(std::map<size_t,std::map<size_t,std::vector<reldata> >> &mymap, std::map<size_t,std::vector<reldata> > &inmap, bam1_t *b,bam_hdr_t *hdr){
void plugin(std::map<size_t,std::map<size_t,std::vector<reldata> >> &mymap, bam1_t *b,bam_hdr_t *hdr){
	reldata point;
	point.d=b;

	mystr = strncpy(mystr,bam_get_qname(b),2048);
	//    fprintf(stderr,"mystr: \'%s\'\n",mystr);
	strtok(mystr,"\n\t:");//machine
	strtok(NULL,"\n\t:");//runname
	strtok(NULL,"\n\t:");//flowcell
	unsigned short int lane = atoi(strtok(NULL,"\n\t:"));
	//int tile = atoi(strtok(NULL,"\n\t:"));

	unsigned short int surf;
	unsigned short int swath;
	unsigned short int tile;
	getTileInfo(strtok(NULL,"\n\t:"), &surf, &swath, &tile);

	// novaseq specific values, inferred from data
	int xlen=32103;
	int ylen=36059;
	unsigned short int nTiles=78;

	//printf("\n\n%d %d %d\n\n",surf,swath,tile);
	//point.xs = atoi(strtok(NULL,"\n\t:"));
	//point.ys = atoi(strtok(NULL,"\n\t:"));
	point.xs = globalX(atoi(strtok(NULL,"\n\t:")),xlen,swath);
	point.ys = globalY(atoi(strtok(NULL,"\n\t:")),ylen,nTiles,tile);
	//globalize(&point.xs, &point.ys, swath, tile);
	point.seqlen = b->core.l_qseq;

	int libid = 0;
	char *lb = bam_get_library(hdr,b);
	if(lb){
		if(char2int.find(lb)==char2int.end())
			char2int[strdup(lb)] = char2int.size();
		libid = char2int.find(lb)->second;
		if(char2int.size()>998){
			fprintf(stderr,"number of libraries is almost above 998, program will exit. Program should be updated\n");
			exit(0);

		}
	}
	//fprintf(stderr,"surface:%d, swath:%d, tile:%d, libid:%d lane:%d rlen:%d xs:%d ys:%d\n",surf,swath,tile,libid,lane,b->core.l_qseq,point.xs,point.ys);

	// 
	// given tile id 1234
	// surf, swath and tile are officially defined as following:
	//
	// 1			2		34
	// -			-		--
	// surface		swath	tile
	// 

	// 
	// READLENGTH+LIBID+LANE+SURFACE
	// XXX+YYY+Z+T 

	//library id; max 999
	size_t key=libid;
	//read length; max 999
	key += b->core.l_qseq*1e3;

	size_t key2=surf;
	key2 += lane*1e1;

	//mymap is the outer map
	//fprintf(stderr,"\nkey:%d, key2:%d surface:%d, swath:%d, tile:%d, libid:%d lane:%d rlen:%d xs:%d ys:%d\n----\n\n",key,key2,surf,swath,tile,libid,lane,b->core.l_qseq,point.xs,point.ys);
	std::map<size_t,std::map<size_t, std::vector<reldata>>>::iterator it= mymap.find(key);
	//key not found
	if(it==mymap.end()){
		std::map<size_t,std::vector<reldata> > inmap;
		//printf("\nNEW KEY1\n");
		std::vector<reldata> rd;
		rd.push_back(point);
		mymap[key][key2]=rd;
	}else{
		std::map<size_t,std::vector<reldata> >::iterator in =it->second.find(key2);

		//key2 not found
		if(in==it->second.end()){
			//printf("\nNEW KEY2\n");
			std::vector<reldata> rd;
			rd.push_back(point);
			it->second[key2]=rd;
		}else{
			//printf("\nKEY1+KEY2 exists, add to list\n");
			in->second.push_back(point);
		}
	}

}

//fp=nocluster,fp2=onlyclusterdup,fp3=pure
/*
   pure should contain no duplicates (but one representative)
   noclust should contain pure+one read from each cluster and all nonclusterdups
   clustdup should contain all clusterdups
   */

//OUTPUT FILES
//
//fp =  noclusterdup
//fp2 = onlyclusterdup

void plugout(std::map<size_t,std::map<size_t,std::vector<reldata> >> &mymap, bam_hdr_t *hdr, samFile *fp, samFile *fp2,std::vector<size_t> &counter){

	//outer iteration:
	//|__ same library id+same read length
	// same mapping position is already a requirement before plugout is called

	//printf("\n\n----> OUTMAPSIZE: %d\n\n",mymap.size());
	size_t dcount;
	int plug=0;
	for(std::map<size_t,std::map<size_t,std::vector<reldata> >>::iterator it=mymap.begin();it!=mymap.end();it++) {
		//dcount= duplicate fragment count, to give preseq
		dcount=0;
		//printf("\n\n----> OUTLOOP: %d\n\n",it->first);
		//
		std::map<size_t,std::vector<reldata>>::iterator in;
		//printf("\n\n----> INMAPSIZE: %d\n\n",it->second.size());


		//inner iteration:
		//|__ same surface+same lane
		for (in = it->second.begin();in !=it->second.end();++in){

			//printf("\n------>INLOOP: %d\n",in->first);
			//printf("\n\nFirst key: %d\n",it->first);
			//printf("\n\nSecond key: %d\n",in->first);

			std::vector<reldata> &rd=in->second;

			//std::vector<reldata> &rd=it->second;
			//fprintf(stderr,"\nrd size: %lu\n",rd.size());
			totaldups+=rd.size();

			//just a pcr duplicate alone in SURFACE+TILE
			if(rd.size()==1){
				//pcrdups++;
				dcount++;

				//noclusterdupcount++;
				if(fp)
					assert(sam_write1(fp, hdr,rd[0].d)>=0);
				continue;
			}


			//#if 0
			for(int i=0;i<rd.size();i++)
				//fprintf(stderr,"\tcc key: out %lu in %lu/%lu val: xs:%d ys:%d pos:%d\n",it->first,in->first,rd.size(),rd[i].xs,rd[i].ys,rd[i].d->core.pos+1);
			//#endif


			if(rd.size()==2){
				double dist = euc_dist(rd[0],rd[1]);
				//fprintf(stderr,"dist is:%f\n",dist);

				if(dist>pxdist){
					//not part of same cluster
					//we have 2 pcr duplicates
					dcount+=2;
					//pcrdups+=2 ;
					//noclusterdupcount+=2;
					if(fp){
						assert(sam_write1(fp, hdr,rd[0].d)>=0);      
						assert(sam_write1(fp, hdr,rd[1].d)>=0);
					}
				}else{
					//same cluster
					//two reads form a cluster
					//we have 1 pcr duplicate and 1 cluster duplicate
					dcount++;
					//pcrdups++;
					//noclusterdupcount++;

					clustdups++ ;
					if(fp)
						assert(sam_write1(fp, hdr,rd[0].d)>=0);
					if(fp2){
						//assert(sam_write1(fp2, hdr,rd[0].d)>=0);
						assert(sam_write1(fp2, hdr,rd[1].d)>=0);
					}
				}
				continue;
				//one of them is the original read
				//-1 to not to report the original one
				//pcrdups--;
			}
			if(rd.size()==3){
				double dist[3] = {euc_dist(rd[0],rd[1]),euc_dist(rd[0],rd[2]),euc_dist(rd[1],rd[2])};
				double d01=dist[0];
				double d02=dist[1];
				double d12=dist[2];
				int val=0; //nr of reads within pxdist
				for(int i=0;i<3;i++)
					if(dist[i]<pxdist)
						val++;

				if(val==0){
					// not part of the same cluster
					// 3 pcr duplicates
					dcount += 3;
					//pcrdups +=3 ;
					//noclusterdupcount+=3;
					if(fp){
						assert(sam_write1(fp, hdr,rd[0].d)>=0);      
						assert(sam_write1(fp, hdr,rd[1].d)>=0);
						assert(sam_write1(fp, hdr,rd[2].d)>=0);
					}
				}else if (val>=2){
					// 3 reads form a cluster
					// one pcr duplicate + 2 cluster duplicates
					dcount++;
					//pcrdups++;
					//noclusterdupcount++;
					clustdups+=2 ;
					if(fp)
						assert(sam_write1(fp, hdr,rd[0].d)>=0);
					if(fp2){
						//assert(sam_write1(fp2, hdr,rd[0].d)>=0);
						assert(sam_write1(fp2, hdr,rd[1].d)>=0);
						assert(sam_write1(fp2, hdr,rd[2].d)>=0);
					}
				}else if (val==1){
					//2 in cluster one outside
					//1 pcr duplicate + (1 pcr duplicate+1cluster duplicate)
					//pcrdups +=2;
					dcount +=2;
					//noclusterdupcount +=2 ;
					clustdups++ ;

					if(d01<pxdist){
						//rd0 and rd1 defines a cluster, rd2 outside
						clustdups++ ;
						if(fp){
							assert(sam_write1(fp, hdr,rd[0].d)>=0);
							assert(sam_write1(fp, hdr,rd[2].d)>=0);
						}
						if(fp2){
							//assert(sam_write1(fp2, hdr,rd[0].d)>=0);
							assert(sam_write1(fp2, hdr,rd[1].d)>=0);
						}
					}
					else if(d02<pxdist){
						//rd0 and rd2 defines a cluster, rd1 outside
						if(fp){
							assert(sam_write1(fp, hdr,rd[0].d)>=0);
							assert(sam_write1(fp, hdr,rd[1].d)>=0);
						}
						if(fp2){
							//assert(sam_write1(fp2, hdr,rd[0].d)>=0);
							assert(sam_write1(fp2, hdr,rd[2].d)>=0);
						}
					}
					else if(d12<pxdist){
						//rd1 and rd2 defines a cluster, rd0 outside
						if(fp){
							assert(sam_write1(fp, hdr,rd[0].d)>=0);
							assert(sam_write1(fp, hdr,rd[1].d)>=0);
						}
						if(fp2){
							//assert(sam_write1(fp2, hdr,rd[1].d)>=0);
							assert(sam_write1(fp2, hdr,rd[2].d)>=0);
						}
					}else{
						fprintf(stderr,"never happens\n");
						exit(0);
					}

				}else{
					fprintf(stderr,"never happens");
					exit(0);
				}
				//one of them is the original read
				//-1 to not to report the original one
				//pcrdups--;
				continue;
			}
			//printf("\n\nRD SIZE > 3\n\n");

			//vector of vectors, containing ids for the reads that cluster together
			std::vector<std::vector<int> > clusters;

			for(int i=0;i<rd.size();i++) {
				//print_clusters(clusters);
				//      fprintf(stderr,"analysing rd:%d\n",i);
				char dingdongsong[rd.size()];//initialize a hit vector that tells us if the current read is close enough to the different clusters
				memset(dingdongsong,0,rd.size());
				for(int j=0;j<clusters.size();j++){//loop over the different clutsters
					std::vector<int> aclust = clusters[j];
					//	fprintf(stderr,"aclust.size():%lu\n",aclust.size());
					for(int jj=0;jj<aclust.size();jj++){//loop over every read for each cluster
						double dist = euc_dist(rd[i],rd[aclust[jj]]);
						//  fprintf(stderr,"\t-> dist(%d,%d):%f\n",i,aclust[jj],dist);
						if(dist<pxdist){
							//	    fprintf(stderr,"hit at dongdong:%d\n",jj);
							dingdongsong[j]=1;
							continue;
						}
					}
				}
#if 0
				for(int i=0;i<rd.size();i++)
					fprintf(stderr,"dingdong i:%d %d\n",i,dingdongsong[i]);
#endif 
				//now dingdongsong contains a 0/1 array indicating which existing clusters it belongs to.
				int nclust=0;//counter for how many clusters
				for(int s=0;s<rd.size();s++){
					//	fprintf(stderr,"analysinz look up table:%d/%lu=%d\n",s,rd.size(),dingdongsong[s]);
					nclust += dingdongsong[s];
				}
				//      fprintf(stderr,"nclust: %d\n",nclust);
				//
				//
				if(nclust==0){//case where it is not within pixel dist to any
					//fprintf(stderr,"\t-> creating new cluster\n");
					std::vector<int> tmp;tmp.push_back(i);
					clusters.push_back(tmp);
					continue;
				}if(nclust==1){//
					//	fprintf(stderr,"only close enough to one cluster put it back in that clusterlist\n");
					for(int j=0;j<rd.size();j++)
						if(dingdongsong[j]==1){
							//	    fprintf(stderr,"pushing back i:%d at j:%d\n",i,j);
							clusters[j].push_back(i);
							continue;
						}
				}if(nclust>1){
					//	fprintf(stderr,"multiclust: nclust:%d\n",nclust);
					//	keep array is number of clusters long, and will contain which clusters to merge
					int keep[nclust];
					int at=0;
					for(int j=0;j<rd.size();j++){
						if(dingdongsong[j]){
							//  fprintf(stderr,"keep[%d]:%d\n",at,j);
							keep[at++] = j;
						}
					}
					for(int j=1;j<nclust;j++){
						clusters[keep[0]].insert(clusters[keep[0]].end(),clusters[j].begin(),clusters[j].end());
					}

					//print_clusters(clusters);
					for(int j=nclust-1;j>0;j--){
						clusters.erase(clusters.begin()+keep[j]);
					}
					// we started with read i, and this caused us to merge clusters because it is within pixel distance to both
					// , now we finally append the ith read into the merged cluster list
					clusters[keep[0]].push_back(i);
					//	print_clusters(clusters);
				}
			}

			//printf("\n\n---------printing clusters:\n\n");
			//print_clusters(clusters);
			//fprintf(stderr,"\t-> Flushing\n");
			//loop over groupings

			for(int i=0;i<clusters.size();i++){

				//fprintf(stderr,"bangbang: %lu tid:%d pos: %d\n",clusters.size(),rd[0].d->core.tid,rd[0].d->core.pos+1);
				std::vector<int> &tmp = clusters[i];
				//fprintf(stderr,"tmp.size():%lu\n",tmp.size());


				//TODO can it be empty tho?
				if(tmp.size()>0){

					//one cluster of cluster duplicates
					//one pcr duplicate [founder]
					//rest is its cluster duplicates

					//pcrdups++;
					//noclusterdupcount++;
					dcount++;

					if(fp)
						assert(sam_write1(fp, hdr,rd[tmp[0]].d)>=0);

					//skip the first read, assume it is the founder pcr duplicate
					//rest of the reads are cluster duplicates
					for(int j=1;j<tmp.size();j++){
						clustdups++;
						if(fp2)
							assert(sam_write1(fp2, hdr,rd[tmp[j]].d)>=0);
					}
				}
			}
			//one of them is the original read
			//-1 to not to report the original one
			//pcrdups--;
		}
		//fprintf(stderr,"\nhere COUNTER:%d\n",dcount);
		counter.push_back(dcount);
		dcount=0;
		//pcrdups--;
		pcrdups=totaldups-clustdups;
		//fprintf(stderr,"\nhere pcrdups:%d\n",pcrdups);
	}
}

void printmap(FILE *fp,std::map<size_t,std::vector<reldata> > &mymap){
#if 1
	fprintf(fp,"std::map.size:%lu\n",mymap.size());
	for(std::map<size_t,std::vector<reldata> >::iterator it=mymap.begin();it!=mymap.end();it++){
		fprintf(fp,"key:%lu\n",it->first);
		std::vector<reldata> &rd=it->second;
		for(int i=0;i<rd.size();i++)
			fprintf(fp,"\tval: xs:%d ys:%d\n",rd[i].xs,rd[i].ys);
	}
#endif
}

void do_magic(queue_t *q,bam_hdr_t *hdr,samFile *fp,samFile *fp2,samFile *nodupFP){

	//fprintf(stderr,"do_magic queue->l:%d queue->m:%d chr:%d pos:%d\n",q->l,q->m,q->d[0]->core.tid,q->d[0]->core.pos);

	//fprintf(stderr,"@@@@@@info\t%d\t%d\n",q->d[0]->core.pos+1,q->l);
	//totaldups += q->l -1;

	std::map<size_t,std::map<size_t,std::vector<reldata>> > mymapF;
	std::map<size_t,std::map<size_t,std::vector<reldata>> > mymapR;

	bam1_t *b = NULL;

	//first loop over all reads(these have the same chr/pos, and group these into queues that are pertile,perlib,pereverything)
	for(int i=0;i<q->l;i++) {
		//fprintf(stderr,"i:%d/%d\n",i,q->l);
		b = q->d[i];
		if(0&&!(b->core.flag &BAM_FDUP)){//never do this,
			if(fp)
				assert(sam_write1(fp, hdr,b)>=0);
			//noclusterdupcount++;
			continue;
		}
		if(bam_is_rev(b))
			plugin(mymapR,b,hdr);
		else
			plugin(mymapF,b,hdr);
	}
	std::vector<size_t> counter;
	plugout(mymapF,hdr,fp,fp2,counter);
	plugout(mymapR,hdr,fp,fp2,counter);
	for (int c=0; c<counter.size();c++){
		if(counter[c]>=histogram_l)
			tsktsk();
		histogram[counter[c]]++;
	}


	//assert(sam_write1(fp3, hdr, q->d[lrand48() %q->l])>=0); //<- this one prints a random read as the represent of the dups
	if(mymapF.size()>0){
		//std::vector<reldata> &re = mymapF.rbegin()->second;
		std::vector<reldata> &re = mymapF.rbegin()->second.rbegin()->second;
		if(nodupFP)
			assert(sam_write1(nodupFP, hdr,re[0].d));
		purecount++;
		//    fprintf(stderr,"%f len:%d purecount:%d\n",CMA,re[0].d->core.l_qseq,purecount);
		CMA = (re[0].d->core.l_qseq+(purecount-1)*CMA)/(1.0*purecount);
		//fprintf(stderr,"%f len:%d purecount:%d\n",CMA,re[0].d->core.l_qseq,purecount);
	}

	if(mymapR.size()>0){
		//std::vector<reldata> &re = mymapR.rbegin()->second;
		std::vector<reldata> &re = mymapR.rbegin()->second.rbegin()->second;
		if(nodupFP)
			assert(sam_write1(nodupFP, hdr,re[0].d));
		purecount++;
		CMA = (re[0].d->core.l_qseq+(purecount-1)*CMA)/(1.0*purecount);
	}

}


int usage(FILE *fp, int is_long_help)
{
	fprintf(fp,
			"\n"
			"Usage: ./superduper [options] <in.bam>|<in.sam>|<in.cram> \n"
			"\n"
			"Options:\n"
			// output options
			"  -b       output BAM\n"
			"  -C       output CRAM (requires -T)\n"
			"  -o FILE  output file name \n"
			"  -p INT   pixeldistance default 5000 \n"
			"  -T FILE  reference in the fastaformat (required from reading and writing crams)\n"
			"  -@ INT   Number of threads to use\n"
			"  -q INT   Mapping quality filter\n"
			"  -m       Discard unmapped reads (default off)\n"
			"  -w       Only calculate statistics (default off)\n"
			"  -v       Verbose mode\n"
			"  -a       Only use the SE part of the bam (default 1)\n"
			"Options for performing extraplation (mirrored from preseq)\n"
			"  -e       maximum extrapolation (default: 1e+10)\n"
			"  -s       step size in extrapolations (default: 1e+06)\n"
			"  -n       number of bootstraps (default: 100)\n"
			"  -c       level for confidence intervals (default: 0.95)\n"
			"  -x       maximum number of terms\n"
			"  -D       defects mode to extrapolate without testing for defects\n"
			"  -r       seed for random number generator\n"
			"\nThe Preseqpaper:\n   Daley, T., Smith, A. Predicting the molecular complexity of sequencing libraries.\n   Nat Methods 10, 325–327 (2013). https://doi.org/10.1038/nmeth.2375\n"
			// read filters
			);
	fprintf(fp,
			"\nNotes:\n"
			"\n"
			"1. This program is usefull for splitting a sorted bam/cram into two files\n"
			"   a) file containing clusterduplicates\n"
			"   a) file without any clusterduplicates, but including otherkinds of duplicates\n"
			"\n"
			"  Example: \n"
			"\t ./superduper input.bam -o outfiles -p 5000\n"
			"\n"
			"  Details:\n"
			"  It loops  over input files, and reads with identical positions\n"
			"  are assumed to be duplicates. It stratifes the duplicates over tile, lane and library,\n"
			"  and uses the euclidian distance (sqrt(da^2+db^2)) to 'find' clusters. Clusters being defined\n"
			"  as a group of reads that are within pxdist to another read within the cluster\n"
			"  program assumes readsnames looks like: \'A00706:12:HGNY3DSXX:3:1110:11930:4867\'\n"
			"  assuming \'discarded:discarded:discared:lanenumber:tileinfo:xpos:ypos\'\n"
			"\n\n"
			"  Limitations(bugs):\n"
			"  1) Does not use strand info\n"
			"  2) Does not use the flag field (external duplicate level e.g flag 1024)\n"
			"  3) Does not work with SAM files.\n"
			"  4) Does not with with regions \n");

	return 0;
}

void parse_sequencingdata(char *fn_out,char *refName,char *fname,int stats_only,int nthreads,int mapped_only,int se_only,int mapq,char *onam3,FILE *fp){
	htsThreadPool p = {NULL, 0};
	samFile *in=NULL;
	samFile *out=NULL;
	samFile *out2=NULL;
	samFile *nodupFP=NULL;

	char onam1[2048]="";
	char onam2[2048]="";
	char onam4[2048]="";

	strcat(onam1,fn_out);
	strcat(onam2,fn_out);
	strcat(onam4,fn_out);

	if(out_mode[1]=='b'){
		strcat(onam1,".noClusterDuplicates.bam");
		strcat(onam2,".onlyClusterDuplicates.bam");
		strcat(onam4,".pure.bam");
	}else{
		strcat(onam1,".noClusterDuplicates.cram");
		strcat(onam2,".onlyClusterDuplicates.cram");
		strcat(onam4,".pure.cram");
	}

	if(refName){
		char *ref =(char*) malloc(10 + strlen(refName) + 1);
		sprintf(ref, "reference=%s", refName);
		hts_opt_add((hts_opt **)&dingding2->specific,ref);
		free(ref);
	}

	if(strstr(fname,".cram")!=NULL &&out_mode[1]=='c'&&refName==NULL){
		fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
		exit(0);
	}

	if(out_mode[1]=='c'&&refName==NULL){
		fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
		exit(0);
	}

	if((in=sam_open_format(fname,"r",dingding2))==NULL ){
		fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
		exit(0);
	}

	if(stats_only==0){
		if ((out = sam_open_format(onam1, out_mode, dingding2)) == 0) {
			fprintf(stderr,"Error opening file for writing\n");
			exit(0);
		}

		if ((out2 = sam_open_format(onam2, out_mode, dingding2)) == 0) {
			fprintf(stderr,"Error opening file for writing\n");
			exit(0);
		}

		if ((nodupFP = sam_open_format(onam4, out_mode, dingding2)) == 0) {
			fprintf(stderr,"Error opening file for writing\n");
			exit(0);
		}

		if(nthreads>1){
			if (!(p.pool = hts_tpool_init(nthreads))) {
				fprintf(stderr, "Error creating thread pool\n");
				exit(0);
			}
			hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
			if (out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
			if (out2) hts_set_opt(out2, HTS_OPT_THREAD_POOL, &p);
		}
	}

	queue_t *queue = init_queue_t(nreads_per_pos);  
	bam_hdr_t  *hdr = sam_hdr_read(in);
	if(stats_only==0){
		assert(sam_hdr_write(out, hdr) == 0);
		assert(sam_hdr_write(out2, hdr) == 0);
		assert(sam_hdr_write(nodupFP, hdr) == 0);
	}

	//purecount=noclusterdupcount=0;
	purecount=0;
	bam1_t *b = bam_init1();

	int ret;
	int refId=-1;
	//int libid;
	while(((ret=sam_read1(in,hdr,b)))>0){
		nproc++;
		if(mapped_only!=0){
			if(b->core.flag&4)
				continue;
		}
		if(se_only==1){
			if(b->core.flag&1)//if paired continue
				continue;
		}



		if(mapq!=-1 && b->core.qual<mapq)
			continue;
		//catch case where there is one read in queue, and the next read is a new position
		//then we simply write it to the output
		if(refId==-1||refId!=b->core.tid){
			refId=b->core.tid;
			//fprintf(stderr,"\t-> Now at Chromosome: %s\n",hdr->target_name[refId]);
		}
		if(queue->l==1 && queue->d[0]->core.pos!=b->core.pos){
			histogram[1]++;
			if(out)
				assert(sam_write1(out, hdr, queue->d[0])>=0); //write into the file containing the pcrdups+normal reads
			if(nodupFP)
				assert(sam_write1(nodupFP, hdr, queue->d[0])>=0);//<- writeinto the file without any dups
			purecount++;
			CMA = (queue->d[0]->core.l_qseq+(purecount-1)*CMA)/(1.0*purecount);
			//noclusterdupcount++;
			queue->l =0;
		}

		if(queue->l>1 &&(queue->d[0]->core.tid!=b->core.tid ||(queue->d[0]->core.pos!=b->core.pos))){
			//fprintf(stderr,"calling do_magic\n");
			do_magic(queue,hdr,out,out2,nodupFP);
			queue->l =0;
		}

		if(queue->l==queue->m)
			realloc_queue(queue);

		assert(bam_copy1(queue->d[queue->l++],b)!=NULL);
		}

		if(queue->l==1){
			histogram[1]++;
			if(out)
				assert(sam_write1(out, hdr, queue->d[0])>=0); //write into the file containing the pcrdups+normal reads
			if(nodupFP)
				assert(sam_write1(nodupFP, hdr, queue->d[0])>=0);//<- writeinto the file without any dups
			purecount++;
			CMA = (queue->d[0]->core.l_qseq+(purecount-1)*CMA)/(1.0*purecount);
			//noclusterdupcount++;
		}else{
			//fprintf(stderr,"calling do_magic\n");
			do_magic(queue,hdr,out,out2,nodupFP);
		}
		queue->l=0;
		if(out)
			assert(sam_close(out)==0);
		if(out2)
			assert(sam_close(out2)==0);
		if(nodupFP)
			assert(sam_close(nodupFP)==0);
		if(in)
			assert(sam_close(in)==0);
		for(int i=0;i<queue->m;i++)
			bam_destroy1(queue->d[i]);

		free(queue->d);
		free(queue);
		//  bam_hdr_destroy(hdr);
		for(aMap::iterator it=char2int.begin();it!=char2int.end();it++)
			free(it->first);

		delete [] mystr;
		bam_destroy1(b);
		hts_opt_free((hts_opt *)dingding2->specific);
		free(dingding2);
		if(stats_only==0)
			fprintf(stderr,"    Dumpingfiles:\t\'%s\'\n\t\t\t\'%s\'\n\t\t\t\'%s\'\n\t\t\t\'%s\'\n",onam1,onam2,onam3,onam4);
		else
			fprintf(stderr,"    Dumpingfiles:\t\'%s\'\n",onam3);
		free(fn_out);
		free(fname);
		fprintf(stderr,
				"    Reads processed: %lu\n"
				"    Total duplicates: %lu\n"
				"    Cluster duplicates: %lu\n"
				"    Pcr duplicates: %d\n"
				"%lu\t%lu\t%f\n"
				,nproc,totaldups,clustdups,pcrdups,purecount,CMA);
				//"    Pcr duplicates: %lu\n"


		fprintf(fp,
				"#    Reads processed: %lu\n"
				"#    Total duplicates: %lu\n"
				"#    Cluster duplicates: %lu\n"
				"#    Pcr duplicates: %lu\n"
				"%lu\t%lu\t%f\n"
				,nproc,totaldups,clustdups,pcrdups,purecount,CMA);

	}


	int main(int argc, char **argv){
		histogram = new size_t [histogram_l];
		for(int i=0;i<histogram_l;i++)
			histogram[i] = 0;
		double max_extrapolation = 1.0e10;
		double step_size = 1e6;
		size_t bootstraps = 100;
		double c_level = 0.95;
		size_t orig_max_terms = 100;
		int DEFECTS = 0;
		int VERBOSE = 0;
		unsigned long int seed = 0;


		clock_t t=clock();
		time_t t2=time(NULL);

		char *fname,*refName;

		FILE *fp = NULL;
		FILE *fphist = NULL;
		FILE *fptable = NULL;
		fname=refName=NULL;
		char *fn_out = NULL;
		int c;
		int nthreads = 1;

		int mapq =-1;
		int mapped_only = 0;
		int se_only = 1;
		int stats_only = 0;
		char *histfile=NULL;
		if(argc==1){
			usage(stdout,0);
			return 0;
		}
		//fix these
		static struct option lopts[] = {
			{"add", 1, 0, 0},
			{"append", 0, 0, 0},
			{"delete", 1, 0, 0},
			{"verbose", 0, 0, 0},
			{"create", 1, 0, 'c'},
			{"file", 1, 0, 0},
			{NULL, 0, NULL, 0}
		};

		while ((c = getopt_long(argc, argv,
						"bCo:T:p:@:q:mwe:s:n:c:x:D:r:a:H:v",
						lopts, NULL)) >= 0) {
			switch (c) {
				case 'b': out_mode[1] = 'b'; break;
				case 'C': out_mode[1] = 'c'; break;
				case 'T': refName = strdup(optarg); break;
				case 'o': fn_out = strdup(optarg); break;
				case 'p': pxdist = atof(optarg); break;
				case '@': nthreads = atoi(optarg); break;
				case 'a': se_only = atoi(optarg); break;
				case 'q': mapq = atoi(optarg); break;
				case 'm': mapped_only = 1; break;
				case 'w': stats_only = 1; break;
				case 'e': max_extrapolation = atof(optarg); break;
				case 's': step_size = atof(optarg); break;
				case 'n': bootstraps = atoi(optarg); break;
				case 'c': c_level = atof(optarg); break;
				case 'x': orig_max_terms = atoi(optarg); break;
				case 'D': DEFECTS = atoi(optarg); break;
				case 'H': histfile = strdup(optarg); break;
				case 'v': VERBOSE = 1; break;
				case '?':
						  if (optopt == '?') {  // '-?' appeared on command line
							  return usage(stdout,0);
						  } else {
							  if (optopt) { // Bad short option
								  fprintf(stdout,"./superduper invalid option -- '%c'\n", optopt);
							  } else { // Bad long option
								  // Do our best.  There is no good solution to finding
								  // out what the bad option was.
								  // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
								  if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
									  fprintf(stdout,"./superduper unrecognised option '%s'\n",argv[optind - 1]);
								  }
							  }
							  return 0;//usage(stderr, 0);
						  }
				default:
						  fname = strdup(optarg);
						  fprintf(stderr,"assinging: %s to fname:%s\n",optarg,fname);
						  break;
			}
		}
		if(optind<argc)
			fname = strdup(argv[optind]);

		if(!fname&&!histfile){
			fprintf(stderr,"\t-> No input file specified\n");
			usage(stdout,0);
			return 0;
		}

		if(!fn_out){
			fprintf(stderr,"\t-> No output file specified\n");
			usage(stdout,0);
			return 0;
		}
		char onamhist[2048]="";
		char onamtable[2048]="";
		char onam3[2048]="";
		strcat(onam3,fn_out);

		snprintf(onamhist,2048,"%s.hist.txt",fn_out);
		snprintf(onamtable,2048,"%s.table.txt",fn_out);
		if ((fphist = fopen(onamhist, "wb")) == NULL) {
			fprintf(stderr,"Error opening file for writing\n");
			return 1;
		}
		strcat(onam3,".dupstat.txt");
		if ((fp = fopen(onam3, "wb")) == NULL) {
			fprintf(stderr,"Error opening file for writing\n");
			return 1;
		}  
		fprintf(stderr,"./superduper refName:%s fname:%s out_mode:%s pxdist:%f nthread:%d mapped_only:%d mapq:%d\nmax_extrap:%f step:%f boot:%lu c_lev:%f max_term:%lu defect:%d verbose:%d seed:%lu se_only:%d\n",refName,fname,out_mode,pxdist,nthreads,mapped_only,mapq,max_extrapolation,step_size,bootstraps,c_level,orig_max_terms,DEFECTS,VERBOSE,seed,se_only);
		fprintf(fp,"./superduper refName:%s fname:%s out_mode:%s pxdist:%f nthread:%d mapped_only:%d mapq:%d\nmax_extrap:%f step:%f boot:%lu c_lev:%f max_term:%lu defect:%d verbose:%d seed:%lu se_only:%d\n",refName,fname,out_mode,pxdist,nthreads,mapped_only,mapq,max_extrapolation,step_size,bootstraps,c_level,orig_max_terms,DEFECTS,VERBOSE,seed,se_only);
#ifdef __WITH_GSL__
		fprintf(stderr,"\t-> Using GSL library functions\n");
		fprintf(fp,"\t-> Using GSL library functions\n");
#else
		fprintf(stderr,"\t-> Using std=c++11 library functions\n");
		fprintf(fp,"\t-> Using std=c++11 library functions\n");
#endif

		if(fname)
			parse_sequencingdata(fn_out,refName,fname,stats_only,nthreads,mapped_only,se_only,mapq,onam3,fp);
		if(histfile){
			FILE *histfile_fp = NULL;
			histfile_fp = fopen(histfile,"rb");
			assert(histfile_fp);
			char histbuf[4096];
			while(fgets(histbuf,4096,histfile_fp)){
				size_t bin = atol(strtok(histbuf,"\t\t\n "));
				size_t value = atol(strtok(NULL,"\t\t\n "));
				if(bin>=histogram_l)
					tsktsk();
				//      fprintf(stderr,"bin: %lu value:%lu\n",bin,value);
				histogram[bin] = value;
			}
			fclose(histfile_fp);
		}

		int last=0;
		for(int i=0;i<histogram_l;i++)
			if(histogram[i])
				last=i;
		std::vector<double> to_preseq;
		for(int i=0;i<=last;i++){
			fprintf(fphist,"%d\t%lu\n",i,histogram[i]);
			//fprintf(stderr,"\n\nHISTOGRAM:\n%d\t%lu\n",i,histogram[i]);
			to_preseq.push_back(histogram[i]);
		}
		fclose(fphist);
		int lc_extrap(std::vector<double> &counts_hist,char *nam,double max_extrapolation, double step_size, size_t bootstraps, double c_level,size_t orig_max_terms, int DEFECTS,int VERBOSE, unsigned long int seed);
		//TODO for downstream enable this
		lc_extrap(to_preseq,onamtable,max_extrapolation,step_size,bootstraps,c_level,orig_max_terms,DEFECTS,VERBOSE,seed);
		fclose(fp);
		fprintf(stderr,
				"\t[ALL done] cpu-time used =  %.2f sec\n"
				"\t[ALL done] walltime used =  %.2f sec\n"
				,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));  
		fprintf(fp,
				"#[ALL done] cpu-time used =  %.2f sec\n"
				"#[ALL done] walltime used =  %.2f sec\n"
				,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));
		return 0;
	}

