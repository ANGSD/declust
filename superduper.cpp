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

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

double pxdist=5000;

typedef struct{
  bam1_t **d;//data
  unsigned l;//at pos
  unsigned m;//maxpos;
}queue_t;

size_t totaldups=0;
size_t pcrdups=0;
size_t clustdups=0;
size_t nproc=0;

int nreads_per_pos=4;//assuming this is the number reads per pos. If larger the program will reallocate

char out_mode[5]="wb";

queue_t *init_queue_t(int l){
  //  fprintf(stderr,"initializing queue with: %d elements\n",l);
  queue_t *ret =(queue_t *) malloc(sizeof(queue_t));
  ret->d =(bam1_t **) malloc(l*sizeof(bam1_t*));
  ret->l=0;
  ret->m=l;
  for(int i=0;i<ret->m;i++){
    //    fprintf(stderr,"queue[%d] init\n",i);
    ret->d[i] = bam_init1();
  }
  return ret;
}



void realloc_queue(queue_t *q){
  //  fprintf(stderr,"reallcing from: q:%d to q:%d\n",q->m,2*q->m);

  for(int i=0;0&&i<q->l;i++)
    fprintf(stderr,"inqueu[%d].pos:%d\n",i,q->d[i]->core.pos);

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
    fprintf(stderr,"onqueu[%d].pos:%d\n",i,q->d[i]->core.pos);
  
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

        // Valid until next query
        static char LB_text[1024];
        for (cp = LB; *cp && *cp != '\t' && *cp != '\n'; cp++)
            ;
        strncpy(LB_text, LB, MIN(cp-LB, 1023));
        LB_text[MIN(cp-LB, 1023)] = 0;

        // Return it; valid until the next query.
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

aMap char2int;

typedef struct{
  bam1_t *d;
  int xs;
  int ys;
}reldata;//<-releavant data


char *mystr =new char[2048];

double euc_dist(reldata &a,reldata &b){
  double dx=a.xs-b.xs;
  double dy=a.ys-b.ys;
  double dist = sqrt(dx*dx+dy*dy);
  

}

void print_clusters(std::vector<std::vector<int> > &clusters){
  fprintf(stderr,"[print clusters]dingding cluster.size():%lu\n",clusters.size());
  for(int i=0;i<clusters.size();i++){
    fprintf(stderr,"[print clusters] cluster %i/%lu\n",i,clusters.size() );
    std::vector<int> &tmp = clusters[i];
    for(int j=0;j<tmp.size();j++)
      fprintf(stderr,"[print clusters] ]%d:%d:%d\n",i,j,tmp[j]);
  }
}

void do_magic(queue_t *q,bam_hdr_t *hdr,samFile *fp,samFile *fp2){
  totaldups += q->l;
  //  fprintf(stderr,"do_magic queue->l:%d queue->m:%d chr:%d pos:%d\n",q->l,q->m,q->d[0]->core.tid,q->d[0]->core.pos);
  std::map<int,std::vector<reldata> > mymap;//<-this contains the library::lane::tile info as key. value is vector of reads,xpos,ypos
  bam1_t *b = NULL;

  //first loop over all reads(these have the same chr/pos, and group these into queues that are pertile,perlib,pereverything)
  for(int i=0;i<q->l;i++){
    b = q->d[i];
    if(0&&!(b->core.flag &BAM_FDUP)){
      assert(sam_write1(fp, hdr,b)>=0);      
      continue;
    }
    //    fprintf(stderr,"Spooling duplicates\n");
    reldata point;
    point.d=b;
    
    mystr = strncpy(mystr,bam_get_qname(b),2048);
    //    fprintf(stderr,"mystr: \'%s\'\n",mystr);
    strtok(mystr,"\n\t:");//machine
    strtok(NULL,"\n\t:");//runname
    strtok(NULL,"\n\t:");//flowcell
    int lane = atoi(strtok(NULL,"\n\t:"));
    int tile = atoi(strtok(NULL,"\n\t:"));
    point.xs = atoi(strtok(NULL,"\n\t:"));
    point.ys = atoi(strtok(NULL,"\n\t:"));
    char *lb = bam_get_library(hdr,b);
    if(char2int.find(lb)==char2int.end())
      char2int[strdup(lb)] = char2int.size();
    int libid = char2int.find(lb)->second;
    //[lib,lane,tile]
    //tile 4digits
    //lane 1digit
    //lib*10000+lane*1000+tile
    int key=libid*10000+lane*1000+tile;
    std::map<int,std::vector<reldata> >::iterator it =mymap.find(key);
    if(it==mymap.end()){
      std::vector<reldata> rd;
      rd.push_back(point);
      mymap[key]=rd;
    }else{
      it->second.push_back(point);
    }
      
  }

#if 0
  fprintf(stderr,"std::map.size:%lu\n",mymap.size());
  for(std::map<int,std::vector<reldata> >::iterator it=mymap.begin();it!=mymap.end();it++){
    fprintf(stderr,"key:%d\n",it->first);
    std::vector<reldata> &rd=it->second;
    for(int i=0;i<rd.size();i++)
      fprintf(stderr,"\tval: xs:%d ys:%d\n",rd[i].xs,rd[i].ys);
  }
#endif

  //looping over different libraries, lanes, etc
  for(std::map<int,std::vector<reldata> >::iterator it=mymap.begin();it!=mymap.end();it++) {
    //    fprintf(stderr,"key:%d\n",it->first);
    std::vector<reldata> &rd=it->second;
    //only one read in lane,lib,tile
    if(rd.size()==1){
      pcrdups++;
      assert(sam_write1(fp, hdr,rd[0].d)>=0);      
      continue;
    }

    
    for(int i=0;0&i<rd.size();i++)
      fprintf(stderr,"\tcc key: %d/%d val: xs:%d ys:%d pos:%d\n",it->first,rd.size(),rd[i].xs,rd[i].ys,rd[i].d->core.pos+1);

    
    //case where two reads come from same lane,lib tile
    if(rd.size()==2){
      double dx=rd[0].xs-rd[1].xs;
      double dy=rd[0].ys-rd[1].ys;
      double dist = euc_dist(rd[0],rd[1]);
      // fprintf(stderr,"dist is:%f\n",dist);
      //      sd=math.sqrt(((coordsm[0]-coordsl[0])**2)+((coordsm[1]-coordsl[1])**2))
      if(dist>pxdist){
	assert(sam_write1(fp, hdr,rd[0].d)>=0);      
	assert(sam_write1(fp, hdr,rd[1].d)>=0);
	pcrdups +=2 ;
      }else{
	assert(sam_write1(fp2, hdr,rd[0].d)>=0);      
	assert(sam_write1(fp2, hdr,rd[1].d)>=0);
	clustdups +=2 ;
      }
      continue;
    }
    
    //case with more than 2 reads. We therefore find the center of all points. And find the difference from the center
    //old version
#if 0
    double dx =0;
    double dy =0;
    for(int i=0;i<rd.size();i++){
      dx += rd[i].xs;
      dy += rd[i].ys;
    }
    double cx = dx/(1.0*rd.size());
    double cy = dy/(1.0*rd.size());
    
    for(int i=0;i<rd.size();i++){
      double dx=rd[i].xs-cx;
      double dy=rd[i].ys-cy;
      double dist = sqrt(dx*dx+dy*dy);
      
      //      fprintf(stderr,"center is: %f %f pos:%d dist:%f\n",cx,cy,rd[0].d->core.pos+1,dist);
      if(dist>(pxdist/2.0)){
	assert(sam_write1(fp, hdr,rd[i].d)>=0);      
      }else
	assert(sam_write1(fp2, hdr,rd[i].d)>=0);      
      
    }
#endif
    
    std::vector<std::vector<int> > clusters;
    //    std::vector<int> tmp;tmp.push_back(0);
    //clusters.push_back(tmp);
    //    print_clusters(clusters);


    //loop over all reads, but the first one, this defines cluster1.
    for(int i=0;i<rd.size();i++) {
      //      print_clusters(clusters);
      //fprintf(stderr,"analysing rd:%d\n",i);
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
	fprintf(stderr,"i:%d %d\n",i,dingdongsong[i]);
#endif 
      //now dingdongsong contains a 0/1 array indicating which existing clusters it belongs to.
      int nclust=0;//counter for how many clusters
      for(int s=0;s<rd.size();s++){
	//	fprintf(stderr,"analysinz look up table:%d/%lu=%d\n",s,rd.size(),dingdongsong[s]);
	nclust += dingdongsong[s];
      }
      //      fprintf(stderr,"nclust: %d\n",nclust);
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
	clusters[keep[0]].push_back(i);
	//	print_clusters(clusters);
      }
    }

    //print_clusters(clusters);
    //    fprintf(stderr,"\t-> Flushing\n");
    for(int i=0;i<clusters.size();i++){
      std::vector<int> &tmp = clusters[i];
      if(tmp.size()>0){
	pcrdups++;
	assert(sam_write1(fp, hdr,rd[tmp[0]].d)>=0);
      }if(tmp.size()<2)
	continue;
      for(int j=1;j<tmp.size();j++){
	clustdups++;
	assert(sam_write1(fp2, hdr,rd[tmp[j]].d)>=0);
      }
    }    
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


int main(int argc, char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);

  char *fname,*refName;
  samFile *in=NULL;
  samFile *out=NULL;
  samFile *out2=NULL;
  FILE *fp = NULL;
  fname=refName=NULL;
  char *fn_out = NULL;
  int c;
  int nthreads = 1;
  htsThreadPool p = {NULL, 0};
  if(argc==1){
    usage(stdout,0);
    return 0;
  }
  
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
			  "bCo:T:p:@:",
			  lopts, NULL)) >= 0) {
    switch (c) {
        case 'b': out_mode[1] = 'b'; break;
        case 'C': out_mode[1] = 'c'; break;
        case 'T': refName = strdup(optarg); break;
        case 'o': fn_out = strdup(optarg); break;
        case 'p': pxdist = atof(optarg); break;
        case '@': nthreads = atoi(optarg); break;
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
      fprintf(stderr,"adsadsfasdf\n");
      fname = strdup(optarg);
      fprintf(stderr,"assinging: %s to fname:%s\n",optarg,fname);
      break;
    }
  }
  fname = strdup(argv[optind]);

  if(!fname){
    fprintf(stderr,"\t-> No input file specified\n");
    usage(stdout,0);
    return 0;
  }

  if(!fn_out){
    fprintf(stderr,"\t-> No output file specified\n");
    usage(stdout,0);
    return 0;
  }
  
  fprintf(stderr,"./superduper refName:%s fname:%s out_mode:%s pxdist:%f nthread:%d\n",refName,fname,out_mode,pxdist,nthreads);
  
  if(refName){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    sprintf(ref, "reference=%s", refName);
    hts_opt_add((hts_opt **)&dingding2->specific,ref);
    free(ref);
  }

  if(strstr(fname,".cram")!=NULL &&out_mode[1]=='c'&&refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    return 0;
  }
  if(out_mode[1]=='c'&&refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    return 0;
  }
  
  if((in=sam_open_format(fname,"r",dingding2))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }

  char onam1[2048]="";
  char onam2[2048]="";
  char onam3[2048]="";
  strcat(onam1,fn_out);
  strcat(onam2,fn_out);
  strcat(onam3,fn_out);
  if(out_mode[1]=='b'){
    strcat(onam1,".noClusterDuplicates.bam");
    strcat(onam2,".onlyClusterDuplicates.bam");
  }else{
    strcat(onam1,".noClusterDuplicates.cram");
    strcat(onam2,".onlyClusterDuplicates.cram");
  }
  strcat(onam3,".dupstat.txt");
  if ((fp = fopen(onam3, "wb")) == NULL) {
    fprintf(stderr,"Error opening file for writing\n");
    return 1;
  }

  fprintf(fp,"#./superduper refName:%s fname:%s out_mode:%s pxdist:%f nthread:%d\n",refName,fname,out_mode,pxdist,nthreads);
  
  if ((out = sam_open_format(onam1, out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    return 1;
  }
  
  if ((out2 = sam_open_format(onam2, out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    return 1;
  }


  if(nthreads>1){
    if (!(p.pool = hts_tpool_init(nthreads))) {
      fprintf(stderr, "Error creating thread pool\n");
      return 0;
    }
    hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
    if (out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    if (out2) hts_set_opt(out2, HTS_OPT_THREAD_POOL, &p);
    
  }


  
  queue_t *queue = init_queue_t(nreads_per_pos);  
  bam_hdr_t  *hdr = sam_hdr_read(in);
  assert(sam_hdr_write(out, hdr) == 0);
  assert(sam_hdr_write(out2, hdr) == 0);
  
  bam1_t *b = bam_init1();

  int ret;
  while(((ret=sam_read1(in,hdr,b)))>0){
    nproc++;
    //catch case where there is one read in queue, and the next read is a new position
    //then we simply write it to the output
    if(queue->l==1 && queue->d[0]->core.pos!=b->core.pos){
      assert(sam_write1(out, hdr, queue->d[0])>=0);
      queue->l =0;
    }
    if(queue->l>1 &&(queue->d[0]->core.tid!=b->core.tid ||(queue->d[0]->core.pos!=b->core.pos))){
      //      fprintf(stderr,"calling do_magic\n");
      do_magic(queue,hdr,out,out2);
      queue->l =0;
    }

    if(queue->l==queue->m)
      realloc_queue(queue);

    bam_copy1(queue->d[queue->l++],b);
  }
  do_magic(queue,hdr,out,out2);
  queue->l;
  assert(sam_close(out)==0);
  assert(sam_close(out2)==0);
  assert(sam_close(in)==0);
  for(int i=0;i<queue->m;i++)
    bam_destroy1(queue->d[i]);
  free(queue->d);
  free(queue);
  bam_hdr_destroy(hdr);
  for(aMap::iterator it=char2int.begin();it!=char2int.end();it++)
    free(it->first);

  delete [] mystr;
  bam_destroy1(b);
  hts_opt_free((hts_opt *)dingding2->specific);
  free(dingding2);
  fprintf(stderr,"    Dumpingfiles:\t\'%s\'\n\t\t\t\'%s\'\n\t\t\t\'%s\'\n",onam1,onam2,onam3);
  free(fn_out);
  free(fname);
  fprintf(stderr,
	  "    reads processed: %lu\n"
          "    total duplicates: %lu\n"
	  "    cluster duplicates: %lu\n"
	  "    noncluster duplicates: %lu\n"
	  ,nproc,totaldups,clustdups,pcrdups);
  fprintf(fp,
	  "    reads processed: %lu\n"
          "    total duplicates: %lu\n"
	  "    cluster duplicates: %lu\n"
	  "    noncluster duplicates: %lu\n"
	  ,nproc,totaldups,clustdups,pcrdups);

  fprintf(stderr,
	  "\t[ALL done] cpu-time used =  %.2f sec\n"
	  "\t[ALL done] walltime used =  %.2f sec\n"
	  ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));  
  fprintf(fp,
	  "\t#[ALL done] cpu-time used =  %.2f sec\n"
	  "\t#[ALL done] walltime used =  %.2f sec\n"
	  ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));  
  fclose(fp);
  return 0;
}
