library(ggplot2)
library(dplyr)

lc_pattern="_lcextrap.txt"
c_pattern="_ccurve.txt"
preseq_dir="/run/user/1000/gvfs/sftp:host=ku/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/superduper/WORKINGDIR/WD-2110/hugh-superduperlist/preseq"

#set lc_nrow to X get the first X M reads from the lc_curve extrapolation result
asDirPath <- function(string){
	if(!endsWith(string,"/"))string=paste0(string,"/")
	if(!dir.exists(string))warning("Directory does not exist, will exit!")
	return(string)
}

multisample_preseq_lc_c <- function(preseq_dir,lc_pattern,c_pattern,lc_nrow=50){

	preseq_dir=asDirPath(preseq_dir)
	lc_files=list.files(path=preseq_dir, pattern=lc_pattern)
	all_data=NULL

	for(lc_file in lc_files){
		lc=read.table(file=lc_file, header=T)
		sample=strsplit(lc_file, lc_pattern)[[1]][1]
		c=read.table(file=list.files(path=preseq_dir, pattern=glob2rx(paste0(sample, "*",c_pattern))), skip=2, col.names=c("TOTAL_READS","DISTINCT_READS_CC"))
		lc<-lc[1:lc_nrow+1,]
		lc_c<-merge(lc, c, by="TOTAL_READS", all.x=T)
		lc_c$sample=sample
		all_data<-rbind(all_data, lc_c)
	}
	return(all_data)
}
preseq_dir="/run/user/1000/gvfs/sftp:host=ku/science/willerslev/users-shared/science-snm-willerslev-pfs488/projects/superduper/WORKINGDIR/WD-2110/hugh-superduperlist/preseq"

sample_id="024045P_sea_HM0575_HM0575USS21A_s2mq_U_ver10_CTFON_210507-GRCh37"

preseq_lc_c <- function(preseq_dir, sample_id, lc_pattern="_lcextrap.txt", c_pattern="_ccurve.txt", lc_nrow=50){

	preseq_dir=asDirPath(preseq_dir)
	lc_file=list.files(path=preseq_dir, pattern=paste0(sample_id, lc_pattern))
	lc=read.table(file=lc_file, header=T)
	sample=strsplit(lc_file, lc_pattern)[[1]][1]
	c=read.table(file=list.files(path=preseq_dir, pattern=glob2rx(paste0(sample,c_pattern))), skip=2, col.names=c("TOTAL_READS","DISTINCT_READS_CC"))
	lc<-lc[1:lc_nrow+1,]
	lc_c<-merge(lc, c, by="TOTAL_READS", all.x=T)
	lc_c$sample=sample

	return(lc_c)
}

plot_lc_c <- function(data){
	sample=unique(data$sample)
	ggplot(data, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT))+
	geom_line(linetype="dotted", size=0.75)+
	geom_line(aes(x=TOTAL_READS, y=DISTINCT_READS_CC), size=1)+
	geom_ribbon(aes(x=TOTAL_READS, ymin=LOWER_0.95CI, ymax=UPPER_0.95CI), alpha=0.3)+
	labs(x="Total number of reads",
		 y="Expected distinct reads",
		 title="Preseq library complexity estimate using lc_extrap and c_curve methods",
		 subtitle=paste0("Sample: ", sample))+
theme_bw()+
theme(plot.title = element_text(face="bold", size="14"),
	  plot.subtitle = element_text(size="10"))+
coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE)
}

plot_lc_c(preseq_lc_c(preseq_dir,sample_id))

plot_lc_c(lc_c)


multisample_preseq_lc_c <- function(preseq_dir, lc_pattern="_lcextrap.txt", c_pattern="_ccurve.txt", lc_nrow=50){

	preseq_dir=asDirPath(preseq_dir)
	lc_files=list.files(path=preseq_dir, pattern=lc_pattern)
	all_data=NULL

	for(lc_file in lc_files){
		lc=read.table(file=lc_file, header=T)
		sample=strsplit(lc_file, lc_pattern)[[1]][1]
		c=read.table(file=list.files(path=preseq_dir, pattern=glob2rx(paste0(sample, "*",c_pattern))), skip=2, col.names=c("TOTAL_READS","DISTINCT_READS_CC"))
		lc<-lc[1:lc_nrow+1,]
		lc_c<-merge(lc, c, by="TOTAL_READS", all.x=T)
		lc_c$sample=sample
		all_data<-rbind(all_data, lc_c)
	}
	return(all_data)
}

d<-multisample_preseq_lc_c(preseq_dir,lc_pattern="_lcextrap.txt",c_pattern="_ccurve.txt")

multi_plot_lc_c <- function(data){
	ggplot(data, aes(x=TOTAL_READS, y=EXPECTED_DISTINCT, colour=sample)) +
	geom_point(size=0.5)+
	geom_line(aes(x=TOTAL_READS, y=DISTINCT_READS_CC), size=0.5) +
	#  geom_ribbon(aes(x=TOTAL_READS, ymin=LOWER_0.95CI, ymax=UPPER_0.95CI), alpha=0.3)+
	scale_colour_brewer(palette = 'Accent') +
	labs(title="Preseq library complexity estimates using lc_extrap and c_curve methods")+
	labs(x="Total number of reads",
		 y="Expected distinct reads",
		 title="Preseq library complexity estimate using lc_extrap and c_curve methods")+
theme_bw()+
theme(plot.title = element_text(face="bold", size="12"))+
theme(legend.position = 'none')
}
multi_plot_lc_c(d)


pdf("preseq_results.pdf")
print(plot_lc_c(lc_c))
print(multi_plot_lc_c(d))
dev.off()

