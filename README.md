## Installation

```bash
git clone https://github.com/ANGSD/superduper.git
cd superduper
make
```

## Usage

```
./superduper [options] <in.bam>|<in.sam>|<in.cram> 

Options:
  -b       Output BAM
  -C       Output CRAM (requires reference fasta; use -T)
  -o FILE  Output file name
  -p INT   Pixel distance (default: 12000)
  -T FILE  Reference in the fasta format (required for reading and writing crams)
  -@ INT   Number of threads to use
  -q INT   Mapping quality filter (default: off)
  -m       Discard unmapped reads (default: off)
  -w       Only calculate statistics (default: off)
  -W       Calculate additional statistics (default: 0, off)
					Output summary table and frequency distribution tables
								MSC - Mean sequence complexity
								MGC - Mean GC content
								MFS - Mean fragment size
								SCD - Sequence complexity distribution
								GCD - GC content distribution
								FSD - Fragment size distribution
					Example: To extract sequence complexity distribution, use:
						`grep ^SCD out.dupstat.txt | cut -f 2-`
  -v       Verbose mode
  -a       Only use the single end part of the bam (default: 1 (enabled), use -a 0 to disable)
  -X INT	Sequence complexity filter, discard read if complexity<INT (0-100, default: off)
  -G INT	Maximum GC content allowed, discard read if GC content>INT (0-100, default: off)
  -l INT	Minimum read length allowed, discard read if read length<INT (default: off)
  -L INT	Maximum read length allowed, discard read if read length>INT (default: off)

Options for performing extraplation (mirrored from preseq)
  -e       maximum extrapolation (default: 1e+10)
  -s       step size in extrapolations (default: 1e+06)
  -n       number of bootstraps (default: 100)
  -c       level for confidence intervals (default: 0.95)
  -x       maximum number of terms
  -D       defects mode to extrapolate without testing for defects
  -r       seed for random number generator

```


Notes:

1. This program is useful for splitting a sorted bam/cram into two files
   1) file containing cluster duplicates
   2) file without any cluster duplicates, but including other kinds of duplicates


  Details:
  It loops  over input files, and reads with identical positions
  are assumed to be duplicates. It stratifes the duplicates over tiles and lanes
  and uses the euclidian distance (sqrt(da^2 +db^2)) to 'find' clusters. Clusters being defined
  as a group of reads that are within pxdist to another read within the cluster
  program assumes read names (QNAME) look like: 'A00706:12:HGNY3DSXX:3:1110:11930:4867'
  assuming 'discarded:discarded:discared:lanenumber:tileinfo:xpos:ypos'

  tileinfo is a 4 digit number including the identifiers for surface, swath and tile
  for given tileinfo 1234, 1=surface, 2=swath, 34=tile
	1			2		34
	-			-		--
	surface		swath	tile

  For more details, see Illumina NovaSeq 6000 Sequencing System Guide 
  Document #1000000019358v14 Material #20023471



## Output
### Summary statistics and frequency distribution tables

MSC - Mean sequence complexity
MGC - Mean GC content
MFS - Mean fragment size
SCD - Sequence complexity distribution
GCD - GC content distribution
FSD - Fragment size distribution

  
To extract the sequence complexity distribution (SCD), use:

`grep ^SCD out.dupstat.txt | cut -f 2-`


## Limitations
  
  1) Does not use strand info
  2) Does not use the flag field (external duplicate level e.g flag 1024)
  3) Does not work with SAM files.
  4) Does not with with regions 


## How to cite

superduper ...

The Preseq paper:
   Daley, T., Smith, A. Predicting the molecular complexity of sequencing libraries.
   Nat Methods 10, 325–327 (2013). https://doi.org/10.1038/nmeth.2375


## License
GPLv3?
