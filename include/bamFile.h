/* bamFile -- interface to binary alignment format files using Heng Li's samtools lib. */

#ifndef BAMFILE_H
#define BAMFILE_H

#include "dnaseq.h"
#include "dystring.h"

#ifdef USE_BAM

// bam.h is incomplete without _IOLIB set to 1, 2 or 3.  2 is used by Makefile.generic:
#ifndef _IOLIB
#define _IOLIB 2
#endif
#ifdef USE_HTS
#include "htslib/sam.h"
typedef samFile samfile_t;
typedef hts_idx_t bam_index_t;
typedef bam_hdr_t bam_header_t;
typedef int (*bam_fetch_f)(const bam1_t *bam, void *data, bam_hdr_t *header) ;
#define samopen(a,b,c) sam_open(a,b)
#define samclose(a) sam_close(a)
#define bam1_qname bam_get_qname
#define bam1_qual bam_get_qual
#define bam1_aux bam_get_aux
#define bam1_cigar bam_get_cigar
#define bam1_seq bam_get_seq
#define bam1_seqi bam_seqi
#define bam_nt16_rev_table seq_nt16_str
#define data_len l_data
#else
#include "bam.h"
#include "sam.h"
#endif

#else // no USE_BAM
typedef struct { } bam1_t;
typedef struct { } bam_index_t;
typedef struct { } samfile_t;
typedef int (*bam_fetch_f)(const bam1_t *b, void *data);

#define COMPILE_WITH_SAMTOOLS "%s: in order to use this functionality you must " \
    "install the samtools library (<A HREF=\"http://samtools.sourceforge.net\" " \
    "TARGET=_BLANK>http://samtools.sourceforge.net</A>) and recompile kent/src with " \
    "USE_BAM=1 in your environment " \
    "(see <A HREF=\"http://genomewiki.ucsc.edu/index.php/Build_Environment_Variables\" " \
    "TARGET=_BLANK>http://genomewiki.ucsc.edu/index.php/Build_Environment_Variables</A>)."

#endif // USE_BAM

struct bamChromInfo
    {
    struct bamChromInfo *next;
    char *name;		/* Chromosome name */
    bits32 size;	/* Chromosome size in bases */
    };

boolean bamFileExists(char *bamFileName);
/* Return TRUE if we can successfully open the bam file and its index file. */

void bamFileAndIndexMustExist(char *fileOrUrl);
/* Open both a bam file and its accompanying index or errAbort; this is what it
 * takes for diagnostic info to propagate up through errCatches in calling code. */

samfile_t *bamOpen(char *fileOrUrl, char **retBamFileName);
/* Return an open bam file as well as the filename of the bam. */

samfile_t *bamMustOpenLocal(char *fileName, char *mode, void *extraHeader);
/* Open up sam or bam file or die trying.  The mode parameter is 
 *    "r" - open SAM to read
 *    "rb" - open BAM to read
 *    "w" - open SAM to write
 *    "wb" - open BAM to write
 * The extraHeader is generally NULL in the read case, and the write case
 * contains a pointer to a bam_header_t with information about the header.
 * The implementation is just a wrapper around samopen from the samtools library
 * that aborts with error message if there's a problem with the open. */

#ifdef USE_HTS
void bamFetchAlreadyOpen(samfile_t *samfile, bam_hdr_t *header,  bam_index_t *idx, char *bamFileName, 
#else
void bamFetchAlreadyOpen(samfile_t *samfile, bam_index_t *idx, char *bamFileName, 
#endif
			 char *position, bam_fetch_f callbackFunc, void *callbackData);
/* With the open bam file, return items the same way with the callbacks as with bamFetch() */
/* except in this case use an already-open bam file and index (use bam_index_load and free() for */
/* the index). It seems a little strange to pass the filename in with the open bam, but */
/* it's just used to report errors. */

void bamFetchPlus(char *fileOrUrl, char *position, bam_fetch_f callbackFunc, void *callbackData,
		 samfile_t **pSamFile, char *refUrl, char *cacheDir);
/* Open the .bam file, fetch items in the seq:start-end position range,
 * and call callbackFunc on each bam item retrieved from the file plus callbackData.
 * This handles BAM files with "chr"-less sequence names, e.g. from Ensembl. 
 * The pSamFile parameter is optional.  If non-NULL it will be filled in, just for
 * the benefit of the callback function, with the open samFile.  
 * refUrl points to the place to grab CRAM reference sequences (if any)
 * cacheDir points to the directory in which CRAM reference sequences are cached */

void bamFetch(char *fileOrUrl, char *position, bam_fetch_f callbackFunc, void *callbackData,
	samfile_t **pSamFile);
/* Open the .bam file, fetch items in the seq:start-end position range,
 * and call callbackFunc on each bam item retrieved from the file plus callbackData.
 * This handles BAM files with "chr"-less sequence names, e.g. from Ensembl. 
 * The pSamFile parameter is optional.  If non-NULL it will be filled in, just for
 * the benefit of the callback function, with the open samFile.  */

void bamClose(samfile_t **pSamFile);
/* Close down a samefile_t */

boolean bamIsRc(const bam1_t *bam);
/* Return TRUE if alignment is on - strand. */

INLINE int bamUnpackCigarElement(unsigned int packed, char *retOp)
/* Given an unsigned int containing a number of bases and an offset into an
 * array of BAM-enhanced-CIGAR ASCII characters (operations), store operation 
 * char into *retOp (retOp must not be NULL) and return the number of bases. */
{
#ifdef USE_BAM
// decoding lifted from samtools bam.c bam_format1_core(), long may it remain stable:
#define BAM_DOT_C_OPCODE_STRING "MIDNSHP=X"
int n = packed>>BAM_CIGAR_SHIFT;
int opcode = packed & BAM_CIGAR_MASK;
if (opcode >= strlen(BAM_DOT_C_OPCODE_STRING))
    errAbort("bamUnpackCigarElement: unrecognized opcode %d. "
	     "(I only recognize 0..%lu [" BAM_DOT_C_OPCODE_STRING "])  "
	     "Perhaps samtools bam.c's bam_format1 encoding changed?  If so, update me.",
	     opcode, (unsigned long)(strlen(BAM_DOT_C_OPCODE_STRING)-1));
*retOp = BAM_DOT_C_OPCODE_STRING[opcode];
return n;
#else // no USE_BAM
errAbort(COMPILE_WITH_SAMTOOLS, "bamUnpackCigarElement");
return 0;
#endif// USE_BAM
}

void bamGetSoftClipping(const bam1_t *bam, int *retLow, int *retHigh, int *retClippedQLen);
/* If retLow is non-NULL, set it to the number of "soft-clipped" (skipped) bases at
 * the beginning of the query sequence and quality; likewise for retHigh at end.
 * For convenience, retClippedQLen is the original query length minus soft clipping
 * (and the length of the query sequence that will be returned). */

void bamUnpackQuerySequence(const bam1_t *bam, boolean useStrand, char *qSeq);
/* Fill in qSeq with the nucleotide sequence encoded in bam.  The BAM format 
 * reverse-complements query sequence when the alignment is on the - strand,
 * so if useStrand is given we rev-comp it back to restore the original query 
 * sequence. */

char *bamGetQuerySequence(const bam1_t *bam, boolean useStrand);
/* Return the nucleotide sequence encoded in bam.  The BAM format 
 * reverse-complements query sequence when the alignment is on the - strand,
 * so if useStrand is given we rev-comp it back to restore the original query 
 * sequence. */

UBYTE *bamGetQueryQuals(const bam1_t *bam, boolean useStrand);
/* Return the base quality scores encoded in bam as an array of ubytes. */

void bamUnpackCigar(const bam1_t *bam, struct dyString *dyCigar);
/* Unpack CIGAR string into dynamic string */

char *bamGetCigar(const bam1_t *bam);
/* Return a BAM-enhanced CIGAR string, decoded from the packed encoding in bam. */

void bamShowCigarEnglish(const bam1_t *bam);
/* Print out cigar in English e.g. "20 (mis)Match, 1 Deletion, 3 (mis)Match" */

void bamShowFlagsEnglish(const bam1_t *bam);
/* Print out flags in English, e.g. "Mate is on '-' strand; Properly paired". */

int bamGetTargetLength(const bam1_t *bam);
/* Tally up the alignment's length on the reference sequence from
 * bam's packed-int CIGAR representation. */

bam1_t *bamClone(const bam1_t *bam);
/* Return a newly allocated copy of bam. */

void bamShowTags(const bam1_t *bam);
/* Print out tags in HTML: bold key, no type indicator for brevity. */

char *bamGetTagString(const bam1_t *bam, char *tag, char *buf, size_t bufSize);
/* If bam's tags include the given 2-character tag, place the value into 
 * buf (zero-terminated, trunc'd if nec) and return a pointer to buf,
 * or NULL if tag is not present. */

void bamUnpackAux(const bam1_t *bam, struct dyString *dy);
/* Unpack the tag:type:val part of bam into dy */

struct bamChromInfo *bamChromList(samfile_t *fh);
/* Return list of chromosomes from bam header. We make no attempty to normalize chromosome names to UCSC format,
   so list may contain things like "1" for "chr1", "I" for "chrI", "MT" for "chrM" etc. */

void bamChromInfoFreeList(struct bamChromInfo **pList);
/* Free a list of dynamically allocated bamChromInfo's */

void samToBed(char *samIn, char *bedOut);
/* samToBed - Convert SAM file to a pretty simple minded bed file.. */

void samToOpenBed(char *samIn, FILE *f);
/* Like samToOpenBed, but the output is the already open file f. */

#endif//ndef BAMFILE_H
