#include "sam_funcs.h"

#undef __MY_DEBUG__
#undef __MY_DEBUG_READ__

#ifdef __MY_DEBUG__
#define __DEBUG_CHR__ "chr11"
#define __DEBUG_POS__ 75279846
#endif

void snpAnalysis(char* bamFileName, char* posFileName, char* methFileName, HashNode** hashTable, char** chrSeqArray, int* chrLen, int* chrDone, int vQualMin, int nLayerMin, int nLayerMax, float vSnpRate, float vSnpPerBase, unsigned int mapqThr) 
{
	int i, j, off, idx, len;
	int m, n, cnt, iread;

	MapRecord* record = (MapRecord*)malloc(sizeof(MapRecord));
	record->qname = (char*)malloc(sizeof(char) * 200);
	record->chrome = (char*)malloc(sizeof(char) * 50);
	record->cigar = (char*)malloc(sizeof(char) * 50);
	record->seq = (char*)malloc(sizeof(char) * 1000);
	record->seqBuf = (char*)malloc(sizeof(char) * 1000);
	record->qual = (char*)malloc(sizeof(char) * 50);
	record->qualBuf = (char*)malloc(sizeof(char) * 200);
	record->comBuf = (char*)malloc(sizeof(char) * 1000);
	
	unsigned short *w_A = NULL;
	unsigned short *w_T = NULL;
	unsigned short *w_C = NULL;
	unsigned short *w_G = NULL;
	unsigned short *c_A = NULL;
	unsigned short *c_T = NULL;
	unsigned short *c_C = NULL;
	unsigned short *c_G = NULL;
	
	unsigned int *w_Aq = NULL;
	unsigned int *w_Tq = NULL;
	unsigned int *w_Cq = NULL;
	unsigned int *w_Gq = NULL;
	unsigned int *c_Aq = NULL;
	unsigned int *c_Tq = NULL;
	unsigned int *c_Cq = NULL;
	unsigned int *c_Gq = NULL;
	
	unsigned short *w_An = NULL;
	unsigned short *w_Tn = NULL;
	unsigned short *w_Cn = NULL;
	unsigned short *w_Gn = NULL;
	unsigned short *c_An = NULL;
	unsigned short *c_Tn = NULL;
	unsigned short *c_Cn = NULL;
	unsigned short *c_Gn = NULL;
	
	unsigned int *w_Q = NULL;
	unsigned int *c_Q = NULL;
	
	unsigned int *w_Mm = NULL;
	unsigned int *w_Mc = NULL;
	unsigned int *w_Mq = NULL;
	unsigned int *c_Mm = NULL;
	unsigned int *c_Mc = NULL;
	unsigned int *c_Mq = NULL;
	
	int rowCnt = 0;
	int finCnt = 0;
	char curChr[50] = "chr123456789";
	uint8_t *s, *t;
    bam_header_t *header;
	
	bamFile in = bam_open(bamFileName, "r");
    if(in == NULL) {
		fprintf(stderr, "Cannot open bam file!\n");
		exit(1);
	}

	FILE* posFptr = fopen(posFileName, "w");
	if(posFptr == NULL) {
		fprintf(stderr, "Could not open file %s!\n", posFileName);
		exit(1);
	}
	fprintf(posFptr, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENOTYPE\tFREQUENCY\tNumber_of_watson[A,T,C,G]\tNumber_of_crick[A,T,C,G]\tMean_Quality_of_Watson[A,T,C,G]\tMean_Quality_of_Crick[A,T,C,G]\n");

	FILE* methFptr = fopen(methFileName, "w");
	if(methFptr == NULL) {
		fprintf(stderr, "Could not open file %s!\n", methFileName);
		exit(1);
	}
	fprintf(methFptr, "#CHROM\tPOS\tCONTEXT\tWatson METH\tWatson COVERAGE\tWatson QUAL\tCrick METH\tCrick COVERAGE\tCrick QUAL\n");
	
    bam1_t* b = bam_init1();
    if(b == NULL) {
		fprintf(stderr, "Cannot init bam structure!\n");
		exit(1);
	}
	
    header = bam_header_read(in);
	while(bam_read1(in, b) >= 0) {
		// Parse record
		if(parse_buffer(header, b, record, mapqThr) == 1) 
			continue;
		
		// Check chrome
		if(strcmp(curChr, record->chrome) != 0) {
			// Save old chrome statics results
			if(rowCnt > 0) {
				// Update
				finCnt = rowCnt;
				// Print
				print_meth(methFptr, len, curChr, w_Mm, w_Mc, w_Mq, c_Mm, c_Mc, c_Mq);
				print_snp(posFptr, chrSeqArray, idx, len, nLayerMin, nLayerMax, vSnpRate, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_An, w_Tn, w_Cn, w_Gn, c_An, c_Tn, c_Cn, c_Gn, w_Q, c_Q);
				// Memory gathering for x_X
				free(w_A);
				free(w_T);
				free(w_C);
				free(w_G);
				free(c_A);
				free(c_T);
				free(c_C);
				free(c_G);
				// Memory gathering for x_Xq
				free(w_Aq);
				free(w_Tq);
				free(w_Cq);
				free(w_Gq);
				free(c_Aq);
				free(c_Tq);
				free(c_Cq);
				free(c_Gq);
				// Memory gathering for x_Xn
				free(w_An);
				free(w_Tn);
				free(w_Cn);
				free(w_Gn);
				free(c_An);
				free(c_Tn);
				free(c_Cn);
				free(c_Gn);
				// Memory gathering for x_Q
				free(w_Q);
				free(c_Q);
				// Memory gathering for x_Mx
				free(w_Mm);
				free(w_Mc);
				free(w_Mq);
				free(c_Mm);
				free(c_Mc);
				free(c_Mq);
			}
            
			// Update current chrome
			strcpy(curChr, record->chrome);
            
			// Get chrome length
			idx = hash_table_lookup(hashTable, curChr);
			if(idx == -1) {
				fprintf(stderr, "%s not found in chrome name array.\n", curChr);
				exit(1);
			}
			len = chrLen[idx];
			
            // Check if chrome has been processed
            if(chrDone[idx] != 0) {
                fprintf(stderr, "%s has already been processed. The bam file is not sorted.\n", curChr);
				exit(1);
            }
            chrDone[idx] = 1;
            
			// Memory alloction for x_X
			if(!(w_A = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_A of %s.\n", curChr);
				exit(1);
			}
			if(!(w_T = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_T of %s.\n", curChr);
				exit(1);
			}
			if(!(w_C = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_C of %s.\n", curChr);
				exit(1);
			}
			if(!(w_G = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_G of %s.\n", curChr);
				exit(1);
			}
			if(!(c_A = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_A of %s.\n", curChr);
				exit(1);
			}
			if(!(c_T = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_T of %s.\n", curChr);
				exit(1);
			}
			if(!(c_C = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_C of %s.\n", curChr);
				exit(1);
			}
			if(!(c_G = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_G of %s.\n", curChr);
				exit(1);
			}
			
			// Memory alloction for x_Xq
			if(!(w_Aq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Aq of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Tq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Tq of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Cq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Cq of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Gq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Gq of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Aq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Aq of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Tq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Tq of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Cq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Cq of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Gq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Gq of %s.\n", curChr);
				exit(1);
			}
			
			// Memory alloction for x_Xn
			if(!(w_An = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_An of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Tn = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_Tn of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Cn = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_Cn of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Gn = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for w_Gn of %s.\n", curChr);
				exit(1);
			}
			if(!(c_An = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_An of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Tn = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_Tn of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Cn = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_Cn of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Gn = (unsigned short*)calloc(len, sizeof(unsigned short)))) {
				fprintf(stderr, "Not enough memory for c_Gn of %s.\n", curChr);
				exit(1);
			}
			
			// Memory alloction for x_Q
			if(!(w_Q = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Q of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Q = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Q of %s.\n", curChr);
				exit(1);
			}
			
			// Memory alloction for x_Mx
			if(!(w_Mm = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Mm of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Mc = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Mc of %s.\n", curChr);
				exit(1);
			}
			if(!(w_Mq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for w_Mq of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Mm = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Mm of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Mc = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Mc of %s.\n", curChr);
				exit(1);
			}
			if(!(c_Mq = (unsigned int*)calloc(len, sizeof(unsigned int)))) {
				fprintf(stderr, "Not enough memory for c_Mq of %s.\n", curChr);
				exit(1);
			}
		}

		/////////////////////////////////////////////////
		// Filter Step 1: Base quality
		/////////////////////////////////////////////////
		for(i = 0; i < record->len; i++) {
			if(record->seq[i] != 'N' && (unsigned short)(record->qual[i]-33) < vQualMin) 
				record->seq[i] = 'N';
		}
		/////////////////////////////////////////////////
		// Additional Step: Check and record CG methy
		/////////////////////////////////////////////////
		off = record->offset - 1;
		if((record->r12 == 1 && record->strand == '+') || (record->r12 == 2 && record->strand == '-')) {
			// Watson chain
			for(i = 0; i < record->len; i++) {
				if(record->seq[i] == 'C') {
					if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] == 'G') {
						w_Mm[off]++;
						w_Mc[off]++;
						w_Mq[off] += record->qual[i] - 33;
					}
				}
				else if(record->seq[i] == 'T') {
					if(chrSeqArray[idx][off] == 'C' && chrSeqArray[idx][off+1] == 'G') {
						w_Mc[off]++;
						w_Mq[off] += record->qual[i] - 33;
					}
				}
				
				off++;
			}
		}
		else {
			// Crick chain
			for(i = 0; i < record->len; i++) {
				if(record->seq[i] == 'G') {
					if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] == 'C') {
						c_Mm[off]++;
						c_Mc[off]++;
						c_Mq[off] += record->qual[i] - 33;
					}
				}
				else if(record->seq[i] == 'A') {
					if(chrSeqArray[idx][off] == 'G' && chrSeqArray[idx][off-1] == 'C') {
						c_Mc[off]++;
						c_Mq[off] += record->qual[i] - 33;
					}
				}
				
				off++;
			}
		}
		/////////////////////////////////////////////////
		// Filter Step 2: SNP per base
		/////////////////////////////////////////////////
		off = record->offset - 1;
		iread = (int)(record->len * vSnpPerBase + 0.5);
		if(iread < 1)
			iread = 1;
		cnt = 0;
		for(i = 0; i < record->len; i++) {
			if(record->strand == '+') {
				switch(record->seq[i]) {
					case 'A':
						if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'G')))
							cnt++;
						break;
					case 'T':
						if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'C')))
							cnt++;
						break;
					case 'C':
						if(chrSeqArray[idx][off] != 'C')
							cnt++;
						break;
					case 'G':
						if(chrSeqArray[idx][off] != 'G')
							cnt++;
						break;
				}
			}
			else {
				switch(record->seq[i]) {
					case 'A':
						if(chrSeqArray[idx][off] != 'A' && (!(record->r12 == 1 && chrSeqArray[idx][off] == 'G')))
							cnt++;
						break;
					case 'T':
						if(chrSeqArray[idx][off] != 'T' && (!(record->r12 == 2 && chrSeqArray[idx][off] == 'C')))
							cnt++;
						break;
					case 'C':
						if(chrSeqArray[idx][off] != 'C')
							cnt++;
						break;
					case 'G':
						if(chrSeqArray[idx][off] != 'G')
							cnt++;
						break;
				}
			}
			
			off++;
		}
		// Snp per read
		if(cnt <= iread) {
			off = record->offset - 1;
			for(i = 0; i < record->len; i++) {
				#ifdef __MY_DEBUG_READ__
				if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
					fprintf(stderr, "%c, %c. [r12:%d, ref:%c, pos:%d].\n", record->seq[i], record->strand, record->r12, chrSeqArray[idx][off], off);
					dispRecord(record);
				}
				#endif 
	
				if(record->strand == '+') {
					switch(record->seq[i]) {
						case 'A':
							if(record->r12 == 2 && chrSeqArray[idx][off] == 'G') {
								c_G[off]++;
								c_Gq[off] += (unsigned int)(record->qual[i] - 33);
								c_Gn[off]++;
								c_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_G.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								if(record->r12 == 1) {
									w_A[off]++;
									w_Aq[off] += (unsigned int)(record->qual[i] - 33);
									w_An[off]++;
									w_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_A.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
								else {
									c_A[off]++;
									c_Aq[off] += (unsigned int)(record->qual[i] - 33);
									c_An[off]++;
									c_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_A.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
							}
							break;
						case 'T':
							if(record->r12 == 1 && chrSeqArray[idx][off] == 'C') {
								w_C[off]++;
								w_Cq[off] += (unsigned int)(record->qual[i] - 33);
								w_Cn[off]++;
								w_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_C.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								if(record->r12 == 1) {
									w_T[off]++;
									w_Tq[off] += (unsigned int)(record->qual[i] - 33);
									w_Tn[off]++;
									w_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_T.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
								else {
									c_T[off]++;
									c_Tq[off] += (unsigned int)(record->qual[i] - 33);
									c_Tn[off]++;
									c_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_T.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
							}
							break;
						case 'C':
							if(record->r12 == 1) {
								w_C[off]++;
								w_Cq[off] += (unsigned int)(record->qual[i] - 33);
								w_Cn[off]++;
								w_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_C.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								c_C[off]++;
								c_Cq[off] += (unsigned int)(record->qual[i] - 33);
								c_Cn[off]++;
								c_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_C.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							break;
						case 'G':
							if(record->r12 == 1) {
								w_G[off]++;
								w_Gq[off] += (unsigned int)(record->qual[i] - 33);
								w_Gn[off]++;
								w_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_G.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								c_G[off]++;
								c_Gq[off] += (unsigned int)(record->qual[i] - 33);
								c_Gn[off]++;
								c_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_G.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							break;
					}
				}
				else {
					switch(record->seq[i]) {
						case 'A':
							if(record->r12 == 1 && chrSeqArray[idx][off] == 'G') {
								c_G[off]++;
								c_Gq[off] += (unsigned int)(record->qual[i] - 33);
								c_Gn[off]++;
								c_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_G.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								if(record->r12 == 2) {
									w_A[off]++;
									w_Aq[off] += (unsigned int)(record->qual[i] - 33);
									w_An[off]++;
									w_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_A.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
								else {
									c_A[off]++;
									c_Aq[off] += (unsigned int)(record->qual[i] - 33);
									c_An[off]++;
									c_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_A.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
							}
							break;
						case 'T':
							if(record->r12 == 2 && chrSeqArray[idx][off] == 'C') {
								w_C[off]++;
								w_Cq[off] += (unsigned int)(record->qual[i] - 33);
								w_Cn[off]++;
								w_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_C.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								if(record->r12 == 2) {
									w_T[off]++;
									w_Tq[off] += (unsigned int)(record->qual[i] - 33);
									w_Tn[off]++;
									w_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_T.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
								else {
									c_T[off]++;
									c_Tq[off] += (unsigned int)(record->qual[i] - 33);
									c_Tn[off]++;
									c_Q[off] += record->mapq;
									#ifdef __MY_DEBUG__
									if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
										fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_T.\n", record->qname, record->r12, record->strand, record->seq[i]);
									}
									#endif
								}
							}
							break;
						case 'C':
							if(record->r12 == 2) {
								w_C[off]++;
								w_Cq[off] += (unsigned int)(record->qual[i] - 33);
								w_Cn[off]++;
								w_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_C.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								c_C[off]++;
								c_Cq[off] += (unsigned int)(record->qual[i] - 33);
								c_Cn[off]++;
								c_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_C.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							break;
						case 'G':
							if(record->r12 == 2) {
								w_G[off]++;
								w_Gq[off] += (unsigned int)(record->qual[i] - 33);
								w_Gn[off]++;
								w_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to w_G.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							else {
								c_G[off]++;
								c_Gq[off] += (unsigned int)(record->qual[i] - 33);
								c_Gn[off]++;
								c_Q[off] += record->mapq;
								#ifdef __MY_DEBUG__
								if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && off+1 == __DEBUG_POS__) {
									fprintf(stderr, "read %s [r12:%d, %c, %c] contributed to c_G.\n", record->qname, record->r12, record->strand, record->seq[i]);
								}
								#endif
							}
							break;
					}
				}
				
				off++;
			}
		}
		else {
			#ifdef __MY_DEBUG__
			if(strcmp(record->chrome, __DEBUG_CHR__) == 0 && record->offset <= __DEBUG_POS__ && record->offset+record->len >= __DEBUG_POS__) {
				fprintf(stderr, "read discarded because snp-ratio overflow. [threshold %d, actually %d.]\n", iread, cnt);
				#ifdef __MY_DEBUG_READ__
				dispRecord(record);
				#endif
			}
			#endif
		}
		// Update row counter
		rowCnt++;
		if(rowCnt % 100000 == 0)
			fprintf(stderr, "%d records have been dealed.\n", rowCnt);
	}
	
	// Last batch
	if(finCnt < rowCnt) {
		// Print
        print_meth(methFptr, len, curChr, w_Mm, w_Mc, w_Mq, c_Mm, c_Mc, c_Mq);
		print_snp(posFptr, chrSeqArray, idx, len, nLayerMin, nLayerMax, vSnpRate, curChr, w_A, w_T, w_C, w_G, c_A, c_T, c_C, c_G, w_Aq, w_Tq, w_Cq, w_Gq, c_Aq, c_Tq, c_Cq, c_Gq, w_An, w_Tn, w_Cn, w_Gn, c_An, c_Tn, c_Cn, c_Gn, w_Q, c_Q);
		// Memory gathering for x_X
		free(w_A);
		free(w_T);
		free(w_C);
		free(w_G);
		free(c_A);
		free(c_T);
		free(c_C);
		free(c_G);
		// Memory gathering for x_Xq
		free(w_Aq);
		free(w_Tq);
		free(w_Cq);
		free(w_Gq);
		free(c_Aq);
		free(c_Tq);
		free(c_Cq);
		free(c_Gq);
		// Memory gathering for x_Xn
		free(w_An);
		free(w_Tn);
		free(w_Cn);
		free(w_Gn);
		free(c_An);
		free(c_Tn);
		free(c_Cn);
		free(c_Gn);
		// Memory gathering for x_Q
		free(w_Q);
		free(c_Q);
		// Memory gathering for x_Mx
		free(w_Mm);
		free(w_Mc);
		free(w_Mq);
		free(c_Mm);
		free(c_Mc);
		free(c_Mq);
	}
			
	bam_header_destroy(header);
	bam_close(in);
	bam_destroy1(b);
	fclose(posFptr);
	fclose(methFptr);
}

void print_meth(FILE* methFptr, int len, char* curChr, unsigned int* w_Mm, unsigned int* w_Mc, unsigned int* w_Mq, unsigned int* c_Mm, unsigned int* c_Mc, unsigned int* c_Mq)
{
	int i;
		
	// Record meth sites
	for(i = 0; i < len; i++) {
		if(w_Mm[i] + c_Mm[i+1] != 0) {
			if(w_Mm[i] == 0)
				fprintf(methFptr, "%s\t%d\tCG\t.\t.\t.\t%d\t%d\t%d\n", curChr, i + 1, c_Mm[i+1], c_Mc[i+1], (unsigned int)(c_Mq[i+1]/c_Mc[i+1]));
			else if(c_Mm[i+1] == 0)
				fprintf(methFptr, "%s\t%d\tCG\t%d\t%d\t%d\t.\t.\t.\n", curChr, i + 1, w_Mm[i], w_Mc[i], (unsigned int)(w_Mq[i]/w_Mc[i]));
			else 
				fprintf(methFptr, "%s\t%d\tCG\t%d\t%d\t%d\t%d\t%d\t%d\n", curChr, i + 1, w_Mm[i], w_Mc[i], (unsigned int)(w_Mq[i]/w_Mc[i]), c_Mm[i+1], c_Mc[i+1], (unsigned int)(c_Mq[i+1]/c_Mc[i+1]));
		}
	}
}

void print_snp(FILE* posFptr, char** chrSeqArray, int idx, int len, int nLayerMin, int nLayerMax, float vSnpRate, char* curChr, unsigned short *w_A, unsigned short *w_T, unsigned short *w_C, unsigned short *w_G, unsigned short *c_A, unsigned short *c_T, unsigned short *c_C, unsigned short *c_G, unsigned int *w_Aq, unsigned int *w_Tq, unsigned int *w_Cq, unsigned int *w_Gq, unsigned int *c_Aq, unsigned int *c_Tq, unsigned int *c_Cq, unsigned int *c_Gq, unsigned short *w_An, unsigned short *w_Tn, unsigned short *w_Cn, unsigned short *w_Gn, unsigned short *c_An, unsigned short *c_Tn, unsigned short *c_Cn, unsigned short *c_Gn, unsigned int *w_Q, unsigned int *c_Q)
{
	int i, j, m;
	int v, n;
	float vRate[3];
		
	// Record snp sites
	for(i = 0; i < len; i++) {
		m = w_A[i]+w_T[i]+w_C[i]+w_G[i]+c_A[i]+c_T[i]+c_C[i]+c_G[i];  
		if(m > nLayerMin && m < nLayerMax) {
			// No snp pt number
			switch(chrSeqArray[idx][i]) {
				case 'A':
					// A>C=(C+WsT)/total reads
					v = w_C[i] + c_C[i] + w_T[i];
					n = m;
					vRate[0] = (float)v/n;
					// A>T=T/total reads
					v = w_T[i] + c_T[i];
					n = m;
					vRate[1] = (float)v/n;
					// A>G=WsG/(WsG+WsA)
					v = w_G[i];
					n = w_G[i] + w_A[i];
					vRate[2] = (float)v/n;
					break;
				case 'T':
					// T>A=A/total reads
					v = w_A[i] + c_A[i];
					n = m;
					vRate[0] = (float)v/n;
					// T>C=CrC/(CrT+CrC)
					v = c_C[i];
					n = c_T[i] + c_C[i];;
					vRate[1] = (float)v/n;
					// T>G=(G+CrA)/total reads
					v = w_G[i] + c_G[i] + c_A[i];
					n = m;
					vRate[2] = (float)v/n;
					break;
				case 'C':
					// C>A=A/total reads
					v = w_A[i] + c_A[i];
					n = m;
					vRate[0] = (float)v/n;
					// C>T=CrT/(CrC+CrT)
					v = c_T[i];
					n = c_C[i] + c_T[i];
					vRate[1] = (float)v/n;
					// C>G=(G+CrA)/total reads
					v = w_G[i] + c_G[i] + c_A[i];
					n = m;
					vRate[2] = (float)v/n;
					break;
				case 'G':
					// G>A=WsA/(WsA+WsG)
					v = w_A[i];
					n = w_A[i] + w_G[i];
					vRate[0] = (float)v/n;
					// G>T=T/total reads
					v = w_T[i] + c_T[i];
					n = m;
					vRate[1] = (float)v/n;
					// G>C=(WsT+C)/total reads
					v = w_T[i] + w_C[i] + c_C[i];
					n = m;
					vRate[2] = (float)v/n;
					break;
			}
			// Filtering
			if(vRate[0] > vSnpRate || vRate[1] > vSnpRate || vRate[2] > vSnpRate) {	
				fprintf(posFptr, "%s\t%d\t%c\t", curChr, i + 1, chrSeqArray[idx][i]);
				fprintf(posFptr, "%u,%u,%u,%u\t", w_A[i], w_T[i], w_C[i], w_G[i]);
				fprintf(posFptr, "%u,%u,%u,%u\t", c_A[i], c_T[i], c_C[i], c_G[i]);
				#ifdef __MY_DEBUG_READ__
				fprintf(posFptr, "%u,%u,%u,%u\t", w_Aq[i], w_Tq[i], w_Cq[i], w_Gq[i]);
				fprintf(posFptr, "%u,%u,%u,%u\n", w_Aq[i], w_Tq[i], w_Cq[i], w_Gq[i]);
				#endif
				fprintf(posFptr, "%u,%u,%u,%u\t", (unsigned int)((float)w_Aq[i]/w_An[i]+0.5), (unsigned int)((float)w_Tq[i]/w_Tn[i]+0.5), (unsigned int)((float)w_Cq[i]/w_Cn[i]+0.5), (unsigned int)((float)w_Gq[i]/w_Gn[i]+0.5));
				fprintf(posFptr, "%u,%u,%u,%u\t", (unsigned int)((float)c_Aq[i]/c_An[i]+0.5), (unsigned int)((float)c_Tq[i]/c_Tn[i]+0.5), (unsigned int)((float)c_Cq[i]/c_Cn[i]+0.5), (unsigned int)((float)c_Gq[i]/c_Gn[i]+0.5));
				fprintf(posFptr, "%u\t", (unsigned int)((float)w_Q[i]/(w_An[i]+w_Tn[i]+w_Cn[i]+w_Gn[i])+0.5));
				fprintf(posFptr, "%u\n", (unsigned int)((float)c_Q[i]/(c_An[i]+c_Tn[i]+c_Cn[i]+c_Gn[i])+0.5));
			}
			else {
				#ifdef __MY_DEBUG__
				if(strcmp(curChr, __DEBUG_CHR__) == 0 && i+1 == __DEBUG_POS__) {
					fprintf(stderr, "pos %d of chr %s is not print because no-proper snp-ratio, actually %f.\n", i+1, curChr, (float)(m-n)/m);
				}
				#endif
			}
		}
		else {
			#ifdef __MY_DEBUG__
			if(strcmp(curChr, __DEBUG_CHR__) == 0 && i+1 == __DEBUG_POS__) {
				fprintf(stderr, "pos %d of chr %s is not print because no-proper layers, actually %d.\n", i+1, curChr, m);
			}
			#endif
		}
	}
}

int parse_buffer(bam_header_t *header, bam1_t *b, MapRecord* record, unsigned int mapqThr)
{
	int i, num, plen, flag;
	int seqPtr, seqBufPtr;
	char nchar;
	char *pdest;
		
	const bam1_core_t *c = &b->core;
	uint8_t *s = bam1_seq(b);
	uint8_t *t = bam1_qual(b);

	/////////////////////////////////
	// Load record
	/////////////////////////////////
	// {{
	// 1 - QNAME
	strcpy(record->qname, bam1_qname(b));
	// 2 - FLAG
	record->flag = c->flag;
	if((record->flag & 0x10) == 0)
		record->strand = '+';
	else
		record->strand = '-';
	if((record->flag & 0x40) != 0)
		record->r12 = 1;
	else
		record->r12 = 2;
	// 3 - RNAME
	strcpy(record->chrome, header->target_name[c->tid]);
	// 4 - POS
	record->offset = c->pos + 1;
	// 5 - MAPQ
	record->mapq = c->qual;
	if(c->qual < mapqThr) 
		return 1;
	// 6 - CIGAR
	for (i = 0; i < c->n_cigar; ++i)
		sprintf(record->cigar, "%d%c", bam1_cigar(b)[i]>>BAM_CIGAR_SHIFT, "MIDNSHP"[bam1_cigar(b)[i]&BAM_CIGAR_MASK]);
	// 10 - SEQ
	for(i = 0; i < c->l_qseq; ++i) 
		record->seqBuf[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
	record->seqBuf[i] = '\0';
	// 11 - QUAL
	if(t[0] == 0xff) {
		strcpy(record->qualBuf, "*");
	}
	else {
		for(i = 0; i < c->l_qseq; ++i) 
			record->qualBuf[i] = (char)(t[i] + 33);
	}
	// }}

	/////////////////////////////////
	// To upper
	/////////////////////////////////
	// {{
	for(i = 0; i < strlen(record->seqBuf); i++)
		record->seqBuf[i] = toupper(record->seqBuf[i]);
	// }}
	
	/////////////////////////////////
	// Cigar analysis
	/////////////////////////////////
	// {{
	pdest = record->cigar;
	nchar = nxtChar(pdest, &plen);
	seqPtr = 0;
	seqBufPtr = 0;
	while(nchar != 'Z') {
		// Digit parsing
		for(i = 0; i < plen; i++) 
			record->comBuf[i] = pdest[i];
		record->comBuf[plen] = '\0';
		num = atoi(record->comBuf);
		// Sequence handling
		switch(nchar) {
			case 'M':
				for(i = 0; i < num; i++) {
					record->seq[seqPtr] = record->seqBuf[seqBufPtr];
					record->qual[seqPtr] = record->qualBuf[seqBufPtr];
					seqPtr++;
					seqBufPtr++;
				}
				break;
			case 'D':
				for(i = 0; i < num; i++) {
					record->seq[seqPtr] = 'N';
					seqPtr++;
				}
				break;
			case 'I':
				seqBufPtr = seqBufPtr + num;
				break;
			case 'S':
				for(i = 0; i < num; i++) {
					record->seq[seqPtr] = 'N';
					seqPtr++;
					seqBufPtr++;
				}
				break;
			default:
				fprintf(stderr, "Error in cigar parsing!\n");
				exit(1);
		}
		
		pdest = pdest + plen + 1;
		nchar = nxtChar(pdest, &plen);
	}
	record->len = seqPtr;
	// }}

	return 0;
}

char nxtChar(char* cigar, int* len)
{
	// "len" is the length of digit part
	int tLen, mLen;
	char mChar;
	char* pdest;

	mLen = 1000;
	mChar = 'Z';

	pdest = strstr(cigar, "M");
	if(pdest != NULL) { 
		tLen = (int)(pdest - cigar);
		if(tLen < mLen) {
			mLen = tLen;
			mChar = 'M';
		}
	}
		
	pdest = strstr(cigar, "D");
	if(pdest != NULL) { 
		tLen = (int)(pdest - cigar);
		if(tLen < mLen) {
			mLen = tLen;
			mChar = 'D';
		}
	}
	
	pdest = strstr(cigar, "I");
	if(pdest != NULL) { 
		tLen = (int)(pdest - cigar);
		if(tLen < mLen) {
			mLen = tLen;
			mChar = 'I';
		}
	}
	
	pdest = strstr(cigar, "S");
	if(pdest != NULL) { 
		tLen = (int)(pdest - cigar);
		if(tLen < mLen) {
			mLen = tLen;
			mChar = 'S';
		}
	}
	
	if(mLen < 1000) 
		*len = mLen;
		
	return mChar;
}

void dispRecord(MapRecord* record)
{
	int i;
	
	for(i = 0; i < 40; i++)
		fprintf(stderr, "*");
	fprintf(stderr, "\n");
	
	fprintf(stderr, "qname: %s\n", record->qname);
	fprintf(stderr, "strand: %c\n", record->strand);
	fprintf(stderr, "chrome: %s\n", record->chrome);
	fprintf(stderr, "flag: %d\n", record->flag);
	fprintf(stderr, "offset: %d\n", record->offset);
	fprintf(stderr, "mapq: %d\n", record->mapq);
	fprintf(stderr, "cigar: %s\n", record->cigar);
	fprintf(stderr, "seq: %s\n", record->seq);
	fprintf(stderr, "qual: %s\n", record->qual);
	fprintf(stderr, "len: %d\n", record->len);
	fprintf(stderr, "r12: %d\n", record->r12);

	for(i = 0; i < 40; i++)
		fprintf(stderr, "*");
	fprintf(stderr, "\n");
}

int seqReverse(char* seq, int seqLen)
{
	int i, j;
	char tmp;
	
	// Step 1, Reverse
	for(i = 0; i < seqLen / 2; i++) {
		tmp = seq[i];
		seq[i] = seq[seqLen - 1 - i];
		seq[seqLen - 1 - i] = tmp;
	}
	// Step 2, Change
	for(i = 0; i < seqLen; i++) {
		switch(seq[i]) {
			case 'A':
				seq[i] = 'T';
				break;
			case 'T':
				seq[i] = 'A';
				break;
			case 'C':
				seq[i] = 'G';
				break;
			case 'G':
				seq[i] = 'C';
				break;
			case 'N':
				seq[i] = 'N';
				break;
			default:
				fprintf(stderr, "Sequence content error!\n");
				for(j = 0; j < seqLen; j++)
					fprintf(stderr, "%c", seq[j]);
				fprintf(stderr, "\n");
				return 1;
		}
	}
	
	return 0;
}
