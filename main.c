#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hash_funcs.h"
#include "chrome_funcs.h"
#include "sam_funcs.h"



int main(int argc, char **argv)
{
	int i, j;
	long cnt;
	FILE *bamFptr, *posFptr;

	//////////////////////////////////////////////////////////////////////////////
	// Parse command args
	//////////////////////////////////////////////////////////////////////////////
	// {{
	if(argc != 13) {
		printf("Not enough args! Argc = %d.\n", argc);
		exit(1);
	}

	char* refSeqFile = argv[1];
	char* bamFileName = argv[2];
	char* snpFileName = argv[3];
	char* methCgFileName = argv[4];
    char* methChgFileName = argv[5];
    char* methChhFileName = argv[6];

    int vQualMin = atoi(argv[7]);
	int nLayerMin = atoi(argv[8]);
	int nLayerMax = atoi(argv[9]);
	float vSnpRate = atof(argv[10]);
	float vSnpPerBase = atof(argv[11]);
	unsigned int mapqThr = atoi(argv[12]);

	fprintf(stderr, "refSeqFile = %s.\n", refSeqFile);
	fprintf(stderr, "bamFileName = %s.\n", bamFileName);
	fprintf(stderr, "snpFileName = %s.\n", snpFileName);
	fprintf(stderr, "methCgFileName = %s.\n", methCgFileName);
    fprintf(stderr, "methChgFileName = %s.\n", methChgFileName);
    fprintf(stderr, "methChhFileName = %s.\n", methChhFileName);

	fprintf(stderr, "vQualMin = %d.\n", vQualMin);
	fprintf(stderr, "nLayerMax = %d.\n", nLayerMax);
	fprintf(stderr, "vSnpRate = %f.\n", vSnpRate);
	fprintf(stderr, "vSnpPerBase = %f.\n", vSnpPerBase);
	fprintf(stderr, "mapqThr = %d.\n", mapqThr);
	// }}

	//////////////////////////////////////////////////////////////////////////////
	// Load reference sequence & Init arrays
	//////////////////////////////////////////////////////////////////////////////
	// {{
	// Init chrome name and length array
    int chrCnt = INITIAL_CHR_NUM;
	int* chrLen = (int*)malloc(sizeof(int) * chrCnt);
	char** chrName = (char**)malloc(sizeof(char*) * chrCnt);
	for(i = 0; i < chrCnt; i++)
		chrName[i] = (char*)malloc(sizeof(char) * 200);
	init_chrome_name_len(refSeqFile, &chrCnt, chrName, chrLen);
	// Realloc
	chrLen = (int*)realloc(chrLen, sizeof(int) * chrCnt);
	for(i = chrCnt; i < INITIAL_CHR_NUM; i++) {
		free(chrName[i]);
	}
	chrName = (char**)realloc(chrName, sizeof(char*) * chrCnt);
	// Cal seq memory
	cnt = 0;
	for(i = 0; i < chrCnt; i++)
		cnt += chrLen[i];
	fprintf(stderr, "Init chrome name-len array completed, total %ld bp of %d chromosomes.\n", cnt, chrCnt);
    // Hash table
	int hash_table_size;
	HashNode** hashTable = (HashNode**)malloc(sizeof(HashNode*) * HASH_TABLE_MAX_SIZE);
	hash_table_init(hashTable, &hash_table_size);
	// Init chrome name-idx hash table
	init_chrome_hash(hashTable, &hash_table_size, chrName, chrCnt);
	fprintf(stderr, "Init chrome name-idx hash table completed.\n");
	// Init chrome seq array
	char** chrSeqArray = (char**)malloc(sizeof(char*) * chrCnt);
	for(i = 0; i < chrCnt; i++)
		chrSeqArray[i] = (char*)malloc(sizeof(char) * chrLen[i] + 1000);
	init_chrome_seq(hashTable, refSeqFile, chrSeqArray, chrLen);
    for(i = 0; i < chrCnt; i++) {
		if(!(chrSeqArray[i] = (char*)realloc((void*)(chrSeqArray[i]), sizeof(char) * chrLen[i]))) {
			fprintf(stderr, "Not enough memory for chrSeqArray.\n");
			exit(1);
		}
		for(j = 0; j < chrLen[i]; j++)
			chrSeqArray[i][j] = toupper(chrSeqArray[i][j]);
	}
	fprintf(stderr, "Init chrome seq array completed.\n");

	int* chrDone = (int*)malloc(sizeof(int) * chrCnt);
	chrDone[i] = 0;
	// }}

	//////////////////////////////////////////////////////////////////////////////
	// SNP process
	//////////////////////////////////////////////////////////////////////////////
	// {{
	snpAnalysis(bamFileName, snpFileName, methCgFileName, methChgFileName, methChhFileName, hashTable, chrSeqArray, chrLen, chrCnt, vQualMin, nLayerMin, nLayerMax, vSnpRate, vSnpPerBase, mapqThr);
	fprintf(stderr, "SNP process completed.\n");
	// }}

	return 0;
}
