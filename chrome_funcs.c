#include "chrome_funcs.h"

void init_chrome_hash(HashNode** hashTable, int* hash_table_size, char* chrLenFile, int* chrCnt)
{
	
	FILE* fPtr = fopen(chrLenFile, "r");
	if(!fPtr) { 
		printf("File %s open error.\n", chrLenFile);
		exit(1);
	}
	fclose(fPtr);
	
	char* token = NULL;
	char readBuffer[500];
	char seps[] = " \t\n\r";
	*chrCnt = 0;
	fPtr = fopen(chrLenFile, "r");
	while(fgets(readBuffer, 500, fPtr)) {
		token = strtok(readBuffer, seps);
		hash_table_insert(hashTable, hash_table_size, token, *chrCnt);
		(*chrCnt)++;
	}
	fclose(fPtr);
	printf("%d chromes loaded.\n", *chrCnt);
}

void init_chrome_name_len(HashNode** hashTable, char* chrLenFile, int chrCnt, char** chrName, int* chrLen)
{
	int i, j, seqPtr;
	
	FILE* fPtr = fopen(chrLenFile, "r");
	if(!fPtr) { 
		printf("File %s open error.\n", chrLenFile);
		exit(1);
	}
	fclose(fPtr);
	
	char* token = NULL;
	char readBuffer[500];
	char seps[] = " \t\n\r";

	fPtr = fopen(chrLenFile, "r");
	for(i = 0; i < chrCnt; i++) {
		fgets(readBuffer, 500, fPtr);
		
		token = strtok(readBuffer, seps);
		seqPtr = hash_table_lookup(hashTable, token);
		strcpy(chrName[seqPtr], token);
		
		token = strtok(NULL, seps);
		j = atoi(token);
		chrLen[seqPtr] = j;
	}
	fclose(fPtr);
}

void init_chrome_seq(HashNode** hashTable, char* refSeqFile, char** chrSeqArray, int* chrLen)
{
	
	FILE* fPtr = fopen(refSeqFile, "r");
	if(!fPtr) { 
		printf("File %s open error.\n", refSeqFile);
		exit(1);
	}
	fclose(fPtr);
	
	char* readBuffer = chrSeqArray[0];
	int readBufferPtr = 0;
	int readSeqPtr = 0;
	int lines = 0;
	
	char* token;
	char seps[] = "\t\n\r";
	char chrName[50];
	
	fPtr = fopen(refSeqFile, "r");
	while(fgets(readBuffer + readBufferPtr, 500, fPtr)) {
		if(*(readBuffer + readBufferPtr) == '>') {
			if(lines > 0) {
				if(readBufferPtr != chrLen[readSeqPtr] ) {
					printf("Length of seq %s error! %d vs %d!\n", chrName, readBufferPtr, chrLen[readSeqPtr]);
					exit(1);
				}
			}
			token = strtok(readBuffer + readBufferPtr + 1, seps);	
			strcpy(chrName, token);
			readSeqPtr = hash_table_lookup(hashTable, chrName);
			if(readSeqPtr == -1) {
				printf("Error find sequence with name %s!\n", chrName);
				exit(1);
			}
			readBuffer = chrSeqArray[readSeqPtr];
			readBufferPtr = 0;	
		}
		else {
			readBufferPtr += strlen(readBuffer + readBufferPtr) - 1;
			if(readBufferPtr > chrLen[readSeqPtr]) {
				printf("Length of seq %s error at line %d! %d vs %d!\n", chrName, lines, readBufferPtr, chrLen[readSeqPtr]);
				exit(1);
			}
		}
		
		lines++;
	}
	if(readBufferPtr != chrLen[readSeqPtr] ) {
		printf("Length of seq %s error! %d vs %d!\n", chrName, readBufferPtr, chrLen[readSeqPtr]);
		exit(1);
	}
	fclose(fPtr);
}

