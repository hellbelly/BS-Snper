#ifndef _CHROME_FUNCS_H
#define _CHROME_FUNCS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "hash_funcs.h"

#define INITIAL_CHR_NUM 10000
#define BUFSIZE 1000000

void init_chrome_hash(HashNode** hashTable, int* hash_table_size, char** chrName, int chrCnt);
void init_chrome_name_len(char* refSeqFile, int* chrCnt, char** chrName, int* chrLen);
void init_chrome_seq(HashNode** hashTable, char* refSeqFile, char** chrSeqArray, int* chrLen);

#endif
