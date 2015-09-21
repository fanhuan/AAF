#ifndef __FILE_READER_RJ_H
#define __FILE_READER_RJ_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "string.h"
#include "vector.h"

/**
 * Sequence IO
 */

typedef struct {
	String name;
	String comment;
	String seq;
	String qual;
} Sequence;

typedef struct {
	FILE *file;
	char *filename;
} fr_file_t;

typedef struct {
	Vector *files;
	uint32_t fidx;
	char *buffer;
	int size;
	int capacity;
	int ptr;
	int last_brk;
	char line_breaker;
	char delimiter;
	String *line;
	String *vline;
	Vector *tabs;
} FileReader;

#define free_sequence(sequence) { if(sequence->name.string) free(sequence->name.string);\
	if(sequence->comment.string) free(sequence->comment.string);\
	if(sequence->seq.string) free(sequence->seq.string);\
	if(sequence->qual.string) free(sequence->qual.string);\
	free(sequence); }

FileReader* fopen_filereader(char *filename);

FileReader* fopen_filereader2(char *prefix, char *postfix);

FileReader* fopen_m_filereader(int n_file, char **file_names);

FileReader* stdin_filereader();

/**
 * Read characters from a copy of string
 */

FileReader* string_filereader(char *string);

void fclose_filereader(FileReader *fr);

int reset_filereader(FileReader *fr);

int fread_line(String *line, FileReader *fr);
int froll_back(FileReader *fr);

int fread_table(FileReader *fr);
#define get_col_vstr(fr, col) ((VirtualString*)get_vec_ref((fr)->tabs, col))
#define get_col_str(fr, col) ((VirtualString*)get_vec_ref((fr)->tabs, col))->string
#define get_col_len(fr, col) ((VirtualString*)get_vec_ref((fr)->tabs, col))->size

typedef struct {
	int is_fq;
	int avg_seq_len;
	int min_seq_len;
	int max_seq_len;
} SeqFileAttr;

void guess_seq_file(FileReader *fr, SeqFileAttr *attr);

#define FASTA_FLAG_NORMAL		0
#define FASTA_FLAG_NO_NAME		1
#define FASTA_FLAG_NO_SEQ		2

int fread_fasta_adv(Sequence **seq, FileReader *fr, int flag);

#define fread_fasta(seq, fr) fread_fasta_adv(seq, fr, FASTA_FLAG_NORMAL)

#define FASTQ_FLAG_NORMAL		0
#define FASTQ_FLAG_NO_NAME		1
#define FASTQ_FLAG_NO_SEQ		2
#define FASTQ_FLAG_NO_QUAL		4

int fread_fastq_adv(Sequence **seq, FileReader *fr, int flag);

#define fread_fastq(seq, fr) fread_fastq_adv(seq, fr, FASTQ_FLAG_NORMAL)

char * fread_all(FileReader *fr);

static inline void print_pretty_seq(FILE *out, String *seq, int line_width){
	char c;
	int i, j;
	i = 0;
	while(i < seq->size){
		j = i + line_width;
		if(j > seq->size) j = seq->size;
		c  = seq->string[j];
		seq->string[j] = '\0';
		fprintf(out, "%s\n", seq->string + i);
		seq->string[j] = c;
		i = j;
	}
}

static inline FILE* open_file_for_read(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	file = fopen(full_name, "r");
	if(file == NULL) fprintf(stderr, "Cannot open file: %s\n", full_name);
	return file;
}

static inline FILE* open_file_for_write(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	file = fopen(full_name, "w+");
	if(file == NULL) fprintf(stderr, "Cannot open file: %s\n", full_name);
	return file;
}

static inline FILE* open_file_for_append(char *name, char *suffix){
	char *full_name;
	FILE *file;
	if(suffix == NULL){
		full_name = name;
	} else {
		full_name = (char*)alloca(strlen(name) + strlen(suffix) + 1);
		memcpy(full_name, name, strlen(name));
		memcpy(full_name + strlen(name), suffix, strlen(suffix) + 1);
	}
	file = fopen(full_name, "a+");
	if(file == NULL) fprintf(stderr, "Cannot open file: %s\n", full_name);
	return file;
}

typedef struct {
	FILE *file;
	void *buffer;
	int buf_off, buf_size, buf_cap;
} BufferedInputFile;

static inline BufferedInputFile* init_bif(FILE *file, int buf_size){
	BufferedInputFile *bif;
	bif = malloc(sizeof(BufferedInputFile));
	bif->file = file;
	bif->buf_off = bif->buf_size = 0;
	bif->buf_cap = buf_size;
	bif->buffer = malloc(buf_size);
	return bif;
}

static inline BufferedInputFile* open_bif(char *filename){
	FILE *file;
	if((file = fopen(filename, "r+")) == NULL){
		return NULL;
	}
	return init_bif(file, 1024);
}

static inline BufferedInputFile* open_bif2(char *filename, char *suffix){
	FILE *file;
	char *name;
	name = alloca(strlen(filename) + strlen(suffix) + 1);
	strcpy(name, filename);
	strcat(name, suffix);
	if((file = fopen(name, "r+")) == NULL){
		return NULL;
	}
	return init_bif(file, 1024);
}

static inline int64_t read_bif(BufferedInputFile *bif, void *data, int64_t size){
	int64_t i, t, ori_size;
	ori_size = size;
	while(size){
		if(bif->buf_size - bif->buf_off >= size){
			for(i=0;i<size;i++) ((unsigned char*)data)[i] = *((unsigned char*)bif->buffer + bif->buf_off + i);
			bif->buf_off += size;
			size = 0;
			break;
		} else if(bif->buf_off < bif->buf_size){
			t = bif->buf_size - bif->buf_off;
			for(i=0;i<t;i++) ((unsigned char*)data)[i] = *((unsigned char*)bif->buffer + bif->buf_off + i);
			data += t;
			size -= t;
			bif->buf_off = bif->buf_size;
		} else {
			bif->buf_size = fread(bif->buffer, 1, bif->buf_cap, bif->file);
			bif->buf_off = 0;
			if(bif->buf_size == 0) break;
		}
	}
	return ori_size - size;
}

static inline void close_bif(BufferedInputFile *bif){
	fclose(bif->file);
	free(bif->buffer);
	free(bif);
}

#endif
