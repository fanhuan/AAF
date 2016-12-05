#include "file_reader.h"
#include <stdlib.h>
#include <string.h>

FileReader* fopen_m_filereader(int n_file, char **filenames){
	FileReader *fr;
	fr_file_t *fc;
	char *cmd;
	int i;
	fr = (FileReader*)malloc(sizeof(FileReader));
	fr->files = init_vec(sizeof(fr_file_t), n_file);
	for(i=0;i<n_file;i++){
		fc = get_next_vec_ref(fr->files);
		if(filenames[i] == NULL || strcmp(filenames[i], "-") == 0){
			fc->file = stdin;
			fc->filename = NULL;
		} else if(strlen(filenames[i]) > 3 && strcmp(filenames[i] + strlen(filenames[i]) - 3, ".gz") == 0){
			cmd = (char*)malloc(strlen(filenames[i]) + 20);
			sprintf(cmd, "gzip -dc %s", filenames[i]);
			fc->filename = (char*)malloc(sizeof(char)* (strlen(filenames[i])+1));
			fc->file = popen(cmd, "r");
			free(cmd);
		} else if((fc->file = fopen(filenames[i], "r")) != NULL){
			fc->filename = (char*)malloc(sizeof(char)* (strlen(filenames[i])+1));
			strcpy(fc->filename, filenames[i]);
		} else {
			return NULL;
		}
	}
	fr->fidx  = 0;
	fr->ptr   = 0;
	fr->last_brk = 0;
	fr->size  = 0;
	fr->capacity = 512;
	fr->buffer = (char*)malloc(fr->capacity + 2);
	fr->line_breaker = '\n';
	fr->delimiter    = '\t';
	fr->line = init_string(81);
	fr->vline = NULL;
	fr->tabs = init_vec(sizeof(VirtualString), 12);
	return fr;
}

FileReader* fopen_filereader(char *filename){
	char *filenames[1];
	filenames[0] = filename;
	return fopen_m_filereader(1, filenames);
}

FileReader* fopen_filereader2(char *prefix, char *postfix){
	char *filename;
	filename = alloca(strlen(prefix) + strlen(postfix) + 1);
	filename[0] = 0;
	strcat(filename, prefix);
	strcat(filename, postfix);
	return fopen_filereader(filename);
}

FileReader* stdin_filereader(){
	return fopen_filereader(NULL);
}

FileReader* string_filereader(char *string){
	FileReader *fr = (FileReader*)malloc(sizeof(FileReader));
	fr->files = init_vec(sizeof(fr_file_t), 1);
	fr->fidx = 0;
	fr->ptr   = 0;
	fr->last_brk = 0;
	fr->size  = strlen(string);
	fr->capacity = fr->size;
	fr->buffer = string;
	fr->line_breaker = '\n';
	fr->delimiter    = '\t';
	fr->line = init_string(81);
	fr->vline = NULL;
	fr->tabs = init_vec(sizeof(VirtualString), 12);
	return fr;
}

void fclose_filereader(FileReader *fr){
	fr_file_t *fc;
	size_t i;
	for(i=0;i<vec_size(fr->files);i++){
		fc = get_vec_ref(fr->files, i);
		if(fc->file && fc->file != stdin){
			if(fc->filename && strlen(fc->filename) > 3 && strcmp(fc->filename + strlen(fc->filename) - 3, ".gz") == 0) pclose(fc->file);
			else if(fc->file != stdin) fclose(fc->file);
		}
		if(fc->filename) free(fc->filename);
	}
	free_vec(fr->files);
	if(fr->buffer != NULL) free(fr->buffer);
	fr->buffer = NULL;
	if(fr->line){ free_string(fr->line); }
	if(fr->vline){ free_string(fr->vline); }
	free_vec(fr->tabs);
	free(fr);
}

static inline int fr_fread(void *buf, size_t e_size, size_t size, FILE *in){
	size_t n;
	int c;
	if(in != stdin || e_size > 1) return fread(buf, e_size, size, in);
	n = 0;
	while(n < size){
		c = getchar();
		if(c == -1) break;
		else ((char*)buf)[n++] = c;
		if(c == '\n') break;
	}
	return n;
}

int fread_line2(String *line, FileReader *fr){
	int ret, last_ptr, n;
	ret = 0;
	last_ptr = fr->ptr;
	while(1){
		if(last_ptr < fr->size){
			while(last_ptr < fr->size){
				if(fr->buffer[last_ptr++] == fr->line_breaker){ ret = 1; break; }
			}
			if(ret == 1) break;
		} else if(fr->fidx < vec_size(fr->files)) {
			if(fr->ptr){
				memmove(fr->buffer, fr->buffer + fr->ptr, fr->size - fr->ptr);
				last_ptr -= fr->ptr;
				fr->size -= fr->ptr;
				fr->ptr = 0;
			}
			if(fr->size == fr->capacity){
				fr->capacity += 4 * 1024;
				fr->buffer = (char*)realloc(fr->buffer, fr->capacity + 2);
			}
			n = fr_fread(fr->buffer + fr->size, sizeof(char), fr->capacity - fr->size, ((fr_file_t*)get_vec_ref(fr->files, fr->fidx))->file);
			if(n == 0){
				fr->fidx ++;
			} else {
				fr->size += n;
			}
		} else {
			break;
		}
	}
	if(last_ptr > fr->ptr){
		append_string(line, fr->buffer + fr->ptr, last_ptr - fr->ptr - ret);
	} else ret = -1;
	fr->last_brk = fr->ptr;
	fr->ptr = last_ptr;
	return ret;
}

int fread_line(String *line, FileReader *fr){
	clear_string(line);
	if(fread_line2(line, fr) < 0){
		return -1;
	} else {
		return line->size;
	}
}

int froll_back(FileReader *fr){
	if(fr->last_brk >= fr->ptr) return 0;
	fr->ptr      = fr->last_brk;
	return 1;
}

int* init_delimiters(char *expr){
	int *delimiters, i, state, len;
	delimiters = (int*)malloc(sizeof(int) * 128);
	memset(delimiters, 0, sizeof(int) * 128);
	len = strlen(expr);
	state = 0;
	for(i=0;i<len;i++){
		if(expr[i] == '\\'){
			if(state){
				delimiters[(int)expr[i]] = 1;
				state = 0;
			} else {
				state = 1;
			}
		} else if(state){
			switch(expr[i]){
				case 't':
					delimiters['\t'] = 1;
					break;
				case 's':
					delimiters[' '] = 1;
					break;
				case 'n':
					delimiters['\n'] = 1;
					break;
				case 'r':
					delimiters['\n'] = 1;
					break;
				default:
					delimiters[(int)expr[i]] = 1;
			}
		} else {
			delimiters[(int)expr[i]] = 1;
		}
	}
	return delimiters;
}

int fread_table(FileReader *fr){
	VirtualString *vstr;
	int i, ret;
	if(fread_line(fr->line, fr) < 0) return -1;
	if(fr->vline == NULL){
		fr->vline = init_string(fr->line->size);
		append_string(fr->vline, fr->line->string, fr->line->size);
	} else {
		clear_string(fr->vline);
		append_string(fr->vline, fr->line->string, fr->line->size);
	}
	clear_vec(fr->tabs);
	ret = split_string(fr->vline, fr->delimiter, fr->tabs);
	for(i=1;i<ret;i++){
		vstr = get_vec_ref(fr->tabs, i);
		vstr->string[-1] = 0;
	}
	return ret;
}

int fread_fasta_adv(Sequence **seq_ptr, FileReader *fr, int fasta_flag){
	Sequence *seq;
	int i, n, flag;
	if(*seq_ptr == NULL){
		seq = (Sequence*)malloc(sizeof(Sequence));
		seq->name.string = seq->comment.string = seq->seq.string = seq->qual.string = NULL;
		seq->name.size = seq->comment.size = seq->seq.size = seq->qual.size = 0;
		seq->name.capacity = seq->comment.capacity = seq->seq.capacity = seq->qual.capacity = 0;
	} else {
		seq = *seq_ptr;
	}
	flag = 0;
	while((n = fread_line(fr->line, fr)) != -1){
		if(n && fr->line->string[0] == '>'){
			if(flag){
				froll_back(fr);
				break;
			}
			flag = 1;
			seq->name.size = 0;
			seq->comment.size = 0;
			seq->seq.size = 0;
			if((fasta_flag & FASTA_FLAG_NO_NAME) == 0){
				for(i=1;i<n;i++){
					switch(fr->line->string[i]){
						case ' ':
						case '\t':
						case '\r':
						case '\n':
						goto BREAK_OUT;
					}
				}
				BREAK_OUT:
				append_string(&(seq->name), fr->line->string + 1, i - 1);
				if(i + 1 < n) append_string(&(seq->comment), fr->line->string + i + 1, n - i - 1);
			}
		} else {
			if((fasta_flag & FASTA_FLAG_NO_SEQ) == 0){
				append_string(&(seq->seq), fr->line->string, n);
			}
			flag = 2;
		}
	}
	if(flag < 2){
		free_sequence(seq);
		*seq_ptr = NULL;
		clear_string(fr->line);
		return 0;
	} else {
		*seq_ptr = seq;
		return 1;
	}
}

int fread_fastq_adv(Sequence **seq_ptr, FileReader *fr, int fastq_flag){
	Sequence *seq;
	int i, n, flag;
	if(*seq_ptr == NULL){
		seq = (Sequence*)malloc(sizeof(Sequence));
		seq->name.string = seq->comment.string = seq->seq.string = seq->qual.string = NULL;
		seq->name.capacity = seq->comment.capacity = seq->seq.capacity = seq->qual.capacity = 0;
	} else {
		seq = *seq_ptr;
	}
	seq->name.size = seq->comment.size = seq->seq.size = seq->qual.size = 0;
	flag = 0;
	while(flag != 4 && (n = fread_line(fr->line, fr)) >= 0){
		switch(flag){
			case 0:
				if(fr->line->string[0] != '@') break;
				flag = 1;
				if(fastq_flag & FASTQ_FLAG_NO_NAME) break;
				for(i=1;i<n;i++) if(fr->line->string[i] == ' ' || fr->line->string[i] == '\t' || fr->line->string[i] == '\n') break;
				append_string(&seq->name, fr->line->string + 1, i - 1);
				if(i < n) append_string(&seq->comment, fr->line->string + i + 1, n - i - 1);
				break;
			case 1:
				flag = 2;
				if(fastq_flag & FASTQ_FLAG_NO_SEQ) break;
				append_string(&seq->seq, fr->line->string, n);
				break;
			case 2:
				if(fr->line->string[0] != '+') break;
				flag = 3;
				break;
			case 3:
				flag = 4;
				if(fastq_flag & FASTQ_FLAG_NO_QUAL) break;
				append_string(&seq->qual, fr->line->string, n);
				break;
		}
	}
	if(flag < 4){
		free_sequence(seq);
		*seq_ptr = NULL;
		clear_string(fr->line);
		return 0;
	} else {
		*seq_ptr = seq;
		return 1;
	}
}

int guess_seq_file_type(FileReader *fr){
	while(fread_line(fr->line, fr) != -1){
		if(fr->line->size == 0) continue;
		if(fr->line->string[0] == '#') continue;
		if(fr->line->string[0] == '>'){
			froll_back(fr);
			return 1;
		} else if(fr->line->string[0] == '@'){
			froll_back(fr);
			return 2;
		} else {
			froll_back(fr);
			return 0;
		}
	}
	return 0;
}

void guess_seq_file(FileReader *fr, SeqFileAttr *attr){
	int n_seq, size;
	Sequence *seq;
	attr->is_fq = (guess_seq_file_type(fr) == 2);
	attr->min_seq_len = 0x7FFFFFFF;
	attr->max_seq_len = -1;
	n_seq = 0;
	size  = 0;
	seq = NULL;
	reset_filereader(fr);
	while(attr->is_fq? fread_fastq(&seq, fr) : fread_fasta(&seq, fr)){
		if(seq->seq.size > attr->max_seq_len) attr->max_seq_len = seq->seq.size;
		if(seq->seq.size < attr->min_seq_len) attr->min_seq_len = seq->seq.size;
		size += seq->seq.size;
		n_seq ++;
		if(n_seq > 10000) break;
	}
	if(seq) free_sequence(seq);
	if(n_seq) attr->avg_seq_len = (size + n_seq / 2) / n_seq;
	else attr->avg_seq_len = -1;
	reset_filereader(fr);
}


char *fread_all(FileReader *fr){
	char *text;
	String *line, *string;
	int num;
	line   = init_string(81);
	string = init_string(1023);
	while((num = fread_line2(line, fr)) >= 0){ add_char_string(line, fr->line_breaker); };
	text = string->string;
	free(line->string);
	free(line);
	free(string);
	return text;
}

int reset_filereader(FileReader *fr){
	uint32_t i;
	fr_file_t *fc;
	for(i=0;i<vec_size(fr->files);i++){
		fc = get_vec_ref(fr->files, i);
		fseek(fc->file, 0, SEEK_SET);
	}
	fr->fidx = 0;
	if(vec_size(fr->files)) fr->size = 0;
	fr->ptr  = 0;
	return 1;
}
