#include "file_reader.h"
#include "hashset.h"
#include <unistd.h>
#include <string.h>

typedef struct {
	char *str;
	int val;
} StrInt;

#define string_hashcode(str1) simp_string_hashcode((str1).str)
#define string_equals(str1, str2) (strcmp((str1).str, (str2).str) == 0)
define_hashset(strhash, StrInt, string_hashcode, string_equals);

typedef struct {
	FileReader *in;
	int has_next;
	int req_next;
	int is_select;
	int filter_flag;
	int print_flag;
	int n_col;
	char *default_col;
	int skip_header;
	int repeated;
	Vector *cmp_cols;
} DataFile;

static strhash *enums = NULL;

static uint32_t enum_val = 1;

int cmp_fields(DataFile *file1, DataFile *file2, Vector *col_types){
	StrInt *h, H;
	uint32_t i, type;
	char *s1, *s2;
	double n1, n2;
	int cmp;
	for(i=0;i<vec_size(col_types);i++){
		type = gget_vec(col_types, i, uint32_t);
		s1 = ((String*)get_vec_ref(file1->in->tabs, gget_vec(file1->cmp_cols, i, uint32_t)))->string;
		s2 = ((String*)get_vec_ref(file2->in->tabs, gget_vec(file2->cmp_cols, i, uint32_t)))->string;
		if(type & 0x01){
			if(type & 0x08){
				H.str = s1;
				h = get_strhash(enums, H);
				if(h == NULL) n1 = 0;
				else n1 = h->val;
				H.str = s2;
				h = get_strhash(enums, H);
				if(h == NULL) n2 = 0;
				else n2 = h->val;
			} else {
				n1 = atof(s1);
				n2 = atof(s2);
			}
			cmp = (n1 == n2)? 0 : ((n1 < n2)? -1 : 1);
		} else {
			cmp = (type&0x04)? strcasecmp(s1, s2) : strcmp(s1, s2);
		}
		if(cmp == 0) continue;
		else return cmp * ((type&0x02)? -1 : 1);
	}
	return 0;
}

int print_filter(Vector *files, FILE *out){
	uint32_t i, k, s, m, n, m2, n2;
	int j;
	DataFile *file, *f;
	k = 0;
	m = n = 0;
	m2 = n2 = 0;
	f = NULL;
	for(i=0;i<vec_size(files);i++){
		file = gget_vec(files, i, DataFile*);
		if(file->is_select){ k ++; f = file; }
		switch(file->filter_flag){
			case 0: if(file->is_select) return 0; break;
			case 1: if(!file->is_select) return 0; break;
			case 2: m ++; if(file->is_select) n++; break;
			case 4: m2 ++; if(file->is_select) n2++; break;
			case 8: break;
		}
	}
	if(m && n == 0) return 0;
	if(m2 && n2 < 2) return 0;
	n = 0;
	m = 0;
	for(i=0;i<vec_size(files);i++){
		file = gget_vec(files, i, DataFile*);
		if(file->repeated) file->req_next = 0;
		if(file->print_flag&0x01) continue;
		if(file->is_select && file->print_flag&0x02) m ++;
		if(file->print_flag&0x04){
			if(file->n_col > 0){
				if(n) fprintf(out, "\t");
				n ++;
				if(file->is_select){
					fprintf(out, "%s", file->in->line->string);
				} else {
					if(file->print_flag & 0x08){
						k = 0;
						for(j=0;j<file->n_col;j++){
							if(j) fprintf(out, "\t");
							s = 0xFFFFFFFFU;
							for(k=0;f&&k<vec_size(file->cmp_cols);k++){
								if(j == (int)gget_vec(file->cmp_cols, k,uint32_t)){ s = k; break; }
							}
							if(s != 0xFFFFFFFFU) fprintf(out, "%s", ((String*)get_vec_ref(f->in->tabs, gget_vec(f->cmp_cols, s, uint32_t)))->string);
							else fprintf(out, "%s", file->default_col);
						}
					} else {
						for(j=0;j<file->n_col;j++){
							if(j) fprintf(out, "\t");
							fprintf(out, "%s", file->default_col);
						}
					}
				}
			} else {
				if(!file->is_select) continue;
				if(n) fprintf(out, "\t");
				n ++;
				fprintf(out, "[%c]\t", i + 'A');
				fprintf(out, "%s", file->in->line->string);
			}
			continue;
		}
		if(!file->is_select){ continue; }
		if((file->print_flag&0x02) && m) continue;
		if(n) fprintf(out, "\t");
		n ++;
		fprintf(out, "%s", file->in->line->string);
	}
	fprintf(out, "\n");
	return 1;
}

int run_filter(Vector *files, Vector *col_types, FILE *out){
	uint32_t i;
	int cmp, ret, n_col;
	DataFile *file, *f;
	Vector *hits;
	for(i=0;i<vec_size(files);i++){
		file = gget_vec(files, i, DataFile*);
		file->req_next = 1;
		if(file->skip_header){
			while(fread_table(file->in) != -1){
				if(get_col_str(file->in, 0)[0] != '#'){
					file->req_next = 0;
					break;
				}
			}
		}
		file->has_next = 1;
	}
	ret = 0;
	hits = init_vec(sizeof(uint32_t), 6);
	while(1){
		clear_vec(hits);
		for(i=0;i<vec_size(files);i++){
			file = gget_vec(files, i, DataFile*);
			file->is_select = 0;
			if(file->req_next){
				if(file->has_next){
					if((n_col = fread_table(file->in)) == -1){
						file->has_next = 0;
						continue;
					} else if(file->n_col == 0) file->n_col = n_col;
				} else continue;
				file->req_next = 0;
			}
			if(vec_size(hits)){
				f = gget_vec(files, gget_vec(hits, 0, uint32_t), DataFile*);
				cmp = cmp_fields(file, f, col_types);
				if(cmp == 0) gpush_vec(hits, i, uint32_t);
				else if(cmp < 0) clear_vec(hits), gpush_vec(hits, i, uint32_t);
			} else {
				gpush_vec(hits, i, uint32_t);
			}
		}
		if(vec_size(hits) == 0) break;
		for(i=0;i<vec_size(hits);i++){
			file = gget_vec(files, gget_vec(hits, i, uint32_t), DataFile*);
			file->req_next = 1;
			file->is_select = 1;
		}
		ret += print_filter(files, out);
	}
	free_vec(hits);
	return ret;
}

int usage(char *prog){
	printf(
"FILTER 1.0 -- comparing sorted files, written by Jue Ruan <ruanjue@gmail.com>\n"
"Usage: %s [options] file1 file2 ...\n"
"Options:\n"
" -k <key_field_types>  A string of 's', 'S', 'c', 'C', 'n', 'N', 'e' and 'E'\n"
"                       's': string;\n"
"                       'c': case-insensitive string;\n"
"                       'n': number;\n"
"                       'e': Enum;\n"
"                        Capital letters: descent.\n"
" -e <enums>            A list of string to define a enum, demiliter by ','\n"
"                       'chr1,chr2,chr3,chrX' indicates chr1 < chr2 < chr3 < chrX\n"
" -b                    Built-in enum\n"
"                       including chromosomes(chr1, chr2, ..., chrX, chrY, chrM)\n"
" -c                    The number of columns is fixed in each file\n"
" -d <str>              Default column string, default: NULL\n"
//" -r                    Unsort mode, just compare line by line\n"
" -A <expr>             Tell how to process file1\n"
"                       <expr> = <flags>,<fields>\n"
"                       <flags> = [RNAOTMPpS], splitted by ','. Default: 'A'\n"
"                       'R': if records in other files are repeated in multple lines,\n"
"                            but unique in this file, please use 'R' to match all of them\n"
"                       'N': cannot exist in this file\n"
"                       'A': must exist in this file\n"
"                       'E': this file is free of existence\n"
"                       'O': must exist at least once in all files with flag 'O'\n"
"                       'T': must exist at least twice in all files with flag 'T'\n"
"                       'M': mark the content from which file in output,\n"
"                            if fix_column set, output fixed columns, \n"
"                            otherwise use  [A-Z] mark each item\n"
"                       'p': don't output content from this file\n"
"                       'P': output this when no other file content\n"
"                       'F': auto fill in key columns if not exists in this file\n"
"                       'S': skip header line\n"
"                       <field_list> = key1,key2,...\n"
"                       col_ids responding to <key_field_types>\n"
"-B <expr>              Tell how to process file2\n"
"                        and so on\n"
"-a <expr>              Tell how to process all files\n"
"\n", prog
);
	return 1;
}

#define MAX_FILES	2000

int main(int argc, char **argv){
	FileReader *in;
	Vector *files, *col_types, *tokens;
	DataFile *file;
	StrInt H;
	VirtualString *str;
	String *astr;
	Vector *strs;
	int c, sorted;
	uint32_t i, j, type, last;
	char name[128];
	sorted = 1;
	files = init_vec(sizeof(DataFile*), MAX_FILES);
	for(i=0;i<MAX_FILES;i++){
		file = malloc(sizeof(DataFile));
		gpush_vec(files, file, DataFile*);
		file->in  = NULL;
		file->has_next  = 1;
		file->req_next  = 1;
		file->is_select = 0;
		file->filter_flag = 1;
		file->print_flag  = 0;
		file->n_col = -1;
		file->default_col = "NULL";
		file->skip_header = 0;
		file->repeated    = 0;
		file->cmp_cols  = init_vec(sizeof(uint32_t), 6);
	}
	col_types = init_vec(sizeof(uint32_t), 6);
	tokens = init_vec(sizeof(VirtualString), 6);
	enum_val = 0;
	enums = init_strhash(13);
	strs = init_vec(sizeof(String*), 6);
	while((c = getopt(argc, argv, "a:bcrd:e:k:A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z:")) != -1){
		if(c == 'c'){
			for(i=0;i<MAX_FILES;i++){
				file = gget_vec(files, i, DataFile*);
				file->n_col = 0;
			}
		} else if(c == 'r'){
			sorted = 0;
		} else if(c == 'd'){
			for(i=0;i<MAX_FILES;i++){
				file = gget_vec(files, i, DataFile*);
				file->default_col = optarg;
			}
		} else if(c == 'k'){
			for(i=0;i<(uint32_t)strlen(optarg);i++){
				type = 0;
				switch(optarg[i]){
					case 's': type = 0; break;
					case 'S': type = 0 + 2; break;
					case 'c': type = 4; break;
					case 'C': type = 4 + 2; break;
					case 'n': type = 1; break;
					case 'N': type = 1 + 2; break;
					case 'e': type = 1 + 8; break;
					case 'E': type = 1 + 8  + 2; break;
					default: fprintf(stderr, "Illegal character '%c' in %s[%d]\n", optarg[i], optarg, i);
				}
				gpush_vec(col_types, type, uint32_t);
			}
		} else if(c == 'e'){
			str = as_vstring(optarg);
			clear_vec(tokens);
			split_vstring(str, ',', tokens, 1);
			free_vstring(str);
			for(i=0;i<vec_size(tokens);i++){
				str = get_vec_ref(tokens, i);
				trim_vstring(str);
				if(str->size == 0) continue;
				gpush_vec(strs, as_string(str->string), String*);
				astr = gget_vec(strs, vec_size(strs) - 1, String*);
				H.str = astr->string;
				H.val = enum_val ++;
				put_strhash(enums, H);
			}
		} else if(c == 'b'){
			for(i=0;i<80;i++){
				sprintf(name, "chr%u", i + 1);
				gpush_vec(strs, as_string(name), String*);
				astr = gget_vec(strs, vec_size(strs) - 1, String*);
				H.str = astr->string;
				H.val = enum_val ++;
				put_strhash(enums, H);
			}
			{
				gpush_vec(strs, as_string("chrX"), String*);
				astr = gget_vec(strs, vec_size(strs) - 1, String*);
				H.str = astr->string;
				H.val = enum_val ++;
				put_strhash(enums, H);
			}
			{
				gpush_vec(strs, as_string("chrY"), String*);
				astr = gget_vec(strs, vec_size(strs) - 1, String*);
				H.str = astr->string;
				H.val = enum_val ++;
				put_strhash(enums, H);
			}
			{
				gpush_vec(strs, as_string("chrM"), String*);
				astr = gget_vec(strs, vec_size(strs) - 1, String*);
				H.str = astr->string;
				H.val = enum_val ++;
				put_strhash(enums, H);
			}
		} else if(c == 'a'){
			str = as_vstring(optarg);
			clear_vec(tokens);
			split_vstring(str, ',', tokens, 1);
			free_vstring(str);
			for(j=0;j<MAX_FILES;j++){
				file = gget_vec(files, j, DataFile*);
				for(i=0;i<vec_size(tokens);i++){
					str = get_vec_ref(tokens, i);
					trim_vstring(str);
					if(str->size == 0) continue;
					switch(str->string[0]){
						case 'R': fprintf(stderr, "Cannot assign 'R' to all files\n"); return usage(argv[0]);
						case 'N': file->filter_flag = 0; break;
						case 'A': file->filter_flag = 1; break;
						case 'E': file->filter_flag = 8; break;
						case 'O': file->filter_flag = 2; break;
						case 'T': file->filter_flag = 4; break;
						case 'M': file->print_flag |= 4; break;
						case 'P': file->print_flag |= 1; break;
						case 'p': file->print_flag |= 2; break;
						case 'S': file->skip_header = 1; break;
						case 'F': file->print_flag |= 8; break;
						default: if(str->string[0] < '0' || str->string[0] > '9'){
									fprintf(stderr, "Cannot parse %s in -%c %s\n", str->string, (char)c, optarg);
									return usage(argv[0]);
								}
								gpush_vec(file->cmp_cols, atoi(str->string)-1, uint32_t);
					}
				}
			}
		} else if(c >= 'A' && c <= 'Z'){
			file = gget_vec(files, c - 'A', DataFile*);
			str = as_vstring(optarg);
			clear_vec(tokens);
			split_vstring(str, ',', tokens, 1);
			free_vstring(str);
			for(i=0;i<vec_size(tokens);i++){
				str = get_vec_ref(tokens, i);
				trim_vstring(str);
				if(str->size == 0) continue;
				switch(str->string[0]){
					case 'R': file->repeated    = 1; break;
					case 'N': file->filter_flag = 0; break;
					case 'A': file->filter_flag = 1; break;
					case 'E': file->filter_flag = 8; break;
					case 'O': file->filter_flag = 2; break;
					case 'T': file->filter_flag = 4; break;
					case 'M': file->print_flag |= 4; break;
					case 'P': file->print_flag |= 1; break;
					case 'p': file->print_flag |= 2; break;
					case 'F': file->print_flag |= 8; break;
					case 'S': file->skip_header = 1; break;
					default: if(str->string[0] < '0' || str->string[0] > '9'){
								 fprintf(stderr, "Cannot parse %s in -%c %s\n", str->string, (char)c, optarg);
								 return usage(argv[0]);
							 }
							 gpush_vec(file->cmp_cols, atoi(str->string)-1, uint32_t);
				}
			}
		}
	}
	free_vec(tokens);
	if(optind == argc) return usage(argv[0]);
	for(i=optind;i<(uint32_t)argc;i++){
		if(i-optind >= MAX_FILES){
			fprintf(stderr, "Too many files(max = %d), ignore the remains\n", MAX_FILES);
			break;
		}
		if((in = fopen_filereader(argv[i])) == NULL){
			fprintf(stderr, "Cannot open %s\n", argv[i]);
			return 1;
		}
		gget_vec(files, i - optind, DataFile*)->in = in;
	}
	set_vec_size(files, i - optind);
	for(i=0;i<vec_size(files);i++){
		file = gget_vec(files, i, DataFile*);
		if(vec_size(file->cmp_cols) >= vec_size(col_types)) continue;
		last = vec_size(file->cmp_cols)? gget_vec(file->cmp_cols, vec_size(file->cmp_cols) - 1, uint32_t) + 1 : 0;
		for(j=vec_size(file->cmp_cols);j<vec_size(col_types);j++){
			gpush_vec(file->cmp_cols, last, uint32_t);
			last ++;
		}
	}
	run_filter(files, col_types, stdout);
	free_vec(col_types);
	for(i=0;i<MAX_FILES;i++){
		file = gget_vec(files, i, DataFile*);
		if(file->in) fclose_filereader(file->in);
		free_vec(file->cmp_cols);
		free(file);
	}
	for(i=0;i<vec_size(strs);i++) free_string(gget_vec(strs, i, String*));
	free_vec(strs);
	free_vec(files);
	free_strhash(enums);
	return 0;
}
