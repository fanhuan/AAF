#include "hashset.h"
#include "file_reader.h"
#include <time.h>
#include <unistd.h>

typedef struct {
	uint64_t kmer:50, count:14;
} kmer_t;

#define MAX_COUNT 0x3FFFU
#define MIN_KSIZE	3
#define MAX_KSIZE	25

static inline uint64_t __lh3_Jenkins_hash_64(uint64_t key){
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}

#define kmer_hashcode(k1) __lh3_Jenkins_hash_64((k1).kmer)
#define kmer_equals(k1, k2) ((k1).kmer == (k2).kmer)
define_hashset(kmerhash, kmer_t, kmer_hashcode, kmer_equals);

static const uint8_t base_bit_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static int dump_id = 0;
static char prefix[120];
static uint64_t mem_limit = 2LLU * 1024 * 1024 * 1024;

typedef struct {
	uint16_t symbol:3, count:13;
	uint8_t  links[4];
} TrieNode;

uint32_t cal_kmer_complexity(uint64_t seq, int seq_size, Vector *nodes){
	int i, j;
	uint8_t v;
	TrieNode *node, *root;
	clear_vec(nodes);
	root = get_next_vec_ref(nodes);
	root->symbol = 4;
	root->count = 0;
	root->links[0] = 0;
	root->links[1] = 0;
	root->links[2] = 0;
	root->links[3] = 0;
	for(i=0;i<seq_size;i++){
		node = root;
		for(j=i;j<seq_size;j++){
			v = (seq >> ((seq_size-j-1)<<1)) & 0x03;
			if(node->links[v]) node = get_vec_ref(nodes, node->links[v]);
			else {
				node->links[v] = vec_size(nodes);
				node = get_next_vec_ref(nodes);
				node->symbol = v;
				node->count = 0;
				node->links[0] = node->links[1] = node->links[2] = node->links[3] = 0;
			}
			node->count ++;
		}
	}
	return vec_size(nodes) - 1;
}

static inline uint64_t get_rev_seq(uint64_t seq, int seq_size){
	seq = ~seq;
	seq = ((seq & 0x3333333333333333LLU)<< 2) | ((seq & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
	seq = ((seq & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq & 0xFF00FF00FF00FF00LLU)>> 8);
	seq = ((seq & 0x0000FFFF0000FFFFLLU)<<16) | ((seq & 0xFFFF0000FFFF0000LLU)>>16);
	seq = ((seq & 0x00000000FFFFFFFFLLU)<<32) | ((seq & 0xFFFFFFFF00000000LLU)>>32);
	return seq >> (64 - (seq_size<<1));
}

void dump_kmers(kmerhash *hash){
	FILE *file;
	char fname[128];
	kmer_t *mer;
	size_t t;
	t = 0;
	sprintf(fname, "%s.%04d", prefix, dump_id);
	if((file = fopen(fname, "w")) == NULL){
		fprintf(stderr, " -- Cannot write %s in %s -- %s:%d --\n", fname, __FUNCTION__, __FILE__, __LINE__);
		abort();
	}
	reset_iter_kmerhash(hash);
	while(ref_iter_kmerhash(hash, &mer)){
		t += fwrite(mer, sizeof(kmer_t), 1, file);
	}
	fclose(file);
	clear_kmerhash(hash);
	dump_id ++;
}

void counting(kmerhash *hash, int kmer_size, uint64_t kmer_mask, Sequence *seq){
	int i, last;
	uint8_t b;
	uint64_t v, vv;
	int exist;
	kmer_t MER, *mer;
	if(seq->seq.size < kmer_size) return;
	v = 0;
	last = -1;
	if((2LLU * hash->size * hash->e_size) >= mem_limit && (uint64_t)(hash->count + seq->seq.size - kmer_size) >= hash->max) dump_kmers(hash);
	MER.count = 0;
	for(i=0;i<seq->seq.size;i++){
		b = base_bit_table[(int)seq->seq.string[i]];
		if(b > 3) last = i;
		v = (v << 2) | b;
		if(i - last < kmer_size) continue;
		v  = v & kmer_mask;
		vv = get_rev_seq(v, kmer_size);
		if(vv > v) vv = v;
		MER.kmer = vv;
		mer = prepare_kmerhash(hash, MER, &exist);
		if(exist){
			if(mer->count < MAX_COUNT) mer->count ++;
		} else {
			mer->kmer = MER.kmer;
			mer->count = 1;
		}
	}
}

int special_filter(uint64_t kmer, int kmer_size){
	int i, nt_cnts[4], cont;
	uint8_t last, cur;
	nt_cnts[0] = nt_cnts[1] = nt_cnts[2] = nt_cnts[3] = cont = 0;
	last = 4;
	for(i=0;i<kmer_size;i++){
		cur = (kmer >> (i << 1)) & 0x03;
		if(cur == last){ cont ++; if(cont > 4) return 0; }
		else { last = cur; cont = 1; }
		nt_cnts[cur] ++;
	}
	for(i=0;i<4;i++){
		if(nt_cnts[i] < 2) return 0;
		if(nt_cnts[i] > 0.5 * kmer_size) return 0;
	}
	return 1;
}

void count_files(kmerhash *hash, int kmer_size, Vector *files, int is_fa, uint32_t n_reads, char *rds_file){
	FileReader *fr;
	Vector *bits;
	Sequence *seq;
	FILE *out;
	uint64_t kmer_mask;
	uint32_t i, nr, total, p;
	kmer_mask = (1LLU << (2 * kmer_size)) - 1LLU;
	bits = NULL;
	srand(time(0));
	out = NULL;
	if(rds_file){
		if((out = fopen(rds_file, "w")) == NULL){
			fprintf(stderr, " -- Cannot open %s for write in %s -- %s:%d --\n", rds_file, __FUNCTION__, __FILE__, __LINE__);
		}
	}
	if(n_reads){
		total = 0;
		for(i=0;i<vec_size(files);i++){
			if((fr = fopen_filereader(gget_vec(files, i, char*))) == NULL){
				fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", gget_vec(files, i, char*), __FUNCTION__, __FILE__, __LINE__);
				abort();
			}
			seq = NULL;
			while((is_fa? fread_fasta_adv(&seq, fr, FASTA_FLAG_NO_NAME) : fread_fastq_adv(&seq, fr, FASTQ_FLAG_NO_NAME | FASTQ_FLAG_NO_QUAL))){
				total ++;
			}
			fclose_filereader(fr);
		}
		if(total == 0) return;
		bits = init_vec(sizeof(uint8_t), (total + 7) / 8);
		if(n_reads >= total) memset(bits->buffer, 0xFF, (total + 7) / 8);
		else {
			p = RAND_MAX * (uint64_t)n_reads / total;
			nr = 0;
			memset(bits->buffer, 0x00, (total + 7) / 8);
			for(i=0;i<total;i++){
				if(rand() < (int)p){
					nr ++;
					((uint8_t*)bits->buffer)[i>>3] |= 1U << (i & 0x07);
					if(nr >= n_reads) break;
				}
			}
			for(i=0;i<total;i++){
				if((((uint8_t*)bits->buffer)[i>>3] >> (i & 0x07)) & 0x01) continue;
				if(nr >= n_reads) break;
				nr ++;
				((uint8_t*)bits->buffer)[i>>3] |= 1U << (i & 0x07);
			}
		}
	}
	nr = 0;
	for(i=0;i<vec_size(files);i++){
		if((fr = fopen_filereader(gget_vec(files, i, char*))) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", gget_vec(files, i, char*), __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		seq = NULL;
		while((is_fa? fread_fasta(&seq, fr) : fread_fastq_adv(&seq, fr, FASTQ_FLAG_NO_QUAL))){
			if(n_reads && !((((uint8_t*)bits->buffer)[nr>>3] >> (nr & 0x07)) & 0x01)){ nr++; continue; }
			//printf(">%d\n%s\n", nr, seq->seq.string);
			nr ++;
			if(out) fprintf(out, ">%s\n%s\n", seq->name.string, seq->seq.string);
			counting(hash, kmer_size, kmer_mask, seq);
		}
		fclose_filereader(fr);
	}
	if(n_reads) free_vec(bits);
	if(out) fclose(out);
}

int cmp_kmer_t(const void *e1, const void *e2){
	kmer_t *k1, *k2;
	k1 = (kmer_t*)e1;
	k2 = (kmer_t*)e2;
	if(k1->kmer == k2->kmer) return 0;
	else if(k1->kmer < k2->kmer) return -1;
	else return 1;
}

uint64_t merge_print_kmers(int kmer_size, uint32_t cutoff, uint64_t n_kmers, int is_special, int is_c, FILE *out){
	FILE *in, **ins;
	int i;
	uint32_t score;
	long long n, t, len;
	uint64_t min, count, nk;
	char fname[128], seqs[32];
	char *flags;
	char *acgt = "ACGT";
	kmer_t *kmers;
	Vector *nodes;
	//fprintf(out, "Sequence\tOccurrence\tComplexity\n");
	if(dump_id < 1) return 0;
	nodes = init_vec(sizeof(TrieNode), 256);
	ins = (FILE**)malloc(sizeof(FILE*) * dump_id);
	for(i=0;i<dump_id;i++){
		sprintf(fname, "%s.%04d", prefix, i);
		if((in = fopen(fname, "r")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", fname, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		fseek(in, 0, SEEK_END);
		len = ftell(in);
		kmers = (kmer_t*)malloc(len);
		fseek(in, 0, SEEK_SET);
		t = 0;
		while((n = fread(kmers + t, sizeof(kmer_t), 1024, in)) > 0){ t += n; }
		fclose(in);
		qsort(kmers, t, sizeof(kmer_t), cmp_kmer_t);
		if((in = fopen(fname, "w")) == NULL){
			fprintf(stderr, " -- Cannot write %s in %s -- %s:%d --\n", fname, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		n = fwrite(kmers, sizeof(kmer_t), t, in);
		fclose(in);
		free(kmers);
		if((ins[i] = fopen(fname, "r")) == NULL){
			fprintf(stderr, " -- Cannot open %s in %s -- %s:%d --\n", fname, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
	}
	kmers = (kmer_t*)malloc(sizeof(kmer_t) * dump_id);
	memset(kmers, 0, sizeof(kmer_t) * dump_id);
	count = 0;
	seqs[kmer_size] = 0;
	nk = 0;
	flags = malloc(dump_id);
	memset(flags, 0, dump_id);
	while(nk < n_kmers){
		min = 0xFFFFFFFFFFFFFFFFLLU;
		count = 0;
		for(i=0;i<dump_id;i++){
			if(flags[i] == 0){ n = fread(kmers + i, sizeof(kmer_t), 1, ins[i]); if(n != 1){ flags[i] = 1; kmers[i].kmer = 0x3FFFFFFFFFFFFLLU; } else flags[i] = 2;  }
			if(flags[i] == 2){
				if(kmers[i].kmer == min) count += kmers[i].count;
				else if(kmers[i].kmer < min) min = kmers[i].kmer, count = kmers[i].count;
			}
		}
		if(min == 0xFFFFFFFFFFFFFFFFLLU) break;
		for(i=0;i<dump_id;i++) if(kmers[i].kmer == min) flags[i] = 0;
		if(count < cutoff) continue;
		if(is_special == 1  && !special_filter(min, kmer_size)) continue;
		if(is_special == -1 && special_filter(min, kmer_size)) continue;
		for(i=0;i<kmer_size;i++){ seqs[i] = acgt[(min >> ((kmer_size - 1 - i) << 1)) & 0x03]; }
		if(is_c){
			score = cal_kmer_complexity(min, (int)kmer_size, nodes);
			fprintf(out, "%s\t%u\t%u\n", seqs, (unsigned)count, score);
		} else {
			fprintf(out, "%s\t%u\n", seqs, (unsigned)count);
		}
		nk ++;
	}
	fflush(out);
	free(flags);
	for(i=0;i<dump_id;i++){
		fclose(ins[i]);
		sprintf(fname, "%s.%04d", prefix, i);
		remove(fname);
	}
	free(ins);
	free(kmers);
	free_vec(nodes);
	return nk;
}

int usage(){
	printf(
"kmer_count 1.0\n"
"Author: Jue Ruan <ruanjue@gmail.com>\n"
"Usage: [options]\n"
"-l <kmer_size>    Size of kmer (3~25), default: 25 (bp)\n"
"-i <input_files>  Reads Files, `-i` can be used multi-times\n"
"-f <file_format>  File format of input files, FA|FQ, default: FA\n"
"-s                Special filter:\n"
"                    1, all nucleotides present at least two times\n"
"                    2, no nucleotide >50%% of fragment\n"
"                    3, no runs of a single nucleotide longer than 4 bases\n"
"-S                Reverse special filter\n"
"-o <output>       Output file\n"
"-C                Calculate sequence complexity\n"
"-N <num_reads>    Randomly subsampling <num_reads> reads\n"
"-k <num_kmers>    Stop load new kmer after threshold\n"
"-r <sample_file>  Dump used reads into sample_file\n"
"-n <freq_cutoff>  Output kmers with freq >= cutoff, default: 0\n"
"-G <max_memory>   Memory limit(Gb), default: 2(Gb)\n"
"-M <max_memory>   Memory limit(Mb), default: 2 * 1024(Mb)\n"
"-K <max_memory>   Memory limit(Kb), default: 2 * 1024 * 1024 (kb)\n"
"\n"
"Example:\n"
" kmer_count -l 25 -f FQ -o freq.txt -r reads.fa -i file1.fq -i file2.fq -i file3.fq -M 256\n"
"\n"
);
	return 1;
}

int main(int argc, char **argv){
	kmerhash *hash;
	Vector *files;
	FILE *out;
	char *out_file, *rds_file;
	int kmer_size;
	uint64_t n_kmers, ret;
	uint32_t cutoff, n_reads;
	int is_fa, is_special, is_c, gz_out, c;
	char *cmd;
	kmer_size = 25;
	is_fa = 1;
	is_c  = 0;
	is_special = 0;
	files = init_vec(sizeof(char*), 6);
	out_file = NULL;
	rds_file = NULL;
	cutoff = 0;
	n_kmers = 0xFFFFFFFFFFFFFFFFLLU;
	n_reads = 0;
	sprintf(prefix, "kmer_count.tmp_file.pid%d", getpid());
	while((c = getopt(argc, argv, "CsSl:n:N:k:i:f:o:r:G:M:K:")) != -1){
		switch(c){
			case 'l': kmer_size = atoi(optarg); break;
			case 'f': is_fa = (strcasecmp(optarg, "FQ") != 0); break;
			case 'o': out_file = optarg; break;
			case 'i': gpush_vec(files, optarg, char*); break;
			case 'n': cutoff = atoi(optarg); break;
			case 'N': n_reads = atoi(optarg); break;
			case 'k': n_kmers = atoll(optarg); break;
			case 'r': rds_file = optarg; break;
			case 's': is_special = 1; break;
			case 'S': is_special = -1; break;
			case 'C': is_c       = 1; break;
			case 'G': mem_limit = atoi(optarg) * 1024LLU * 1024 * 1024; break;
			case 'M': mem_limit = atoi(optarg) * 1024LLU * 1024; break;
			case 'K': mem_limit = atoi(optarg) * 1024LLU; break;
		}
	}
	if(vec_size(files) == 0) return usage();
	if(out_file == NULL) return usage();
	if(kmer_size > MAX_KSIZE) return usage();
	if(kmer_size < MIN_KSIZE) return usage();
	if(strlen(out_file) > 3 && strcmp(out_file + strlen(out_file) - 3, ".gz") == 0){
		cmd = malloc(strlen(out_file) + 60);
		sprintf(cmd, "gzip -c >%s", out_file);
		if((out = popen(cmd, "w")) == NULL){
			fprintf(stderr, " -- Cannot invoke '%s' in %s -- %s:%d --\n", cmd, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		gz_out = 1;
	} else {
		if((out = fopen(out_file, "w")) == NULL){
			fprintf(stderr, " -- Cannot write %s in %s -- %s:%d --\n", out_file, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		gz_out = 0;
	}
	hash = init_kmerhash(13);
	count_files(hash, kmer_size, files, is_fa, n_reads, rds_file);
	free_vec(files);
	dump_kmers(hash);
	free_kmerhash(hash);
	ret = merge_print_kmers(kmer_size, cutoff, n_kmers, is_special, is_c, out);
	fprintf(stdout, "Total: %llu kmers\n", (unsigned long long)ret);
	if(gz_out) pclose(out);
	else fclose(out);
	return 0;
}

