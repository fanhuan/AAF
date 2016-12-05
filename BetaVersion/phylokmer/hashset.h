#ifndef __HASH_SET_RJ
#define __HASH_SET_RJ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifndef HASH_FLAG_MACROS
#define HASH_FLAG_MACROS
#define is_entity_null(flags, idx)    ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x01)
#define is_entity_del(flags, idx)     ((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x02)
#define exists_entity(flags, idx)     (!((flags)[(idx)>>4]>>(((idx)&0x0f)<<1)&0x03))
#define set_entity_null(flags, idx)   ((flags)[(idx)>>4] |= (0x01u<<(((idx)&0x0f)<<1)))
#define set_entity_del(flags, idx)    ((flags)[(idx)>>4] |= (0x02u<<(((idx)&0x0f)<<1)))
#define clear_entity_null(flags, idx) ((flags)[(idx)>>4] &= ~(0x01u<<(((idx)&0x0f)<<1)))
#define clear_entity_del(flags, idx)  ((flags)[(idx)>>4] &= ~(0x02u<<(((idx)&0x0f)<<1)))
#endif

#define int64_hashcode(key) (uint32_t)((key)>>33^(key)^(key)<<11)

static inline uint32_t __ac_X31_hash_string(const char *s){
	uint32_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}


#define init_hashset_macro(hash_type, hash_key_type) \
typedef struct { void *array;  uint32_t *flags; size_t e_size; size_t size; size_t count; size_t max; float load_factor; size_t iter_ptr; } hash_type; \
static inline int hash_type##_is_prime(size_t num){                            \
	size_t i, max;                                                             \
	if(num < 4) return 1;                                                      \
	if(num % 2 == 0) return 0;                                                 \
	max = (size_t)sqrt((float)num);                                            \
	for(i=3;i<max;i+=2){ if(num % i == 0) return 0; }                          \
	return 1;                                                                  \
}                                                                              \
static inline size_t hash_type##_find_next_prime(size_t num){                  \
	if(num % 2 == 0) num ++;                                                   \
	while(1){ if(hash_type##_is_prime(num)) return num; num += 2; }            \
}                                                                              \
static inline hash_type* init_##hash_type(){                                   \
	hash_type *set;                                                            \
	set = (hash_type*)malloc(sizeof(hash_type));                               \
	set->e_size = sizeof(hash_key_type);                                       \
	set->size   = 13;                                                          \
	set->count  = 0;                                                           \
	set->load_factor = 0.73f;                                                  \
	set->max    = set->size * set->load_factor;                                \
	set->iter_ptr    = 0;                                                      \
	set->array       = calloc(set->size, set->e_size);                         \
	set->flags       = malloc((set->size + 15)/16 * 4);                        \
	memset(set->flags, 0x55, (set->size + 15) / 16 * 4);                       \
	return set;                                                                \
}

#define get_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline hash_key_type* get_##hash_type(hash_type *set, hash_key_type key){\
	hash_key_type *e;                                                          \
	size_t hc;                                                                 \
	hc = hash_code_macro(key) % set->size;                                     \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			return NULL;                                                       \
		} else if(is_entity_del(set->flags, hc)){                              \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)) return e;                            \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return NULL;                                                               \
}

#define prepare_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline void encap_##hash_type(hash_type *set, size_t num);              \
static inline hash_key_type* prepare_##hash_type(hash_type *set, hash_key_type key, int *exists){\
	hash_key_type *e;                                                          \
	size_t hc, d;                                                              \
	encap_##hash_type(set, 1);                                                 \
	hc = hash_code_macro((key)) % set->size;                                  \
	d = set->size;                                                             \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			if(d == set->size){                                                \
				clear_entity_null(set->flags, hc);                             \
			} else {                                                           \
				hc = d;                                                        \
				clear_entity_del(set->flags, hc);                              \
			}                                                                  \
			*exists = 0;                                                       \
			set->count ++;                                                     \
			e = ((hash_key_type*)set->array) + hc;                             \
			return e;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
			if(d == set->size) d = hc;                                         \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro((*e), (key))){                                \
				*exists = 1;                                                   \
				return e;                                                      \
			}                                                                  \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return NULL;                                                               \
}

#define exists_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline int exists_##hash_type(hash_type *set, hash_key_type key){       \
	hash_key_type *e;                                                          \
	size_t hc;                                                                 \
	hc = hash_code_macro(key) % set->size;                                     \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			return 0;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)) return 1;                            \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return 0;                                                                  \
}

#define add_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline hash_key_type* add_##hash_type(hash_type *set, hash_key_type key){       \
	hash_key_type *e;                                                          \
	size_t d, hc;                                                              \
	hc = hash_code_macro(key) % set->size;                                     \
	d  = set->size;                                                            \
	do{                                                                        \
		if(is_entity_null(set->flags, hc)){                                    \
			if(d == set->size){                                                \
				d = hc;                                                        \
				clear_entity_null(set->flags, d);                              \
			} else {                                                           \
				clear_entity_del(set->flags, d);                               \
			}                                                                  \
			e = ((hash_key_type*)set->array) + d;                              \
			*e = key;                                                          \
			set->count ++;                                                     \
			return e;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
			if(d == set->size) d = hc;                                         \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)){                                     \
				return e;                                                      \
			}                                                                  \
		}                                                                      \
		if(hc + 1 == set->size) hc = 0;                                        \
		else hc = hc + 1;                                                      \
	} while(1);                                                                \
	return NULL;                                                                  \
}

#define put_hashset_macro(hash_type, hash_key_type) \
static inline hash_key_type* put_##hash_type(hash_type *set, hash_key_type key){         \
	encap_##hash_type(set, 1);                                                 \
	return add_##hash_type(set, key);                                          \
}

#define remove_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro) \
static inline int remove_##hash_type(hash_type *set, hash_key_type key){       \
	hash_key_type *e;                                                          \
	size_t hc;                                                                 \
	hc = hash_code_macro(key) % set->size;                                     \
	while(1){                                                                  \
		if(is_entity_null(set->flags, hc)){                                    \
			return 0;                                                          \
		} else if(is_entity_del(set->flags, hc)){                              \
		} else {                                                               \
			e = ((hash_key_type*)set->array) + hc;                             \
			if(hash_equal_macro(*e, key)){                                     \
				set->count --;                                                 \
				set_entity_del(set->flags, hc);                                \
				return 1;                                                      \
			}                                                                  \
		}                                                                      \
		hc ++;                                                                 \
		hc %= set->size;                                                       \
	}                                                                          \
	return 0;                                                                  \
}

#define reset_iter_hashset_macro(hash_type) static inline void reset_iter_##hash_type(hash_type *set){ set->iter_ptr = 0; }

#define iter_hashset_macro(hash_type, hash_key_type) \
static inline int iter_##hash_type(hash_type *set, hash_key_type *ret){        \
	if(set->iter_ptr >= set->size) return 0;                                   \
	while(set->iter_ptr < set->size){                                          \
		if(exists_entity(set->flags, set->iter_ptr)){                          \
			*ret = *(((hash_key_type*)set->array) + set->iter_ptr);            \
			set->iter_ptr ++;                                                  \
			return 1;                                                          \
		}                                                                      \
		set->iter_ptr ++;                                                      \
	}                                                                          \
	return 0;                                                                  \
}

#define ref_iter_hashset_macro(hash_type, hash_key_type) \
static inline int ref_iter_##hash_type(hash_type *set, hash_key_type **ret){        \
	if(set->iter_ptr >= set->size) return 0;                                   \
	while(set->iter_ptr < set->size){                                          \
		if(exists_entity(set->flags, set->iter_ptr)){                          \
			*ret = (((hash_key_type*)set->array) + set->iter_ptr);             \
			set->iter_ptr ++;                                                  \
			return 1;                                                          \
		}                                                                      \
		set->iter_ptr ++;                                                      \
	}                                                                          \
	return 0;                                                                  \
}

#define iter_remove_hashset_macro(hash_type, hash_key_type) \
static inline int iter_remove_##hash_type(hash_type *set, hash_key_type* ret){ \
	if(set->iter_ptr >= set->size) return 0;                                   \
	while(set->iter_ptr < set->size){                                          \
		if(exists_entity(set->flags, set->iter_ptr)){                          \
			set_entity_del(set->flags, set->iter_ptr);                         \
			*ret = *(((hash_key_type*)set->array) + set->iter_ptr);            \
			set->iter_ptr ++;                                                  \
			return 1;                                                          \
		}                                                                      \
		set->iter_ptr ++;                                                      \
	}                                                                          \
	return 0;                                                                  \
}

#define count_hashset_macro(hash_type) static inline int64_t count_##hash_type(hash_type *set){ return set->count; }

#define clear_hashset_macro(hash_type) \
static inline void clear_##hash_type(hash_type *set){                          \
	memset(set->flags, 0x55, (set->size + 15) / 16 * 4);                       \
	set->count = 0;                                                            \
	set->iter_ptr = 0;                                                         \
}

#define ffwrite(ptr, e_size, size, file) (e_size * fwrite(ptr, e_size, size, file))
#define ffread(ptr, e_size, size, file) (e_size * fread(ptr, e_size, size, file))

#define dump_hashset_macro(hash_type) \
static inline size_t sizeof_##hash_type(hash_type *set){                       \
	return sizeof(size_t) * 3 + sizeof(float) + set->e_size * set->size        \
				+ sizeof(uint32_t) * ((set->size + 15) / 16);                  \
}                                                                              \
static inline size_t dump_##hash_type(hash_type *set, FILE *out){              \
	size_t len, i, n;                                                          \
	n =  ffwrite(&set->e_size, sizeof(size_t), 1, out);                        \
	n += ffwrite(&set->size, sizeof(size_t), 1, out);                          \
	n += ffwrite(&set->count, sizeof(size_t), 1, out);                         \
	n += ffwrite(&set->load_factor, sizeof(float), 1, out);                    \
	for(i=0;i<set->size;i+=1024){                                              \
		len = set->size - i;                                                   \
		if(len > 1024) len = 1024;                                             \
		n += ffwrite(set->array + i * set->e_size, set->e_size, len, out);     \
	}                                                                          \
	for(i=0;i<(set->size+15)/16;i+=1024){                                      \
		len = (set->size + 15) / 16 - i;                                       \
		if(len > 1024) len = 1024;                                             \
		n += ffwrite(set->flags + i, sizeof(uint32_t), len, out);              \
	}                                                                          \
	fflush(out);                                                               \
	return n;                                                                  \
}

#define load_hashset_macro(hash_type) \
static inline hash_type* load_##hash_type(FILE *in){                           \
	hash_type *set;                                                            \
	size_t n;                                                                  \
	set = (hash_type*)malloc(sizeof(hash_type));                               \
	n =  ffread(&set->e_size, sizeof(size_t), 1, in);                          \
	n += ffread(&set->size, sizeof(size_t), 1, in);                            \
	n += ffread(&set->count, sizeof(size_t), 1, in);                           \
	n += ffread(&set->load_factor, sizeof(float), 1, in);                      \
	set->max   = set->size * set->load_factor;                                 \
	set->array = malloc(set->size * set->e_size);                              \
	n += ffread(set->array, set->e_size, set->size, in);                       \
	set->flags = (uint32_t*)malloc((set->size + 15) / 16 * 4);                 \
	n += ffread(set->flags, sizeof(uint32_t), (set->size + 15) / 16, in);      \
	return set;                                                                \
}

#define free_hashset_macro(hash_type) \
static inline void free_##hash_type(hash_type *set){                           \
	free(set->array);                                                          \
	free(set->flags);                                                          \
	free(set);                                                                 \
}

#define encap_hashset_macro(hash_type, hash_key_type, hash_code_macro) \
static inline void encap_##hash_type(hash_type *set, size_t num){             \
	uint32_t *flags, *f;                                                      \
	uint64_t i, n, size, hc;                                                  \
	hash_key_type key, tmp;                                                   \
	if(set->count + num <= set->max) return;                                  \
	n = set->size;                                                            \
	do{ n = hash_type##_find_next_prime(n * 2); } while(n * set->load_factor < set->count + num);    \
	set->array = realloc(set->array, n * set->e_size);                        \
	if(set->array == NULL){                                                   \
		fprintf(stderr, "-- Out of memory --\n");                             \
		abort();                                                              \
	}                                                                         \
	flags = malloc((n+15)/16 * 4);                                            \
	memset(flags, 0x55, (n+15)/16 * 4);                                       \
	size = set->size;                                                         \
	set->size = n;                                                            \
	set->max = n * set->load_factor;                                          \
	f = set->flags;                                                           \
	set->flags = flags;                                                       \
	flags = f;                                                                \
	for(i=0;i<size;i++){                                                      \
		if(!exists_entity(flags, i)) continue;                                \
		key = ((hash_key_type*)set->array)[i];                                \
		set_entity_del(flags, i);                                             \
		while(1){                                                             \
			hc = hash_code_macro(key) % set->size;                            \
			while(!is_entity_null(set->flags, hc)){ hc = (hc + 1) % set->size; }        \
			clear_entity_null(set->flags, hc);                                \
			if(hc < size && exists_entity(flags, hc)){                        \
				tmp = key;                                                    \
				key = ((hash_key_type*)set->array)[hc];                       \
				((hash_key_type*)set->array)[hc] = tmp;                       \
				set_entity_del(flags, hc);                                    \
			} else {                                                          \
				((hash_key_type*)set->array)[hc] = key;                       \
				break;                                                        \
			}                                                                 \
		}                                                                     \
	}                                                                         \
	free(flags);                                                              \
}



// ---------------------- Define your own hashset ----------------------------------
// Example: 
// typedef struct { int group; int user; } Info;
// #define my_hashcode(val) (val)->group
// #define my_hashequal(v1, v2) (((v1)->group == (v2)->group) && ((v1)->user == (v2)->user))
// define_hashset(myhash, Info, my_hashcode, my_hashequal);

#define define_hashset(hash_type, hash_key_type, hash_code_macro, hash_equal_macro)    \
	init_hashset_macro(hash_type, hash_key_type);                              \
	get_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro);    \
	prepare_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro);\
	exists_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro); \
	add_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro);    \
	put_hashset_macro(hash_type, hash_key_type);                               \
	remove_hashset_macro(hash_type, hash_key_type, hash_code_macro, hash_equal_macro); \
	iter_hashset_macro(hash_type, hash_key_type);                              \
	ref_iter_hashset_macro(hash_type, hash_key_type);                          \
	iter_remove_hashset_macro(hash_type, hash_key_type);                       \
	reset_iter_hashset_macro(hash_type);                                       \
	count_hashset_macro(hash_type);                                            \
	clear_hashset_macro(hash_type);                                            \
	dump_hashset_macro(hash_type);                                             \
	load_hashset_macro(hash_type);                                             \
	free_hashset_macro(hash_type);                                             \
	encap_hashset_macro(hash_type, hash_key_type, hash_code_macro);

/* ------------------ Useful functions ------------------------------------- */
static inline uint32_t simp_string_hashcode(const char *s){
	uint32_t  h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}

#endif
