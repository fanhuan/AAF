#ifndef __VECTOR_H_RJ
#define __VECTOR_H_RJ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct Vector {
	void *buffer;
	size_t size;
	size_t cap;
	unsigned int e_size;
} Vector;

typedef int (*cmp_vec_fun)(const void *k1, const void *k2);

static inline void init_memvec(Vector *vec, unsigned int e_size, unsigned int init_size){
	vec->e_size = e_size;
	vec->size   = 0;
	vec->cap    = init_size;
	vec->buffer = malloc(((unsigned int)init_size) * e_size);
	if(vec->buffer == NULL){
		fprintf(stderr, " -- Out of memory in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		abort();
	}
	memset(vec->buffer, 0, ((unsigned int)init_size) * e_size);
}

static inline Vector* init_vec(unsigned int e_size, unsigned int init_size){
	Vector *vec = (Vector*)malloc(sizeof(Vector));
	if(vec == NULL){
		fprintf(stderr, " -- Out of memory in %s -- %s:%d --\n", __FUNCTION__, __FILE__, __LINE__);
		abort();
	}
	init_memvec(vec, e_size, init_size);
	return vec;
}

// Faster than gcc's memcpy
static inline void vec_memcpy(void *dst, void *src, size_t size){
	register size_t i;
	for(i=0;i<size;i++){ ((uint8_t*)dst)[i] = ((uint8_t*)src)[i]; }
}

#define vec_size(v) ((v)->size)

static inline int encap_vec(Vector *vec, unsigned int add_size){
	size_t size;
	if(add_size + vec->size > vec->cap){
		size = add_size + vec->size;
		while(size > vec->cap){
			if(vec->cap < 0xFFFFFU){
				if(vec->cap) vec->cap <<= 1;
				else vec->cap = 8;
			} else {
				vec->cap += 0xFFFFFU;
			}
		}
		vec->buffer = realloc(vec->buffer, vec->cap * vec->e_size);
		if(vec->buffer == NULL){
			fprintf(stderr, " -- Out of memory, try alloc %ld bytes in %s -- %s:%d --\n", (long int)vec->cap * vec->e_size, __FUNCTION__, __FILE__, __LINE__);
			abort();
		}
		memset(vec->buffer + vec_size(vec) * vec->e_size, 0, (vec->cap - vec_size(vec)) * vec->e_size);
	}
	return 1;
}

static inline int add_vec_size(Vector *vec, size_t add_size){
	encap_vec(vec, add_size);
	vec->size += add_size;
	return 1;
}

static inline int reduce_vec_size(Vector *vec, size_t size){
	if(size > vec_size(vec)) return 0;
	vec->size -= size;
	return 1;
}

static inline int set_vec_size(Vector *vec, size_t size){
	return vec->size = size;
}

static inline void push_vec(Vector *vec, void *e){
	encap_vec(vec, 1);
	vec_memcpy(vec->buffer + (vec_size(vec)) * vec->e_size, e, vec->e_size);
	vec->size ++;
}

#define gpush_vec(vec, v, data_type) (encap_vec(vec, 1), (((data_type *)(vec)->buffer)[vec_size(vec)] = (v)), (vec)->size ++)

static inline int pop_vec(Vector *vec, void *e){
	if(reduce_vec_size(vec, 1) == 0) return 0;
	vec_memcpy(e, vec->buffer + vec_size(vec) * vec->e_size, vec->e_size);
	return 1;
}

#define gpop_vec(vec, v, data_type) (reduce_vec_size(vec, 1)? (v = ((data_type *)(vec)->buffer)[vec_size(vec)], 1) : 0)

static inline void set_vec(Vector *vec, size_t idx, void *e){
	vec_memcpy(vec->buffer + idx * vec->e_size, e, vec->e_size);
}

#define gset_vec(vec, idx, v, data_type) ((data_type*)(vec)->buffer)[idx] = v

static inline int get_vec(Vector *vec, size_t idx, void *e){
	vec_memcpy(e, vec->buffer + idx * vec->e_size, vec->e_size);
	return 1;
}

#define gget_vec(vec, idx, data_type) ((data_type *)(vec)->buffer)[idx]

static inline void* get_vec_ref(Vector *vec, size_t idx){
	return vec->buffer + idx * vec->e_size;
}

#define gpeer_vec(vec, data_type) (vec_size(vec)? ((data_type *)(vec)->buffer)[vec_size(vec) - 1] : 0)

static inline void* get_last_vec_ref(Vector *vec){
	if(vec_size(vec)) return vec->buffer + (vec_size(vec) - 1) * vec->e_size;
	else return NULL;
}

static inline void* get_next_vec_ref(Vector *vec){
	add_vec_size(vec, 1);
	return vec->buffer + (vec_size(vec) - 1) * vec->e_size;
}

static inline void qsort_vec(Vector *vec, cmp_vec_fun fun){
	qsort(vec->buffer, vec_size(vec), vec->e_size, fun);
}

static inline void* bsearch_vec(Vector *vec, void *q, cmp_vec_fun fun){
	return bsearch(q, vec->buffer, vec_size(vec), vec->e_size, fun);
}

#define search_array(uniq_flag, array, size, key, val_macro, ret) long i##uniq_flag, j##uniq_flag, m##uniq_flag;\
	i##uniq_flag = 0; j##uniq_flag = size; while(i##uniq_flag < j##uniq_flag){\
	m##uniq_flag = i##uniq_flag + (j##uniq_flag - i##uniq_flag) / 2;\
	if(val_macro((array)[m##uniq_flag]) < val_macro(key)) i##uniq_flag = m##uniq_flag + 1;\
	else j##uniq_flag = m##uniq_flag;\
	}\
	if(i##uniq_flag < (long)size && val_macro((array)[i##uniq_flag]) == val_macro(key)) ret = i##uniq_flag; \
	else ret = - (i##uniq_flag + 1)

static inline void reverse_vec(Vector *vec){
	size_t i, j;
	void *buf;
	if(vec_size(vec) == 0) return;
	buf = malloc(vec->e_size);
	i = 0;
	j = vec_size(vec) - 1;
	while(i < j){
		vec_memcpy(buf, vec->buffer + i * vec->e_size, vec->e_size);
		vec_memcpy(vec->buffer + i * vec->e_size, vec->buffer + j * vec->e_size, vec->e_size);
		vec_memcpy(vec->buffer + j * vec->e_size, buf, vec->e_size);
		i ++;
		j --;
	}
	free(buf);
}

static inline int cat_vec(Vector *dst, Vector *src){
	if(dst->e_size != src->e_size) return -1;
	add_vec_size(dst, vec_size(src));
	vec_memcpy(dst->buffer + (vec_size(dst) - vec_size(src)) * dst->e_size, src->buffer, vec_size(src) * src->e_size);
	return 0;
}

#define clear_vec(vec) set_vec_size(vec, 0)

static inline void reset_vec(Vector *vec){
	set_vec_size(vec, 0);
	memset(vec->buffer, 0, vec->e_size * vec->cap);
}

static inline void free_vec(Vector *vec){
	if(vec == NULL) return;
	free(vec->buffer);
	free(vec);
}

#endif
