/*
 * 
 * Copyright (c) 2011, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#ifndef __STRING_RJ_H
#define __STRING_RJ_H

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include "vector.h"

/**
 * String
 */

#ifndef SWAP_TMP
#define SWAP_TMP
#define swap_tmp(a, b, t) t = a; a = b; b = t
#endif

typedef struct {
	char *string;
	int size;
	int capacity;
} String;

typedef struct {
	char *string;
	int  size;
} VirtualString;

#define uc(ch) (((ch) >= 'a' && (ch) <= 'z')? (ch) + 'A' - 'a' : (ch))
#define lc(ch) (((ch) >= 'A' && (ch) <= 'Z')? (ch) + 'a' - 'A' : (ch))

static inline String* init_string(int cap){
	String *str;
	str = (String*)malloc(sizeof(String));
	str->size = 0;
	str->capacity = (cap&0x1)? cap:cap+1;
	str->string = (char*)malloc(sizeof(char) * (str->capacity + 1));
	str->string[0] = 0;
	return str;
}

static inline void encap_string(String *str, int inc){
	if(inc + str->size >= str->capacity){
		if(inc < str->size) str->capacity = str->size * 2 + 1;
		else str->capacity = inc * 2 + 1;
		str->string = (char*)realloc(str->string, str->capacity + 1);
	}
}

static inline void uc_string(String *str){
	int i;
	for(i=0;i<str->size;i++){
		if(str->string[i] >= 'a' && str->string[i] <= 'z') str->string[i] = str->string[i] + 'A' - 'a';
	}
}

static inline void lc_string(String *str){
	int i;
	for(i=0;i<str->size;i++){
		if(str->string[i] >= 'A' && str->string[i] <= 'Z') str->string[i] = str->string[i] + 'a' - 'A';
	}
}

static inline char* substr(char *string, int start, int end, char *dst){
	int i, size;
	char *str;
	size = strlen(string);
	if(start > size) start = size;
	else if(start < 0) start = 0;
	if(end > size) end = size;
	else if(end < 0) end = 0;
	size = end - start;
	if(size < 0) size = 0;
	if(dst != NULL) str = dst;
	else str = (char*)malloc(sizeof(char) * (size + 1));
	for(i=start;i<end;i++){
		str[i-start] = string[i];
	}
	str[size] = '\0';
	return str;
}

static inline char* catstr(int n_str, ...){
	char *str, *s;
	int i, len;
	va_list params;
	
	len = 0;
	str = NULL;
	va_start(params, n_str);
	for(i=0;i<n_str;i++){
		s = va_arg(params, char*);
		len += strlen(s);
		str = realloc(str, len + 1);
		if(i == 0) str[0] = 0;
		strcat(str, s);
	}
	va_end(params);
	return str;
}

static inline void chomp_string(String *str){
	if(str->size && str->string[str->size - 1] == '\n'){
		str->size --;
		str->string[str->size] = 0;
	}
}

static inline void chomp_vstring(VirtualString *str){
	if(str->size && str->string[str->size - 1] == '\n'){
		str->size --;
	}
}

static inline void trim_string(String *str){
	int i, j;
	i = str->size - 1;
	while(i >= 0 && (str->string[i] == '\n' || str->string[i] == '\t' || str->string[i] == ' ')) i--; 
	str->size = i + 1;
	i = 0;
	while(i < str->size && (str->string[i] == '\n' || str->string[i] == '\t' || str->string[i] == ' ')) i++;
	if(i){
		for(j=i;j<str->size;j++){ str->string[j-i] = str->string[j]; }
		str->size -= i;
	}
	str->string[str->size] = 0;
}

static inline void trim_vstring(VirtualString *str){
	int i;
	i = str->size - 1;
	while(i >= 0 && (str->string[i] == '\n' || str->string[i] == '\t' || str->string[i] == ' ')) i--; 
	str->size = i + 1;
	i = 0;
	while(i < str->size && (str->string[i] == '\n' || str->string[i] == '\t' || str->string[i] == ' ')) i++;
	str->string += i;
}

static inline void append_string(String *str, char *src, int offlen){
	int i;
	encap_string(str, offlen);
	for(i=0;i<offlen;i++) str->string[str->size + i] = src[i];
	str->size += offlen;
	str->string[str->size] = 0;
}

static inline String* as_string(char *chs){
	int len;
	String *str;
	len = strlen(chs);
	str = init_string(len);
	append_string(str, chs, len);
	return str;
}

static inline VirtualString* as_vstring(char *chs){
	int len;
	VirtualString *str;
	len = strlen(chs);
	str = malloc(sizeof(VirtualString));
	str->string = chs;
	str->size = len;
	return str;
}

static inline void add_char_string(String *str, char ch){
	encap_string(str, 1);
	str->string[str->size] = ch;
	str->size ++;
	str->string[str->size] = 0;
}

static inline void clear_string(String *str){ str->size = 0; }

static inline int split_string(String *str, char separator, Vector *virtual_strings){
	VirtualString *vstr;
	int n_tab, i, s;
	n_tab = 0;
	i = 0;
	s = 0;
	while(i <= str->size){
		if(i == str->size || str->string[i] == separator){
			vstr = get_next_vec_ref(virtual_strings);
			vstr->string = str->string + s;
			n_tab ++;
			vstr->size = i - s;
			s = i + 1;
		}
		i ++;
	}
	return n_tab;
}

static inline int split_vstring(VirtualString *str, char separator, Vector *virtual_strings, int cut){
	VirtualString *vstr;
	int n_tab, i, s;
	n_tab = 0;
	i = 0;
	s = 0;
	while(i <= str->size){
		if(i == str->size || str->string[i] == separator){
			if(cut) str->string[i] = '\0';
			vstr = get_next_vec_ref(virtual_strings);
			vstr->string = str->string + s;
			n_tab ++;
			vstr->size = i - s;
			s = i + 1;
		}
		i ++;
	}
	return n_tab;
}

static inline void reverse_string(String *str){
	int i, j;
	char c;
	i = 0;
	j = str->size - 1;
	while(i < j){
		swap_tmp(str->string[i], str->string[j], c);
		i ++;
		j --;
	}
}

static inline void trunc_string(String *str, int size){
	if(size >= str->size) return;
	str->size = size;
	str->string[size] = 0;
}

static inline String* clone_string(String *str){
	String *clone;
	clone = init_string(str->size);
	append_string(clone, str->string, str->size);
	return clone;
}

static inline void free_string(String *str){ free(str->string); free(str); }

static inline void free_vstring(VirtualString *str){ free(str); }

#endif
