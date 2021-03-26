/*
 * $Id: flexi.c,v 1.4 2004/08/04 16:25:57 ewalter Exp $
 */

/****************************************************************************
*  'FlexiLib' - Flexible Input Library                                      *
*  Copyright (C) 2000  Gerhard J. Theurich (gjt)                            *
*                                                                           *
*  This program is free software; you can redistribute it and/or modify     *
*  it under the terms of the GNU General Public License as published by     *
*  the Free Software Foundation; either version 2 of the License, or        *
*  (at your option) any later version.                                      *
*                                                                           *
*  This program is distributed in the hope that it will be useful,          *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
*  GNU General Public License for more details.                             *
*                                                                           *
*  You should have received a copy of the GNU General Public License        *
*  along with this program; if not, write to the Free Software              *
*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA  *
****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "flexi.h"


typedef struct keycell{

  char key[80];             /* key */
  struct keycell *next;     /* next keycell */
  struct keycell *last;     /* last keycell */
  int n;                    /* number of paramters */
  int resolved;             /* resolved flag: 0 not resolved, 1 resolved */
  int required;             /* required flag: 0 not required, 1 required */
  char **param_t;           /* parameter type array */ 
  void **param_p;           /* parameter pointer array */  
  
}keycell;


static keycell *keyqueue;   /* pointer onto the current keyqueue */
static int number_of_keys;  /* number of keys in current keyqueue */


static int scan_key(FILE *fp, char *key);
  /* Scan next key in 'fp', return in 'key'. Returns 0 if no problem was   */
  /* detected and 1 if no key was found.                                   */

static keycell *find_key(char *key);       
  /* Find 'key' in current keyqueue. Returns the cell where 'key' was      */
  /* found if it was found and NULL if 'key' was not found.                */


int flexi_request_key(char *key, int required, char *fmt, ...){
  /* Request to add 'key' with format string 'fmt' into the current        */
  /* keyqueue. Returns 0 if no problem was detected and 1 if there were    */
  /* less arguments given than 'fmt' required.                             */
  
  int i, n;
  char *c;
  char *t;
  keycell *keyc;
  keycell *cell;
  char **param_t;
  void **param_p;
  va_list ap;
  
  /* scan through 'fmt' and count '%' characters */
  n=0;
  for (c = fmt; *c; c++)
    if (*c == '%') n++;
  
  /* don't bother if nothing was requested for this key */
  if (n == 0) return 0;
  
  /* check if this key already exists in the current keyqueue */
  if ((keyc=find_key(key))){
    /* copy previously requested arguments */
    param_t = keyc->param_t;
    param_p = keyc->param_p;
    keyc->param_t = (char **)malloc((keyc->n+n)*sizeof(char *));
    keyc->param_p = (void **)malloc((keyc->n+n)*sizeof(void *));
    for (i=0; i<keyc->n; i++){
      keyc->param_t[i] = param_t[i];
      keyc->param_p[i] = param_p[i];
    }
    free(param_t);
    free(param_p);
    keyc->n += n;
    if (!keyc->required && required){
      keyc->required = required;
      ++number_of_keys;
    }else if (keyc->required && !required){
      keyc->required = required;
      --number_of_keys;
    } 
  }else{
    /* prepare the keycell */
    keyc = (keycell *)malloc(sizeof(keycell)); 
    keyc->n = n;
    strncpy(keyc->key, key, 80);
    keyc->param_t = (char **)malloc(keyc->n*sizeof(char *));
    keyc->param_p = (void **)malloc(keyc->n*sizeof(void *));
    keyc->required = required;
    if (keyc->required)
      ++number_of_keys;
  }
  
  va_start(ap, fmt); /* get started */
  
  /* loop over the arguments */  
  c = fmt;
  while (*c!='%') c++;
  for (i=keyc->n-n; i<keyc->n; i++){
    keyc->param_t[i] = (char *)malloc(5*sizeof(char));
    keyc->param_p[i] = va_arg(ap, void *);
    if (keyc->param_p[i]==NULL) return 1;  /* not enough arguments in call */ 
    t = keyc->param_t[i];
    /* copy format string */
    do *t++ = *c++;
    while (*c!='%' && *c);
    *t = '\0';
    
    /*    printf("%p, type: %s, pointer: %p\n", 
	  keyc, keyc->param_t[i], keyc->param_p[i]); */
    
  }
  
  va_end(ap); /* clean up */
  
  /* keycell is complete -> add it to current keyqueue if not yet present */
  if (!find_key(key)){
    if (keyqueue == NULL){
      keyqueue = keyc;
      keyc->next = NULL;
      keyc->last = NULL;
    }else{
      cell = keyqueue;
      while (cell->next != NULL) cell = cell->next;
      cell->next = keyc;
      keyc->next = NULL;
      keyc->last = cell;
    }
  }
  return 0;
}


int flexi_gather_keys(FILE *fp){
  /* Gather all requested keys in current keyqueue. The corresponding      */
  /* arguments are set to those found at the last matching key in 'fp'.    */
  /* Returns the number of required keys in current keyqueue that could    */
  /* not be resolved in 'fp'.                                              */
  
  int i;
  int keys_gathered=0;
  int c;
  char key[80];
  keycell *keyc;
  
  while (!scan_key(fp, key)){
    if ((keyc=find_key(key))){
      for (i=0; i<keyc->n; i++){
        do{
          /* scanning until found some real stuff */
          while(((c=getc(fp)) != EOF) && ((c =='\n') || (c==' ') || (c=='\t')));
          if (c=='#'){
            /* skip the rest of this line */
            while(((c=getc(fp)) != EOF) && (c!='\n'));
          }
        }while ((c=='\n') || (c==' ') || (c=='\t') || (c=='#'));
        /* worth scanning for */        
        ungetc(c, fp);
        if (fscanf(fp, keyc->param_t[i], (void *)keyc->param_p[i])==EOF)
          break;
      }
      if (i==keyc->n){
        keyc->resolved = 1;
        if (keyc->required) ++keys_gathered;
      }
    }
  }
   
  return number_of_keys-keys_gathered;
}


int flexi_resolved_key(char *key){
  /* Returns 1 if 'key' was resolved by a previous 'flexi_gather_keys()'   */
  /* call, 0 if 'key' was not resolved and -1 if 'key' does not exist in   */
  /* current keyqueue.                                                     */
  
  keycell *keyc = find_key(key);
  
  if (keyc)
    return keyc->resolved;
  else
    return -1;
}


void flexi_clear_keys(void){
  /* Clear the entire current keyqueue.                                    */
  
  int i;
  keycell *cell = keyqueue;

  if (cell){
  
    while (cell->next != NULL) cell = cell->next;
  
    while (keyqueue){
      for (i=0; i<cell->n; i++)
        free(cell->param_t[i]);
      free(cell->param_t);
      free(cell->param_p);
      if (cell->last != NULL){
        cell = cell->last;
        free(cell->next);
      }else{
        free(cell);
        keyqueue = NULL;
      } 
    }
  }
  number_of_keys = 0;
  return;
}


static int scan_key(FILE *fp, char *key){
  /* Scan next key in 'fp', return in 'key'. Returns 0 if no problem was   */
  /* detected and 1 if no key was found.                                   */

  int i;
  int c;
  
  key[0] = '\0';  /* erase string */
  /* repeat until '[' found */
  do{
    while(((c=getc(fp)) != EOF) && (c!='[') && (c!='#'));  
    if (c==EOF) return 1;                    /* reached end of file */
    if (c=='#') while(((c=getc(fp)) != EOF) && (c!='\n')); /* skip comments */
  }while(c!='[');
  i=0;
  while(((c=getc(fp)) != EOF) && (c!=']') && (c!='[')){
    key[i++]=c;
  }
  key[i]='\0';
  return 0;  
}


static keycell *find_key(char *key){
  /* Find 'key' in current keyqueue. Returns the cell where 'key' was      */
  /* found if it was found and NULL if 'key' was not found.                */
  
  keycell *cell = keyqueue;
  
  if (cell==NULL) return NULL;
  
  do{
    if (strcasecmp(cell->key,key)==0) break;
    cell = cell->next;
  }while (cell != NULL);
  
  return cell;
}

