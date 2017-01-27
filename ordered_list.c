#include <stdio.h>
#include <stdlib.h>

#include "ordered_list.h"

#define err_exit(msg) \
  { fprintf(stderr, "ERROR: %s (%s:%d)\n", msg, __FILE__, __LINE__); exit(1); }

void init_list(ordered_list * list)
{
  if (list == NULL)
    err_exit("No list given!");

  list->begin = NULL;
  list->length = 0;
}

void insert_list(ordered_list * list, size_t val)
{
  ordered_list_element * it, * tmp;

  if (list == NULL)
    err_exit("No list given!");

  /* Find insert position */
  it = list->begin;
  if (it != NULL)
  {
    while (it->next != NULL && it->next->val <= val)
    {
      it = it->next;
    }
  }

  /* insert at the beginning */
  if (it == list->begin && (it == NULL || it->val > val) )
  {
    tmp = list->begin;
    list->begin = (ordered_list_element *) malloc(sizeof(ordered_list_element));
    if (list->begin == NULL) err_exit("Insert into ordered list failed!");
    list->begin->val = val;
    list->begin->next = tmp;
    ++list->length;
  }
  /* insert after current element */
  else if (it->val != val)
  {
    tmp = it->next;
    it->next = (ordered_list_element *) malloc(sizeof(ordered_list_element));
    if (it->next == NULL) err_exit("Insert into ordered list failed!");
    it->next->val = val;
    it->next->next = tmp;
    ++list->length;
  }
}

void print_list(ordered_list const * list)
{
  ordered_list_element const * it;

  if (list == NULL)
    err_exit("No list given!");

  it = list->begin;
  if (it != NULL)
  {
    while(it != NULL)
    {
      printf("%d ", (int)it->val);
      it = it->next;
    }
  }
  printf("\n\n");
}

void clear_list(ordered_list * list)
{
  ordered_list_element *it, *tmp;

  it = list->begin;
  if (it != NULL)
  {
    while (it != NULL)
    {
      tmp = it;
      it = it->next;
      free(tmp);
    }
  }

  list->begin = NULL;
  list->length = 0;
}
