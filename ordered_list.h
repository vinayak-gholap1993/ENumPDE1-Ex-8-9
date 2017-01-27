#ifndef ORDERED_LIST_H
#define ORDERED_LIST_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * An element in the ordered list
 */
typedef struct ol_element {
  size_t val;
  struct ol_element * next;
} ordered_list_element;

/*
 * The ordered list object
 */
typedef struct  {
  ordered_list_element * begin;
  size_t length;
} ordered_list;

/*
 * Initializes the list object to an empty list
 */
void init_list(ordered_list * list);

/*
 * Inserts a new element with the given value in the correct position.
 * Does nothing if the element already exists.
 */
void insert_list(ordered_list * list, size_t val);

/*
 * Prints the ordered list to stdout
 */
void print_list(ordered_list const * list);

/*
 * Deallocates all data structures of the ordered list (but not the ordered
 * list object itself)
 */
void clear_list(ordered_list * list);

#ifdef __cplusplus
}
#endif

#endif
