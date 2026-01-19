/*
 * CSE 351 Lab 1b (Manipulating Bits in C)
 *
 * Name(s):  Luke Verlangieri
 * NetID(s): lverl23
 *
 * This is a file for managing a store of various aisles, represented by an
 * array of 64-bit integers. See aisle_manager.c for details on the aisle
 * layout and descriptions of the aisle functions that you may call here.
 *
 * Written by Porter Jones (pbjones@cs.washington.edu)
 */

#include <stddef.h>  // To be able to use NULL
#include "aisle_manager.h"
#include "store_client.h"
#include "store_util.h"

// Number of aisles in the store
#define NUM_AISLES 10

// Number of sections per aisle
#define SECTIONS_PER_AISLE 4

// Number of items in the stockroom (2^6 different id combinations)
#define NUM_ITEMS 64

// Number of items in a section
#define ITEMS_PER_SECTION 10

// Global array of aisles in this store. Each unsigned long in the array
// represents one aisle.
unsigned long aisles[NUM_AISLES];

// Array used to stock items that can be used for later. The index of the array
// corresponds to the item id and the value at an index indicates how many of
// that particular item are in the stockroom.
int stockroom[NUM_ITEMS];


/* Starting from the first aisle, refill as many sections as possible using
 * items from the stockroom. A section can only be filled with items that match
 * the section's item id. Prioritizes and fills sections with lower addresses
 * first. Sections with lower addresses should be fully filled (if possible)
 * before moving onto the next section.
 */
void refill_from_stockroom() {
  for(int i = 0; i < NUM_AISLES; i++) { //Lowest address section is 
    unsigned long* aisle = &aisles[i]; //Grab the address of a the right aisle
    for(int j = 0; j < SECTIONS_PER_AISLE; j++){ //go through sections
      unsigned short id = get_id(aisle, j);                // pull the id
      unsigned short items = num_items(aisle, j);         // pull the number of items in an aisle
      unsigned short slots_avaliable = (ITEMS_PER_SECTION - items); //find out how many spaces are avaliable

      int* items_in_stock = &stockroom[id];               //pointer to the items that are in stock

      if(*items_in_stock < slots_avaliable) { //if there are more slots than items in stock, add what we have
        add_items(aisle, j, *items_in_stock);
        *items_in_stock = 0; 
      } else { //if there are more items in stock than slots avaliable, only add what is needed
        add_items(aisle, j, slots_avaliable);
        *items_in_stock -= slots_avaliable; 
      }
    }
  }
}

/* Remove at most num items from sections with the given item id, starting with
 * sections with lower addresses, and return the total number of items removed.
 * Multiple sections can store items of the same item id. If there are not
 * enough items with the given item id in the aisles, first remove all the
 * items from the aisles possible and then use items in the stockroom of the
 * given item id to finish fulfilling an order. If the stockroom runs out of
 * items, you should remove as many items as possible.
 */
int fulfill_order(unsigned short id, int num) {
  int items_removed = 0; 
  for(int i = 0; i < NUM_AISLES; i++){ //start at lowest aisle addr
    unsigned long* aisle = &aisles[i];
    for(int j = 0; j < SECTIONS_PER_AISLE; j++){ //start at lowest section addr per aisle
      if(id == get_id(aisle, j)){
        unsigned short items = num_items(aisle, j);
        if(num >= items) { //remove all items if there is more requested than avaliable
          remove_items(aisle, j, items);
          items_removed += items; 
          num -= items; 
        } else { //if more items in section than requested remove as much as we can
          remove_items(aisle, j, num);
          items_removed += num;
          num = 0; 
        }
      }
      if(num == 0) return items_removed; // if we hit our quota we can return immedaitely 
    }
  }

  //because we return when we get num == 0 in the loop we will always have num > 0 here
  int* items_in_stock = &stockroom[id]; //reference to the stockroom
  if(*items_in_stock >= num){ //if we have enough in the stock room
    items_removed += num; 
    *items_in_stock -= num; 
  } else {                    //if we dont have enough in the stock room
    items_removed += *items_in_stock;
    *items_in_stock = 0;
  }

  return items_removed;
}

/* Return a pointer to the first section in the aisles with the given item id
 * that has no items in it or NULL if no such section exists. Only consider
 * items stored in sections in the aisles (i.e., ignore anything in the
 * stockroom). Break ties by returning the section with the lowest address.
 */
unsigned short* empty_section_with_id(unsigned short id) {
  
  for(int i = 0; i < NUM_AISLES; i++){
    unsigned long* aisle = &aisles[i];
    for(int j = 0; j < SECTIONS_PER_AISLE; j++){
      if(id == get_id(aisle, j)){                 //find matching ID
        if(num_items(aisle, j) == 0){             //check if empty
          return ((unsigned short*)aisle + j);    //return a pointer to the corresponding section (first found is lowest addr)
        } 
      }
    }
  }

  return NULL;
}

/* Return a pointer to the section with the most items in the store. Only
 * consider items stored in sections in the aisles (i.e., ignore anything in
 * the stockroom). Break ties by returning the section with the lowest address.
 */
unsigned short* section_with_most_items() {

  // if num_items are equal keep the original pointer (has the lowest addr)

  //init pointer and counter of most_items
  unsigned short* section = (unsigned short*)aisles; //points to the first section  
  unsigned short most_items = 0;

  for(int i = 0; i < NUM_AISLES; i++){
    unsigned long* aisle = &aisles[i];
    for(int j = 0; j < SECTIONS_PER_AISLE; j++){  
      unsigned short items = num_items(aisle, j);     //find matching ID
      if(items > most_items){                         //check if more items (preserve the lower addr pointer)
        section = ((unsigned short*)aisle + j);       //update pointer to the corresponding section
        most_items = items;
      }
    }
  }

  return section;
}
