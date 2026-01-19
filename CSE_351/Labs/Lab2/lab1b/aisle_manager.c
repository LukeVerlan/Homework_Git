/*
 * CSE 351 Lab 1b (Manipulating Bits in C)
 *
 * Name(s): Luke Verlangieri
 * NetID(s): lverl23
 *
 * ----------------------------------------------------------------------------
 * Overview
 * ----------------------------------------------------------------------------
 *  This is a program to keep track of the items in a small aisle of a store.
 *
 *  A store's aisle is represented by a 64-bit long integer, which is broken
 *  into 4 16-bit sections representing one type of item each. Note that since
 *  a section is 16-bits, it fits nicely into C's short datatype.
 *
 *  Aisle Layout:
 *
 *    Within an aisle, sections are indexed starting with the least-significant
 *    section being at index 0 and continuing up until one less than the number
 *    of sections.
 *
 *    Example aisle:
 *
 *                MSB                                                       LSB
 *                  +-------------+-------------+-------------+-------------+
 *                  |  Section 3  |  Section 2  |  Section 1  |  Section 0  |
 *                  +-------------+-------------+-------------+-------------+
 *                  |             |             |             |             |
 *      bit offset: 64            48            32            16            0
 *
 *  Section Layout:
 *
 *    A section in an aisle is broken into 2 parts. The 6 most significant bits
 *    represent a unique identifier for the type of item stored in that
 *    section. The rest of the bits in a section (10 least significant)
 *    indicate individual spaces for items in that section. For each of the 10
 *    bits/spaces, a 1 indicates that an item of the section's type is stored
 *    there and a 0 indicates that the space is empty.
 *
 *    Example aisle section: 0x651A
 *
 *                MSB                               LSB
 *                  +-----------+-------------------+
 *                  |0 1 1 0 0 1|0 1 0 0 0 1 1 0 1 0|
 *                  +-----------+-------------------+
 *                  | item id   | section spaces    |
 *      bit offset: 16          10                  0
 *
 *      In this example, the item id is 0b011001, and there are currently 4
 *      items stored in the section (at bit offsets 8, 4, 3, and 1) and 6
 *      vacant spaces.
 *
 *  Written by Porter Jones (pbjones@cs.washington.edu)
 */

#include "aisle_manager.h"
#include "store_util.h"

// the number of total bits in a section
#define SECTION_SIZE 16

// The number of bits in a section used for the item spaces
#define NUM_SPACES 10

// The number of bits in a section used for the item id
#define ID_SIZE 6

// The number of sections in an aisle
#define NUM_SECTIONS 4

// TODO: Fill in these with the correct hex values

// Mask for extracting a section from the least significant bits of an aisle.
// (aisle & SECTION_MASK) should preserve a section's worth of bits at the
// lower end of the aisle and set all other bits to 0. This is essentially
// extracting section 0 from the example aisle shown above.
#define SECTION_MASK 0xFFFF // A section is all 16 bits

// Mask for extracting the spaces bits from a section.
// (section & SPACES_MASK) should preserve all the spaces bits in a section and
// set all non-spaces bits to 0.
#define SPACES_MASK 0x3FF // Spaces take up the 10 LSB from a section

// Mask for extracting the ID bits from a section.
// (section & ID_MASK) should preserve all the id bits in a section and set all
// non-id bits to 0.
#define ID_MASK 0xFC00  // ID takes up indexes 10-15 


/* Given a pointer to an aisle and a section index, return the section at the
 * given index of the given aisle.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
unsigned short get_section(unsigned long* aisle, int index) {
  //Dereference to get the binary value of the aisle 
  //Shift right the desired index of the section times the size of each section 
  // This value is implicitly casted to a unsigned short, removing the other bits
  return ((*aisle) >> (index * SECTION_SIZE));
}

/* Given a pointer to an aisle and a section index, return the spaces of the
 * section at the given index of the given aisle. The returned short should
 * have the least 10 significant bits set to the spaces and the 6 most
 * significant bits set to 0.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
unsigned short get_spaces(unsigned long* aisle, int index) {
  unsigned short section = get_section(aisle, index); //Get the section at the right index
  return section & SPACES_MASK; //Mask the values to return only the spaces
}

/* Given a pointer to an aisle and a section index, return the id of the
 * section at the given index of the given aisle. The returned short should
 * have the least 6 significant bits set to the id and the 10 most significant
 * bits set to 0.
 *
 * Example: if the section is 0b0110010100011010, return 0b0000000000011001.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
unsigned short get_id(unsigned long* aisle, int index) {
  unsigned short section = get_section(aisle,index); //get the section at the right index
  //Mask the values and shift over by the number of spaces to put the ID in the rightmost index 
  return (section & ID_MASK) >> NUM_SPACES; 
}

/* Given a pointer to an aisle, a section index, and a short representing a new
 * bit pattern to be used for section, set the section at the given index of
 * the given aisle to the new bit pattern.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
void set_section(unsigned long* aisle, int index, unsigned short new_section) {

  //Step 1 - clear the area to change with 0s. 
  // cast the mask to an unsigned long to match the size of the aisle, left shift that 
  // value over by the section index * the section bit size, we want to preserve everything
  // but what is inside the mask, so we need to not this value. Then "and" that value with the
  // dereferenced aisle pointer value to clear that section, and set the dereferenced aisle pointer value
  // to be the one with the cleared section to not allocate new memory
  *aisle = (*aisle & (~((unsigned long)SECTION_MASK << (index * SECTION_SIZE))));

  //Step 2 - reassign the area with the new section value 
  // cast the short to an unsigned long to make it an equivalent amount of bits from 16 to 64 (padded with 0s)
  // left shift it to the same section using the same principle as above. 
  // because its all zeros, or the bits onto the now cleared section 
  // set the dereferenced aisle pointer to be that value
  *aisle = (*aisle | ((unsigned long)new_section << (index * SECTION_SIZE)));
}

/* Given a pointer to an aisle, a section index, and a short representing a new
 * bit pattern to be used for the spaces of the section, set the spaces for the
 * section at the given index of the given aisle to the new bit pattern.
 * 
 * If new_spaces is invalid (if it has 1s outside of its 10 least significant
 * bits), this method should return without modifying anything.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
void set_spaces(unsigned long* aisle, int index, unsigned short new_spaces) {
  //Check to make sure that new spaces only has bits in the spaces section
  if((new_spaces & (~SPACES_MASK)) == 0) {
    //same as above, only this time clearing spaces and not the whole section
    *aisle = (*aisle & (~((unsigned long)SPACES_MASK << (index * SECTION_SIZE))));

    //same idea as above function 
    *aisle = (*aisle | ((unsigned long)new_spaces << (index * SECTION_SIZE)));
  }
}

/* Given a pointer to an aisle, a section index, and a short representing a new
 * bit pattern to be used for the id of the section, set the id for the section
 * at the given index of the given aisle to the new bit pattern.
 * 
 * If new_id is invalid (if it has 1s outside of its 6 least significant
 * bits), this method should return without modifying anything.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
void set_id(unsigned long* aisle, int index, unsigned short new_id) {
  // defensive cast to unsigned to make a logical shift
  if((new_id & (~((unsigned)ID_MASK >> NUM_SPACES))) == 0){
    //same as above, only this time clearing the ID and not the whole section
    *aisle = (*aisle & (~(((unsigned long)ID_MASK) << (index * SECTION_SIZE))));
    
    //Need to logical shift left to put in the right index
    *aisle = (*aisle | ((unsigned long)(new_id << NUM_SPACES) << (index * SECTION_SIZE)));
  }
}

/* Given a pointer to an aisle, a section index, and a space index, toggle the
 * item in the given space index of the section at the given section index in
 * the given aisle. Toggling means that if the space was previously empty, it
 * should now be filled with an item, vice-versa.
 *
 * Can assume the section index is a valid index (0-3 inclusive).
 * Can assume the spaces index is a valid index (0-9 inclusive).
 */
void toggle_space(unsigned long* aisle, int index, int space_index) {
  unsigned short spaces = get_spaces(aisle, index); //get the section
  spaces = spaces ^ (0x1 << space_index); //mask the right bits and flip the unmasked bit
  set_spaces(aisle, index, spaces); // reset the spaces
}

/* Given a pointer to an aisle and a section index, return the number of items
 * in the section at the given index of the given aisle.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
unsigned short num_items(unsigned long* aisle, int index) {
  unsigned short spaces = get_spaces(aisle, index);
  char items = 0; 
  for(char i = 0; i < NUM_SPACES; i++) items += ((spaces >> i) & 0x1); //extract each bit value, any bit holding a 1 will be added
  return items;
}

/* Given a pointer to an aisle, a section index, and the desired number of
 * items to add, add at most the given number of items to the section at the
 * given index in the given aisle. Items should be added to the least
 * significant spaces possible. If n is larger than or equal to the number of
 * empty spaces in the section, then the section should appear full after the
 * method finishes.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
void add_items(unsigned long* aisle, int index, int n) {
  unsigned short spaces = get_spaces(aisle, index); //get the spaces
  for(char i = 0; i < NUM_SPACES; i++) { //iterate through the spaces
    if (n == 0) break;                   //break if out of items
    if (((spaces >> i) & 0x01) == 0) {   //check if index is open
      spaces = spaces | (0x01 << i);     //add if open
      n--;                               //deincrement an item
    }
  }
  set_spaces(aisle,index,spaces);       //set the new spaces
}


/* Given a pointer to an aisle, a section index, and the desired number of
 * items to remove, remove at most the given number of items from the section
 * at the given index in the given aisle. Items should be removed from the
 * least significant spaces possible. If n is larger than or equal to the
 * number of items in the section, then the section should appear empty after
 * the method finishes.
 *
 * Can assume the index is a valid index (0-3 inclusive).
 */
void remove_items(unsigned long* aisle, int index, int n) {
  unsigned short spaces = get_spaces(aisle, index); //get the spaces
  for(char i  = 0; i < NUM_SPACES; i++) { //iterate through the spaces
    if (n == 0) break;                   //break if removed enough items
    if (((spaces >> i) & 0x01) == 0x01) {   //check if index is full
      spaces = spaces ^ (0x01 << i);     //remove if full (1 ^ 1 = 0)
      n--;                               //deincrement an item
    }
  }
  set_spaces(aisle,index,spaces);       //set the new spaces
}

/* Given a pointer to an aisle, a section index, and a number of slots to
 * rotate by, rotate the items in the section at the given index of the given
 * aisle to the left by the given number of slots.
 *
 * Example: if the spaces are 0b0111100001, then rotating left by 2 results
 *          in the spaces     0b1110000101
 *
 * Can assume the index is a valid index (0-3 inclusive).
 * Can NOT assume n < NUM_SPACES (hint: find an equivalent rotation).
 */
void rotate_items_left(unsigned long* aisle, int index, int n) {
  n = n % NUM_SPACES; // clean to an equivalent rotation (going around NUM_SPACES times is a circle)
  unsigned short spaces = get_spaces(aisle, index); //get spaces
  unsigned short rollover = 0;                      //init a rollover 
  for(char i = 1; i < NUM_SPACES; i++) {            //iterate through the spaces
    if(n == 0) break; 
    spaces = (spaces << 1);                              // Rotate left
    char rollover_bit = ((spaces >> NUM_SPACES) & 0x01); // Capture the rolled over bit
    rollover = (rollover << 1) | rollover_bit;           // add the rolled over bit to the LSB of the rollover bucket
    n--; 
  }

  //The there will be 0s equal to the rotation amount % 10 from the LSB in spaces, therefore just | whatever was put into the rollover. 
  //Need to and away any left over bits (if n % 10 > 6 bits in spaces will be dropped off the face of the earth)
  set_spaces(aisle, index, ((spaces | rollover) & SPACES_MASK)); 
}

/* Given a pointer to an aisle, a section index, and a number of slots to
 * rotate by, rotate the items in the section at the given index of the given
 * aisle to the right by the given number of slots.
 *
 * Example: if the spaces are 0b1000011110, then rotating right by 2 results
 *          in the spaces     0b1010000111
 *
 * Can assume the index is a valid index (0-3 inclusive).
 * Can NOT assume n < NUM_SPACES (hint: find an equivalent rotation).
 */
void rotate_items_right(unsigned long* aisle, int index, int n) {
  n = n % NUM_SPACES; // clean to an equivalent rotation (going around NUM_SPACES times is a circle)
  unsigned short spaces = get_spaces(aisle, index); //get spaces
  unsigned short rollover = 0;                      //init a rollover 
  for(char i = 1; i < NUM_SPACES; i++) {            //iterate through the spaces
    if(n == 0) break; 
    unsigned short rollover_bit = (spaces & 0x01);       // capture the rollover bit (declaring unsigned short instead of casting later)
    spaces = (spaces >> 1);                              // Rotate right 1
    rollover = (rollover >> 1);                          // scoot everything right in rollover
    rollover = rollover | (rollover_bit << (NUM_SPACES - 1)); // shift to the last spaces bit and | it into that place
    n--;
  }
  // Similar idea to above, & not needed becuase right-shifting spaces just drops bits off the number
  set_spaces(aisle, index, (spaces | rollover)); 
}

