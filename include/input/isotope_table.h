#ifndef __CCEX_INPUT_ISOTOPE_TABLE_H_
#define __CCEX_INPUT_ISOTOPE_TABLE_H_

/**
 * @file isotope_table.h
 * @brief The nuclear isotopes CCEX knows without a gyrofile (src/input/isotope_table.cpp, generated from the EasySpin
 *        isotope database by example/tools/IsotopeTable/make_isotope_table.py).
 */

/** One isotope: name (mass number + symbol, e.g. "13C"), spin, gyromagnetic ratio (rad/ms/G). */
struct IsotopeEntry {
    const char* name;
    float spin;
    double gyro;
};

/** The isotope of this name, or nullptr when the table does not have it. */
const IsotopeEntry* findIsotope(const char* name);

/** Number of isotopes in the table. */
int isotopeTableSize();

#endif // __CCEX_INPUT_ISOTOPE_TABLE_H_
