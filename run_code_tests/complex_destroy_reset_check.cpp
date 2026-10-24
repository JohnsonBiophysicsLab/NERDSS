/*! \file complex_destroy_reset_check.cpp
 * \brief Checks that Complex::destroy() leaves nothing in the slot that the
 *        next occupant would inherit as its own ownership.
 *
 * A destroyed complex's slot goes on Complex::emptyComList, and
 * break_interaction() builds a dissociation product in it field by field.
 * nerdss_mpi decides who owns a new complex in claim_unowned_complexes(),
 * which only claims ownerRank == -1, so a slot that kept the destroyed
 * complex's owner handed the product to that rank, whoever made it.  No
 * model reaches it today (remove_empty_slots() empties the list at the end of
 * every MPI step, before the next dissociation), which is why this is a unit
 * check and not a run.
 *
 * Built and run by `make serial && make checks`.  Exits non-zero on failure.
 */
#include "classes/class_Molecule_Complex.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;

int main()
{
    int failures { 0 };
    std::vector<Molecule> moleculeList {};
    std::vector<Complex> complexList(2);
    for (int i = 0; i < 2; i++) {
        complexList[i].index = i;
        complexList[i].id = 40 + i;
        complexList[i].isEmpty = false;
        complexList[i].ownerRank = 3;
    }
    Complex::emptyComList.clear();

    complexList[1].destroy(moleculeList, complexList);

    if (complexList[1].ownerRank != -1) {
        std::cerr << "FAIL: a destroyed complex keeps ownerRank "
                  << complexList[1].ownerRank << ", which the next complex built in "
                  << "its slot inherits\n";
        ++failures;
    }
    if (!complexList[1].isEmpty || complexList[1].id != -1) {
        std::cerr << "FAIL: a destroyed complex is not marked empty\n";
        ++failures;
    }
    if (std::find(Complex::emptyComList.begin(), Complex::emptyComList.end(), 1)
        == Complex::emptyComList.end()) {
        std::cerr << "FAIL: the destroyed slot is not on emptyComList\n";
        ++failures;
    }
    if (complexList[0].ownerRank != 3) {
        std::cerr << "FAIL: destroying one complex changed another's owner\n";
        ++failures;
    }

    if (failures)
        return 1;
    std::cout << "complex_destroy_reset_check: ok\n";
    return 0;
}
