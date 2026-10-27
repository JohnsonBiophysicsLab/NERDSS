/*! \file mpi_column_split_check.cpp
 * \brief Checks how nerdss_mpi divides the x columns of the box between ranks.
 *
 * MpiContext::column_share() must, for every column count and every rank
 * count up to max_ranks_for_columns():
 *   - give every rank at least one column, and the ranks' shares in order
 *     must tile [0, xTotal) exactly;
 *   - give every rank between two others at least minMiddleColumns, because
 *     with fewer its first half-step reaches the column its right neighbour
 *     binds into during the same half-step (see region_mid_index());
 *   - leave the even split it always made untouched wherever that split
 *     already satisfies the rule, so that np>1 runs that were valid before
 *     keep their decomposition.
 * And one rank more than max_ranks_for_columns() must be impossible, i.e.
 * the bound is tight.
 *
 * Built and run by `make serial && make checks`.  Exits non-zero on failure.
 */
#include "mpi_context.hpp"

#include <iostream>
#include <vector>

#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;

int main()
{
    int failures { 0 };
    auto fail = [&](int xTotal, int nprocs, const char* what) {
        if (failures < 20)
            std::cerr << "FAIL: " << xTotal << " columns, " << nprocs << " ranks: " << what << '\n';
        ++failures;
    };

    for (int xTotal = 1; xTotal <= 80; xTotal++) {
        const int most = MpiContext::max_ranks_for_columns(xTotal);
        // tight: one more rank cannot satisfy ends >= 1 and middles >= 3
        const int more = most + 1;
        if (more <= xTotal && more >= 3 && xTotal - MpiContext::minMiddleColumns * (more - 2) >= 2)
            fail(xTotal, more, "max_ranks_for_columns() is not the largest valid rank count");
        for (int nprocs = 1; nprocs <= most; nprocs++) {
            MpiContext context(nprocs, 16);
            int next = 0;
            for (int rank = 0; rank < nprocs; rank++) {
                int first = -1, count = -1;
                context.column_share(xTotal, rank, first, count);
                if (first != next)
                    fail(xTotal, nprocs, "shares do not tile the columns in order");
                if (count < 1)
                    fail(xTotal, nprocs, "a rank owns no column");
                if (rank > 0 && rank < nprocs - 1 && count < MpiContext::minMiddleColumns)
                    fail(xTotal, nprocs, "a middle rank owns fewer than minMiddleColumns");
                // the split used before minMiddleColumns existed
                const int div = xTotal / nprocs, mod = xTotal % nprocs;
                const int oldFirst = div * rank + (rank < mod ? rank : mod);
                const int oldCount = div + (rank < mod ? 1 : 0);
                const bool oldValid = nprocs < 3 || div >= MpiContext::minMiddleColumns;
                if (oldValid && (first != oldFirst || count != oldCount))
                    fail(xTotal, nprocs, "a split that was already valid changed");
                int tempRank = rank;
                context.init_x_domain_and_offset(xTotal, tempRank);
                if (context.startCell != first || context.endCell != first + count - 1)
                    fail(xTotal, nprocs, "init_x_domain_and_offset() disagrees with column_share()");
                next = first + count;
            }
            if (next != xTotal)
                fail(xTotal, nprocs, "shares do not cover the last column");
        }
    }

    if (failures) {
        std::cerr << failures << " failure(s)\n";
        return 1;
    }
    std::cout << "mpi_column_split_check: ok\n";
    return 0;
}
