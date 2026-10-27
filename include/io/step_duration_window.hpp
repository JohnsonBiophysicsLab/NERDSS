#pragma once

#include <chrono>
#include <cstddef>
#include <vector>

/*! \ingroup IO
 * \brief The last few timestep durations, for the "Avg timestep duration" estimate.
 *
 * The step loop used to keep these in a std::vector and drop the oldest with
 * erase(begin()) every step: a shift of the whole window (8 KB for 1000
 * entries) per step, to feed a line printed once per timeWrite.  On Rockfish's
 * x86 nodes that overlapping memmove was 39% of clathrin_pioneer ffl and 7-13%
 * of many other fast-step models.  A ring overwrites the oldest entry instead,
 * and sum() adds the entries oldest first, the order std::accumulate() walked
 * the shifted vector in.
 */
class StepDurationWindow {
public:
    using duration = std::chrono::duration<double>;

    StepDurationWindow(std::size_t size, duration initial)
        : values(size, initial)
    {
    }

    //! Replaces the oldest duration with this one.
    void push(duration newest)
    {
        if (values.empty())
            return;
        values[oldest] = newest;
        if (++oldest == values.size())
            oldest = 0;
    }

    //! The sum of the window, oldest first.
    duration sum() const
    {
        duration total { 0 };
        for (std::size_t i { oldest }; i < values.size(); ++i)
            total += values[i];
        for (std::size_t i { 0 }; i < oldest; ++i)
            total += values[i];
        return total;
    }

    std::size_t size() const { return values.size(); }

private:
    std::vector<duration> values;
    std::size_t oldest { 0 }; //!< index of the oldest entry, the next to be replaced
};
