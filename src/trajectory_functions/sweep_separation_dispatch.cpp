/*! \file sweep_separation_dispatch.cpp
 * \brief Picks the box or sphere overlap sweep for each kind of complex.
 *
 * These four functions were four files, each holding a dozen live lines behind
 * two hundred or so lines of its own former body, commented out.  Nothing is
 * left in any of them but the geometry choice, and the choice is the same in
 * all four - which is only visible with them side by side.  The fiber sweep is
 * the odd one: there is no spherical version of it, so it forwards
 * unconditionally.
 *
 * The commented-out bodies are gone; `git log` has them.
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

//! \brief Overlap sweep for a complex in solution.
void sweep_separation_complex_rot(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    if (membraneObject.isSphere)
        sweep_separation_complex_rot_sphere(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
    else
        sweep_separation_complex_rot_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
}

//! \brief Overlap sweep for a complex bound to the membrane.
void sweep_separation_complex_rot_memtest(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    if (membraneObject.isSphere)
        sweep_separation_complex_rot_memtest_sphere(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
    else
        sweep_separation_complex_rot_memtest_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
}

//! \brief Membrane sweep that resolves overlap a cluster at a time.
void sweep_separation_complex_rot_memtest_cluster(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    if (membraneObject.isSphere)
        sweep_separation_complex_rot_memtest_cluster_sphere(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
    else
        sweep_separation_complex_rot_memtest_cluster_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
}

//! \brief Overlap sweep for a complex on a fiber.  No spherical version exists.
void sweep_separation_complex_rot_fiber(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    sweep_separation_complex_rot_fiber_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
}
