#include "reactions/association/association.hpp"

Vec3D create_arbitrary_vector(Vec3D& vec)
{
    vec.normalize();

    /* The x-axis branch this function is documented as having is unreachable,
     * and always has been.  It was guarded by
     *
     *     vec.dot_theta(x_axis) != 0 && vec.dot_theta(x_axis) != M_PI
     *
     * where `x_axis` was a fresh `Vector(1, 0, 0)` whose `magnitude` cache was
     * never written.  `dot_theta` returned 0 for any operand with a cached
     * magnitude below 1E-8, so the first comparison was `0 != 0` on every call
     * and the y-axis arm always ran - while printing a spurious "angle between
     * vectors with at least one of magnitude 0" warning on the way past.
     *
     * Reinstating the x-axis branch would change which rotation axis
     * theta_rotation() picks whenever sigma and the com-iface vector are
     * parallel, so the arm that actually ran is what is kept here.  What is
     * lost is the warning.
     */
    return vec.unit_cross(Vec3D { 0, 1, 0 });
}
