#include "math/Faddeeva.hpp"
#include "tracing.hpp"

double pirr_pfree_ratio_psF_1D(
    double rCurr, double r0, double tCurr, double Dtot, double bindrad, double ka, double ps_prev)
{
    if (bindrad == 0.0){
        /*
        This follows the solution of "semi-permeable boundary condition" for co-localized proteins.
        The PDE should be rewritten as:
            dp/dt = D d2p/dx2 - kact delta(x) p(x,t)
        with initial condition 
            p(x,0) = delta(x - r0)
        The solution is
            p(x,t) = (1/sqrt(4 pi D t)) [ exp(-(x - r0)^2/4 D t) ]
                 - (ka/2D) exp(ka (|x| + x0)/D + (ka^2 t)/D) erfc( (|x| + x0)/sqrt(4 D t) + ka sqrt(t/D) )
        The normalized p_free is
            p_free_norm = (1/sqrt(4 pi D t)) [ exp(-(x - r0)^2/4 D t) ]
        which is simply the free diffusion solution with initial condition p(x,0) = delta(x - r0).
        */
        // p_irr
        const double fDt = 4.0 * Dtot * tCurr;
        const double sq_fDt = sqrt(fDt);
        const double sq_pifDt = sqrt(4.0 * M_PI * Dtot * tCurr);
        const double dx_m = rCurr - r0;
        const double dx_p = rCurr + r0;
        const double sq_scale_t{ka * sqrt(tCurr / Dtot)};  
        double ep1 = exp(-dx_m * dx_m / fDt);
        double ep2 = exp(-dx_p * dx_p / fDt);
        double p_irr = (ep1) / sq_pifDt - ka / Dtot * ep2 * Faddeeva::erfcx( dx_p / sq_fDt + sq_scale_t );
        // normalized p_free
        double p_free_norm = (ep1) / sq_pifDt; 
        // reweighting ratio
        double ratio = p_irr / ps_prev / p_free_norm;
        return ratio;
    } else {
        /*
        This follows the solution of "radiation boundary condition"
        The PDE should be rewritten as:
            dp/dt = D d2p/dx2
        with initial condition 
            p(x,0) = delta(x - r0)
        and boundary condition at x = bindrad
            D dp/dx |(x=bindrad) = ka p(bindrad, t)
        */
        // combined variables
        const double fDt = 4.0 * Dtot * tCurr;
        const double sq_fDt = sqrt(fDt);
        const double sq_pifDt = sqrt(4.0 * M_PI * Dtot * tCurr);
        const double dx_m = rCurr - r0;
        const double dx_p = rCurr + r0;
        const double dx_sigma = rCurr + r0 - 2.0 * bindrad;
        const double x_sigma_m = r0 - bindrad;
        const double x_sigma_p = r0 + bindrad;
        const double sq_scale_t{ka * sqrt(tCurr / Dtot)};

        // normalized p_free
        double ep1 = exp(-dx_m * dx_m / fDt);
        double ep2 = exp(-dx_p * dx_p / fDt);
        double I_free = sqrt(M_PI * Dtot * tCurr) *
                (erf(x_sigma_m / sq_fDt) + erf(x_sigma_p / sq_fDt));
        double p_free_norm = (ep1 - ep2) / I_free;

        // p_irr
        double ep3 = exp(-dx_sigma * dx_sigma / fDt);
        double sep = dx_sigma / sq_fDt + sq_scale_t;
        double p_irr = (ep1 + ep3) / sq_pifDt - ka / Dtot * ep3 * Faddeeva::erfcx(sep);

        // reweighting ratio
        double ratio = p_irr / ps_prev / p_free_norm;

        return ratio;
    }
    
}
