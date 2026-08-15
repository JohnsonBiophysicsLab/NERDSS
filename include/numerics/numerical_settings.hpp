/*! \file numerical_settings.hpp
 * \brief User-configurable numerical policy and its defaults.
 */
#pragma once

#include "numerics/comparison.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

struct NumericalSettings {
    struct Integration {
        // GSL tolerances used to construct the 2D reaction tables.
        double tableAbsoluteError { 1e-7 };
        double tableRelativeError { 1e-7 };
        double fallbackError { 1e-6 };
        double tailCutoff { 1e-10 };

        // The normalization table historically used a separate criterion.
        double normalizationAbsoluteError { 1e-6 };
        double normalizationRelativeError { 1e-6 };
    } integration;

    struct TableLookup {
        ComparisonTolerance reactionRate { 1e-8, 0.0 };
        ComparisonTolerance diffusionCoefficient { 1e-4, 0.0 };
    } tableLookup;

    struct Classification {
        // These remain separate because the existing explicit- and
        // implicit-lipid algorithms deliberately use different cutoffs.
        double explicitLipidFlatDiffusion { 1e-10 };
        double implicitLipidFlatDiffusion { 1e-15 };
    } classification;

    struct AssociationAngles {
        // General equality for periodic/end-point association angles.
        ComparisonTolerance sameAngle { 1e-4, 0.0 };

        // A tighter criterion for deciding whether a rotation has converged.
        double rotationConvergenceTolerance { 1e-8 };

        // Prevents numerical noise at 0 and pi from changing a dihedral sign.
        double endpointSignTolerance { 1e-11 };
    } associationAngles;

    void validate() const
    {
        require_positive("numericsIntegrationAbsError", integration.tableAbsoluteError);
        require_positive("numericsIntegrationRelError", integration.tableRelativeError);
        require_positive("numericsIntegrationFallbackError", integration.fallbackError);
        require_positive("numericsIntegrationTailCutoff", integration.tailCutoff);
        require_positive("numericsNormalizationAbsError", integration.normalizationAbsoluteError);
        require_positive("numericsNormalizationRelError", integration.normalizationRelativeError);

        require_nonnegative("numericsTableRateAbsTolerance", tableLookup.reactionRate.absolute);
        require_nonnegative("numericsTableRateRelTolerance", tableLookup.reactionRate.relative);
        require_nonnegative("numericsTableDiffusionAbsTolerance", tableLookup.diffusionCoefficient.absolute);
        require_nonnegative("numericsTableDiffusionRelTolerance", tableLookup.diffusionCoefficient.relative);
        if (tableLookup.reactionRate.absolute == 0.0 && tableLookup.reactionRate.relative == 0.0)
            throw std::invalid_argument("Table-rate absolute and relative tolerances cannot both be zero.");
        if (tableLookup.diffusionCoefficient.absolute == 0.0
            && tableLookup.diffusionCoefficient.relative == 0.0)
            throw std::invalid_argument("Table-diffusion absolute and relative tolerances cannot both be zero.");

        require_nonnegative(
            "numericsExplicitLipidFlatDiffusion", classification.explicitLipidFlatDiffusion);
        require_nonnegative(
            "numericsImplicitLipidFlatDiffusion", classification.implicitLipidFlatDiffusion);

        require_nonnegative("numericsAssociationSameAngleAbsTolerance", associationAngles.sameAngle.absolute);
        require_nonnegative("numericsAssociationSameAngleRelTolerance", associationAngles.sameAngle.relative);
        if (associationAngles.sameAngle.absolute == 0.0 && associationAngles.sameAngle.relative == 0.0)
            throw std::invalid_argument(
                "Association-angle absolute and relative tolerances cannot both be zero.");
        require_positive(
            "numericsAssociationRotationTolerance", associationAngles.rotationConvergenceTolerance);
        require_positive(
            "numericsAssociationEndpointSignTolerance", associationAngles.endpointSignTolerance);
    }

private:
    static void require_positive(const char* name, double value)
    {
        if (!std::isfinite(value) || value <= 0.0)
            throw std::invalid_argument(std::string(name) + " must be finite and greater than zero.");
    }

    static void require_nonnegative(const char* name, double value)
    {
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument(std::string(name) + " must be finite and non-negative.");
    }
};
