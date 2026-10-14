/*! \file step_count_check.cpp
 * \brief Checks parse_step_count() and resolve_step_count(), the readers for nItr, timeWrite and
 * the other step counts of the parameters block.
 *
 * Built and run by `make serial && make checks`.  Exits non-zero on the first failure.
 */
#include "parser/parser_functions.hpp"

#include <climits>
#include <cstdlib>
#include <iostream>

#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;

static int failures { 0 };

static void expect_value(const std::string& text, long long expected)
{
    long long value { -12345 };
    const StepCountStatus status = parse_step_count(text, value);
    if (status != StepCountStatus::ok || value != expected) {
        std::cerr << "FAIL: \"" << text << "\" should read as " << expected << " (status "
                  << static_cast<int>(status) << ", value " << value << ")\n";
        ++failures;
    }
}

static void expect_status(const std::string& text, StepCountStatus expected)
{
    long long value { 0 };
    const StepCountStatus status = parse_step_count(text, value);
    if (status != expected) {
        std::cerr << "FAIL: \"" << text << "\" should give status " << static_cast<int>(expected) << ", got "
                  << static_cast<int>(status) << '\n';
        ++failures;
    }
}

static void expect_reading(const std::string& text, const StepCountRule& rule, long long value, bool fatal,
    bool hasMessage)
{
    const StepCountReading reading = resolve_step_count(text, rule);
    if (reading.value != value || reading.fatal != fatal || reading.message.empty() == hasMessage) {
        std::cerr << "FAIL: " << rule.keyword << " = \"" << text << "\" -> value " << reading.value << ", fatal "
                  << reading.fatal << ", message \"" << reading.message << "\"\n";
        ++failures;
    }
}

int main()
{
    // Plain integers, the forms the sample inputs use.
    expect_value("100000", 100000);
    expect_value("500000000", 500000000);
    expect_value("10000000000", 10000000000LL);
    expect_value("+10", 10);
    expect_value("-1", -1);
    expect_value("0", 0);
    expect_value("-0", 0);
    expect_value("007", 7);

    // Exponent forms stay supported.
    expect_value("1e6", 1000000);
    expect_value("1E6", 1000000);
    expect_value("5E5", 500000);
    expect_value("6e5", 600000);
    expect_value("1e9", 1000000000);
    expect_value("1e+3", 1000);
    expect_value("2.5e6", 2500000);
    expect_value("2.50e6", 2500000);
    expect_value("0.5e1", 5);
    expect_value(".5e1", 5);
    expect_value("1.e3", 1000);
    expect_value("12300e-2", 123);
    expect_value("1e18", 1000000000000000000LL);
    expect_value("0e99999", 0);
    expect_value("0.000e-99999", 0);

    // A decimal point that writes a whole number is fine.
    expect_value("1000000.0", 1000000);
    expect_value("1000000.000", 1000000);
    expect_value("1.0e3", 1000);

    // A value that is too long for a double to hold exactly is read exactly.
    expect_value("9007199254740993", 9007199254740993LL); // 2^53 + 1
    expect_value("1234567890123456789", 1234567890123456789LL);
    expect_value("9223372036854775807", LLONG_MAX);
    expect_value("-9223372036854775808", LLONG_MIN);

    // Decimals are rejected, not truncated.
    expect_status("1000.5", StepCountStatus::notWholeNumber);
    expect_status("1000000.7", StepCountStatus::notWholeNumber);
    expect_status("1.5", StepCountStatus::notWholeNumber);
    expect_status("0.5", StepCountStatus::notWholeNumber);
    expect_status("1e-1", StepCountStatus::notWholeNumber);
    expect_status("1.5e0", StepCountStatus::notWholeNumber);
    expect_status("12345e-3", StepCountStatus::notWholeNumber);
    expect_status("1e-99999", StepCountStatus::notWholeNumber);

    // Not numbers: nan, inf, empty, junk, incomplete exponents, hex.
    expect_status("nan", StepCountStatus::notANumber);
    expect_status("NaN", StepCountStatus::notANumber);
    expect_status("inf", StepCountStatus::notANumber);
    expect_status("-inf", StepCountStatus::notANumber);
    expect_status("INFINITY", StepCountStatus::notANumber);
    expect_status("", StepCountStatus::notANumber);
    expect_status("-", StepCountStatus::notANumber);
    expect_status(".", StepCountStatus::notANumber);
    expect_status("e5", StepCountStatus::notANumber);
    expect_status("1e", StepCountStatus::notANumber);
    expect_status("1e+", StepCountStatus::notANumber);
    expect_status("1000;", StepCountStatus::notANumber);
    expect_status("12abc", StepCountStatus::notANumber);
    expect_status("0x10", StepCountStatus::notANumber);
    expect_status("1.2.3", StepCountStatus::notANumber);
    expect_status("1,000", StepCountStatus::notANumber);
    expect_status("1 000", StepCountStatus::notANumber);
    expect_status("--1", StepCountStatus::notANumber);

    // Out of range for a long long, rather than overflowing.
    expect_status("9223372036854775808", StepCountStatus::outOfRange);
    expect_status("-9223372036854775809", StepCountStatus::outOfRange);
    expect_status("12345678901234567890", StepCountStatus::outOfRange);
    expect_status("1e19", StepCountStatus::outOfRange);
    expect_status("1e30", StepCountStatus::outOfRange);
    expect_status("1e400", StepCountStatus::outOfRange);
    expect_status("1e99999999999999999999", StepCountStatus::outOfRange);

    // nItr: positive whole numbers only, nothing to fall back on.
    expect_reading("1e6", StepCountRules::nItr, 1000000, false, false);
    expect_reading("0", StepCountRules::nItr, 0, true, true);
    expect_reading("-5", StepCountRules::nItr, 0, true, true);
    expect_reading("1000.5", StepCountRules::nItr, 0, true, true);
    expect_reading("nan", StepCountRules::nItr, 0, true, true);
    expect_reading("inf", StepCountRules::nItr, 0, true, true);
    expect_reading("1e30", StepCountRules::nItr, 0, true, true);

    // Write intervals: a positive value is a period, 0 or below is off and passes as written with a
    // note; -1 on checkPoint and transitionWrite means nItr / 10.  An unreadable value falls back
    // to the default with a warning.
    expect_reading("1000", StepCountRules::timeWrite, 1000, false, false);
    expect_reading("0", StepCountRules::timeWrite, 0, false, true);
    expect_reading("-1", StepCountRules::timeWrite, -1, false, true);
    expect_reading("2.5", StepCountRules::timeWrite, 10, false, true);
    expect_reading("abc", StepCountRules::timeWrite, 10, false, true);
    expect_reading("5000", StepCountRules::trajWrite, 5000, false, false);
    expect_reading("-1", StepCountRules::trajWrite, -1, false, true);
    expect_reading("0", StepCountRules::trajWrite, 0, false, true);
    expect_reading("-7", StepCountRules::trajWrite, -7, false, true);
    expect_reading("1.5", StepCountRules::trajWrite, 10, false, true);
    expect_reading("0", StepCountRules::restartWrite, 0, false, true);
    expect_reading("inf", StepCountRules::restartWrite, 10, false, true);
    expect_reading("-1", StepCountRules::pdbWrite, -1, false, true);
    expect_reading("0", StepCountRules::pdbWrite, 0, false, true);
    expect_reading("1e400", StepCountRules::pdbWrite, -1, false, true);
    expect_reading("-1", StepCountRules::checkPoint, -1, false, true);
    expect_reading("-2", StepCountRules::checkPoint, -2, false, true);
    expect_reading("0", StepCountRules::checkPoint, 0, false, true);
    expect_reading("100.0", StepCountRules::checkPoint, 100, false, false);
    expect_reading("0", StepCountRules::transitionWrite, 0, false, true);
    expect_reading("nan", StepCountRules::transitionWrite, -1, false, true);
    expect_reading("0", StepCountRules::bondedComplexWrite, 0, false, true);
    expect_reading("1e3", StepCountRules::bondedComplexWrite, 1000, false, false);

    // The -1 note on checkPoint names nItr / 10, the 0 note does not.
    if (resolve_step_count("-1", StepCountRules::checkPoint).message.find("nItr / 10") == std::string::npos
        || resolve_step_count("0", StepCountRules::checkPoint).message.find("nItr / 10") != std::string::npos) {
        std::cerr << "FAIL: checkPoint notes for -1 and 0 should differ\n";
        ++failures;
    }

    if (failures != 0) {
        std::cerr << failures << " step-count check(s) failed\n";
        return EXIT_FAILURE;
    }
    std::cout << "step_count_check: all checks passed\n";
    return EXIT_SUCCESS;
}
