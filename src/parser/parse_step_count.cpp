/*! \file parse_step_count.cpp
 * \ingroup Parser
 * \brief Reads the step counts of the parameters block (nItr, timeWrite, ...) as exact whole numbers.
 *
 * These used to go through std::stod and a cast to long long, "otherwise scientific notation breaks".
 * That accepted 1000000.7 as 1000000, turned an nItr too long for a double into a rounded or
 * overflowed count, and let "nan" and "inf" through.  The parser here keeps the exponent forms
 * (1e6, 5E5, 2.5e6) and does the arithmetic on the digit string, so a value is either the whole
 * number the user wrote or rejected.
 */
#include "parser/parser_functions.hpp"

#include <cctype>
#include <climits>
#include <sstream>

StepCountStatus parse_step_count(const std::string& text, long long& value)
{
    const std::size_t n = text.size();
    std::size_t pos { 0 };

    bool negative { false };
    if (pos < n && (text[pos] == '+' || text[pos] == '-')) {
        negative = text[pos] == '-';
        ++pos;
    }

    // Mantissa: digits with at most one decimal point.  The point is dropped and remembered as a
    // shift, so 2.5e6 is the digit string "25" moved six places right and two places left.
    std::string digits;
    long long pointShift { 0 };
    bool seenPoint { false };
    for (; pos < n; ++pos) {
        const char c = text[pos];
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits += c;
            if (seenPoint)
                ++pointShift;
        } else if (c == '.' && !seenPoint) {
            seenPoint = true;
        } else {
            break;
        }
    }
    if (digits.empty())
        return StepCountStatus::notANumber; // "", ".", "-", "nan", "inf", "e5", ...

    // Exponent.  Its magnitude is capped: anything past 10^4 is already far outside a long long
    // (or, negative, far below 1), and the cap keeps a page of exponent digits from overflowing.
    long long exponent { 0 };
    if (pos < n && (text[pos] == 'e' || text[pos] == 'E')) {
        ++pos;
        bool exponentNegative { false };
        if (pos < n && (text[pos] == '+' || text[pos] == '-')) {
            exponentNegative = text[pos] == '-';
            ++pos;
        }
        std::size_t exponentDigits { 0 };
        for (; pos < n && std::isdigit(static_cast<unsigned char>(text[pos])); ++pos, ++exponentDigits) {
            if (exponent < 10000)
                exponent = exponent * 10 + (text[pos] - '0');
        }
        if (exponentDigits == 0)
            return StepCountStatus::notANumber; // "1e", "1e+"
        if (exponentNegative)
            exponent = -exponent;
    }
    if (pos != n)
        return StepCountStatus::notANumber; // trailing characters: "1000;", "12abc", "0x10"

    const std::size_t firstNonZero = digits.find_first_not_of('0');
    if (firstNonZero == std::string::npos) {
        value = 0; // 0, 0.0, 0e9, -0
        return StepCountStatus::ok;
    }
    digits.erase(0, firstNonZero);

    const long long shift = exponent - pointShift; // places the decimal point moves to the right
    if (shift < 0) {
        // The digits that would end up behind the point must all be zero: 1000000.0 and 1.0e3
        // are whole numbers, 1000.5 and 1.5 are not.
        const std::size_t drop = static_cast<std::size_t>(-shift);
        if (drop >= digits.size())
            return StepCountStatus::notWholeNumber; // 0 < |value| < 1
        if (digits.find_first_not_of('0', digits.size() - drop) != std::string::npos)
            return StepCountStatus::notWholeNumber;
        digits.erase(digits.size() - drop);
    } else if (shift > 0) {
        if (shift > 19)
            return StepCountStatus::outOfRange; // a non-zero digit times 10^20 is past 2^63
        digits.append(static_cast<std::size_t>(shift), '0');
    }
    if (digits.size() > 19)
        return StepCountStatus::outOfRange;

    // At most 19 digits, so the magnitude fits an unsigned long long (max 1.8e19).
    unsigned long long magnitude { 0 };
    for (char c : digits)
        magnitude = magnitude * 10 + static_cast<unsigned long long>(c - '0');

    const unsigned long long limit = negative ? static_cast<unsigned long long>(LLONG_MAX) + 1ULL
                                              : static_cast<unsigned long long>(LLONG_MAX);
    if (magnitude > limit)
        return StepCountStatus::outOfRange;
    if (negative)
        value = magnitude == limit ? LLONG_MIN : -static_cast<long long>(magnitude);
    else
        value = static_cast<long long>(magnitude);
    return StepCountStatus::ok;
}

// What each step-count keyword must be, and what stands in for a value it cannot be.  The
// meaning of 0 or below is the rule every periodic write follows (`interval > 0 && simItr %
// interval == 0` in EXEs/nerdss.cpp): off.  checkPoint and transitionWrite read -1 as their
// documented default of nItr / 10 instead.  A null offMeaning means a non-positive value is a
// mistake rather than "off".
const StepCountRule StepCountRules::nItr { "nItr", 0, true, nullptr, nullptr };
const StepCountRule StepCountRules::timeWrite { "timeWrite", 10, false,
    "turns the periodic statistics output off; only the first and the last rows are written", nullptr };
const StepCountRule StepCountRules::trajWrite { "trajWrite", 10, false,
    "turns the trajectory write off; no frame is written, not even the initial one", nullptr };
const StepCountRule StepCountRules::restartWrite { "restartWrite", 10, false,
    "turns the periodic restart write off; the restart file at the end of the run is still written", nullptr };
const StepCountRule StepCountRules::pdbWrite { "pdbWrite", -1, false,
    "turns the PDB write off; no PDB file is written, not even the initial one", nullptr };
const StepCountRule StepCountRules::checkPoint { "checkPoint", -1, false,
    "turns checkpoint restart files off", "means a checkpoint restart file every nItr / 10 timesteps" };
const StepCountRule StepCountRules::transitionWrite { "transitionWrite", -1, false,
    "turns the periodic transition-matrix write off", "means the transition matrix every nItr / 10 timesteps" };
const StepCountRule StepCountRules::bondedComplexWrite { "bondedComplexWrite", -1, false,
    "turns the bonded-complex JSON files off", nullptr };

StepCountReading resolve_step_count(const std::string& text, const StepCountRule& rule)
{
    StepCountReading reading { rule.fallback, false, "" };
    std::ostringstream message;

    long long value { 0 };
    const StepCountStatus status = parse_step_count(text, value);

    if (status == StepCountStatus::ok && value > 0) {
        reading.value = value;
        return reading;
    }

    // 0 or below is documented for every write interval and stored as written; the note says what
    // it does, since users set -1 on one keyword expecting what it means on another.
    if (status == StepCountStatus::ok && rule.offMeaning != nullptr) {
        reading.value = value;
        const char* meaning = (value == -1 && rule.minusOneMeaning != nullptr) ? rule.minusOneMeaning
                                                                               : rule.offMeaning;
        message << "Note: " << rule.keyword << " = " << value << " " << meaning << ".";
        reading.message = message.str();
        return reading;
    }

    message << (rule.fatal ? "Error: " : "Warning: ") << rule.keyword << " = \"" << text << "\" ";
    switch (status) {
    case StepCountStatus::ok:
        message << "must be a positive number of timesteps";
        break;
    case StepCountStatus::notWholeNumber:
        message << "is not a whole number of timesteps";
        break;
    case StepCountStatus::outOfRange:
        message << "is outside the range of a 64-bit integer";
        break;
    case StepCountStatus::notANumber:
        message << "is not a number (accepted forms: 100000, 1e5, 2.5E5)";
        break;
    }
    if (rule.fatal) {
        message << ", and " << rule.keyword << " has no default. Exiting.";
        reading.fatal = true;
    } else {
        message << "; using the default, " << rule.keyword << " = " << rule.fallback << ".";
    }
    reading.message = message.str();
    return reading;
}
