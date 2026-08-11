/*! \file attributes.hpp
 * \brief Portable spellings of the compiler attributes used in the headers.
 */

#pragma once

/*!
 * \brief Warn when the result of a pure, read-only function is discarded.
 *
 * `[[nodiscard]]` is C++17, and this project is built as C++11 (`-std=c++0x` in
 * the Makefile, `CMAKE_CXX_STANDARD 11` in CMakeLists.txt).  Spelling it
 * `[[nodiscard]]` anyway makes clang warn about a C++17 extension at every
 * include, so the GNU attribute is used instead until the standard is bumped;
 * both compilers this project is built with support it.
 */
#ifndef NERDSS_NODISCARD
#if defined(__cplusplus) && __cplusplus >= 201703L
#define NERDSS_NODISCARD [[nodiscard]]
#elif defined(__GNUC__) || defined(__clang__)
#define NERDSS_NODISCARD __attribute__((warn_unused_result))
#else
#define NERDSS_NODISCARD
#endif
#endif
