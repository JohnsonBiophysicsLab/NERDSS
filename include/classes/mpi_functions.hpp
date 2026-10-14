/*! \file mpi_functions.hpp
 *
 */

#pragma once

#include <array>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#ifdef mpi_
#include "mpi.h"
#endif

// Serialize number into array of char called arrayRank
// starting from nArrayRank byte
// and increase nArrayRank by number of serialized bytes after:
#define PUSH(variable)                                            \
  *((__typeof__(variable) *)(arrayRank + nArrayRank)) = variable; \
  nArrayRank += sizeof(variable);
// TODO: instead of int, use

// Serialize number into toRank[iRank] array of char
// by storing received expression as variable into variable var,
// and appending toRank[iRank] with each byte of var:
#define PUSH_TO(variable, iRank)                                         \
  {                                                                      \
    __typeof__(variable) var = variable;                                 \
    int nBytes = sizeof(var);                                            \
    for (int nChr = 0; nChr < nBytes; nChr++)                            \
      mpiContext.toRank[iRank].push_back(((unsigned char *)&var)[nChr]); \
  }

// Deserialize from array of char called arrayRank starting from nArrayRank byte
// into number variable
// and increase nArrayRank by size of variable in bytes after:
#define POP(variable)                                             \
  variable = *((__typeof__(variable) *)(arrayRank + nArrayRank)); \
  nArrayRank += sizeof(variable);

// Serialize number into array of char
// starting from length byte
// and increase length by number of serialized bytes after:
#define PUSH_INTO(variable, array, length)                \
  *((__typeof__(variable) *)(array + length)) = variable; \
  length += sizeof(variable);
// TODO: instead of int, use

// Deserialize from array of char starting from length byte
// into number variable
// and increase length by size of variable in bytes after:
#define POP_FROM(variable, array, length)                 \
  variable = *((__typeof__(variable) *)(array + length)); \
  length += sizeof(variable);

void serialize_string(std::string to_serialize, unsigned char *arrayRank,
                      int &nArrayRank);

void deserialize_string(std::string &to_deserialize, unsigned char *arrayRank,
                        int &nArrayRank);

void serialize_vector_strings(std::vector<std::string> to_serialize,
                              unsigned char *arrayRank, int &nArrayRank);

void deserialize_vector_strings(std::vector<std::string> &to_deserialize,
                                unsigned char *arrayRank, int &nArrayRank);

template <typename T>
void serialize_primitive_vector(std::vector<T> to_serialize,
                                unsigned char *arrayRank, int &nArrayRank) {
  // serialize to_serialize vector of primitive type T
  *((int *)(arrayRank + nArrayRank)) =
      to_serialize
          .size();  // store vector size in first bytes needed for an int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (auto &it : to_serialize) {  // for all vector elements
    *((T *)&(arrayRank[nArrayRank])) =
        it;                   // read one primitive type T, e.g. int
    nArrayRank += sizeof(T);  // increase number of bytes read
  }
}

template <typename T>
void deserialize_primitive_vector(std::vector<T> &to_deserialize,
                                  unsigned char *arrayRank, int &nArrayRank) {
  to_deserialize.clear();
  int vector_size = *(
      (int *)(arrayRank +
              nArrayRank));  // extract number of vector elements
                             // from first bytes of arrayRank needed for one int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (int i = 0; i < vector_size; i++) {  // for all vector elements
    T temp = *((T *)&(
        arrayRank[nArrayRank]));     // extract an element of type T (e.g. int)
    nArrayRank += sizeof(T);         // increase number of bytes read
    to_deserialize.push_back(temp);  // push element to the back of the vector
  }
}

template <typename T>
void serialize_abstract_vector(std::vector<T> to_serialize,
                               unsigned char *arrayRank, int &nArrayRank,
                               bool verbose = false) {
  // serialize to_serialize vector of abstract type T
  *((int *)(arrayRank + nArrayRank)) =
      to_serialize
          .size();  // store vector size in first bytes needed for an int
  if (verbose)
    std::cout << "Serialize: to_serialize.size() = " << to_serialize.size()
              << std::endl;
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (auto &it : to_serialize)  // for all vector elements
    it.serialize(arrayRank,
                 nArrayRank);  // read one abstract type T, e.g. Vec3D
                               // and increase number of bytes read
                               // Note: T has to implement serialize method
  if (verbose)
    std::cout << "Serialize: nArrayRank = " << nArrayRank << std::endl;
}

template <typename T>
void deserialize_abstract_vector(std::vector<T> &to_deserialize,
                                 unsigned char *arrayRank, int &nArrayRank,
                                 bool verbose = false) {
  // Note: if -std=c++14 or -std=gnu++14 is used, typename T would not have to
  // be passed, but auto could be used
  to_deserialize.clear();
  int vector_size = *(
      (int *)(arrayRank +
              nArrayRank));  // extract number of vector elements
                             // from first bytes of arrayRank needed for one int
  if (verbose)
    std::cout << "Deserialize: vector size = " << vector_size << std::endl;
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (int i = 0; i < vector_size; i++) {  // for all vector elements
    T temp;  // Create local object of type T (e.g. Vec3D)
    temp.deserialize(
        arrayRank, nArrayRank);  // extract one vector element (e.g. Vec3D)
                                 // and increase number of bytes read
                                 // Note: T has to implement deserialize method
    to_deserialize.push_back(temp);  // push element to the back of the vector
    // if(verbose) std::cout << "push_back finished." << std::endl;
  }
  if (verbose)
    std::cout << "Deserialize: nArrayRank = " << nArrayRank << std::endl;
}

template <typename T>
void serialize_primitive_matrix(std::vector<std::vector<T> > to_serialize,
                                unsigned char *arrayRank, int &nArrayRank) {
  // serialize to_serialize matrix of abstract type T
  *((int *)(arrayRank + nArrayRank)) =
      to_serialize
          .size();  // store number or rows in first bytes needed for an int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (auto &it : to_serialize)  // for all rows
    serialize_primitive_vector<T>(
        it, arrayRank,
        nArrayRank);  // read abstract vector T
                      // and increase number of bytes read
                      // Note: T has to implement serialize method
}

template <typename T>
void deserialize_primitive_matrix(std::vector<std::vector<T> > &to_deserialize,
                                  unsigned char *arrayRank, int &nArrayRank) {
  // Note: if -std=c++14 or -std=gnu++14 is used, typename T would not have to
  // be passed, but auto could be used
  to_deserialize.clear();
  int vector_size = *(
      (int *)(arrayRank +
              nArrayRank));  // extract number of rows
                             // from first bytes of arrayRank needed for one int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (int i = 0; i < vector_size; i++) {  // for all rows
    std::vector<T> temp;                   // Create local vector of type T
    deserialize_primitive_vector<T>(
        temp, arrayRank,
        nArrayRank);  // extract one vector
                      // and increase number of bytes read
                      // Note: T has to implement deserialize method
    to_deserialize.push_back(
        temp);  // push resulting vector to the back of the matrix
  }
}

template <typename T, std::size_t S>
void serialize_vector_array(std::vector<std::array<int, S> > to_serialize,
                            unsigned char *arrayRank, int &nArrayRank) {
  // serialize to_serialize matrix of abstract type T
  *((int *)(arrayRank + nArrayRank)) =
      to_serialize
          .size();  // store number or rows in first bytes needed for an int
  // std::cout << "to_serialize.size()=" << to_serialize.size() << std::endl;
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (auto &it : to_serialize) {  // for all rows
    for (size_t k = 0; k < S; k++) {
      *((T *)&(arrayRank[nArrayRank])) = it[k];
      nArrayRank += sizeof(T);
    }
  }
}

template <typename T, std::size_t S>
void deserialize_vector_array(std::vector<std::array<T, S> > &to_deserialize,
                              unsigned char *arrayRank, int &nArrayRank) {
  // Note: if -std=c++14 or -std=gnu++14 is used, typename T would not have to
  // be passed, but auto could be used
  to_deserialize.clear();
  int vector_size =
      *((int *)(arrayRank + nArrayRank));  // extract number of rows
  // std::cout << "vector_size=" << vector_size << std::endl;
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (int i = 0; i < vector_size; i++) {  // for all rows
    std::array<T, S> temp;                 // Create local vector of type T
    for (size_t k = 0; k < S; k++) {
      temp[k] = *((T *)&(
          arrayRank[nArrayRank]));  // extract an element of type T (e.g. int)
      nArrayRank += sizeof(T);      // increase number of bytes read
    }
    to_deserialize.push_back(
        temp);  // push resulting vector to the back of the matrix
  }
}

template <typename T>
void serialize_abstract_matrix(std::vector<std::vector<T> > to_serialize,
                               unsigned char *arrayRank, int &nArrayRank) {
  // serialize to_serialize matrix of abstract type T
  *((int *)(arrayRank + nArrayRank)) =
      to_serialize
          .size();  // store number or rows in first bytes needed for an int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (auto &it : to_serialize)  // for all rows
    serialize_abstract_vector<T>(
        it, arrayRank,
        nArrayRank);  // read abstract vector T
                      // and increase number of bytes read
                      // Note: T has to implement serialize method
}

template <typename T>
void deserialize_abstract_matrix(std::vector<std::vector<T> > &to_deserialize,
                                 unsigned char *arrayRank, int &nArrayRank) {
  // Note: if -std=c++14 or -std=gnu++14 is used, typename T would not have to
  // be passed, but auto could be used
  to_deserialize.clear();
  int vector_size = *(
      (int *)(arrayRank +
              nArrayRank));  // extract number of rows
                             // from first bytes of arrayRank needed for one int
  nArrayRank +=
      sizeof(int);  // set current number of bytes read to the size of int
  for (int i = 0; i < vector_size; i++) {  // for all rows
    std::vector<T> temp;                   // Create local vector of type T
    deserialize_abstract_vector<T>(
        temp, arrayRank,
        nArrayRank);  // extract one vector
                      // and increase number of bytes read
                      // Note: T has to implement deserialize method
    to_deserialize.push_back(
        temp);  // push resulting vector to the back of the matrix
  }
}

#ifdef mpi_

template <typename T>
bool test_object_serialization(T to_test, unsigned char *array1,
                               bool verbose = false) {
  if (verbose)
    std::cout << "--------------------- Testing serialization and "
                 "deserialization ---------------------"
              << std::endl;

  int sizeSeserialized = 0;
  to_test.serialize(array1, sizeSeserialized);
  T deserialized_object;
  int sizeDeserialized = 0;
  deserialized_object.deserialize(array1, sizeDeserialized);
  unsigned char *array2 = (unsigned char *)malloc(sizeSeserialized);
  int sizeSeserializedAgain = 0;
  deserialized_object.serialize(array2, sizeSeserializedAgain);

  if (verbose)
    std::cout << "Tseserialized vector size = " << sizeSeserialized << " bytes."
              << std::endl;
  if (verbose)
    std::cout << "Tdeseserialized vector size = " << sizeDeserialized
              << " bytes." << std::endl;
  if (sizeSeserialized != sizeDeserialized) {
    std::cout
        << "Test abstract vector serialization failed. Sizes do not match."
        << std::endl;
    free(array2);
    return false;  // not equal as they are supposed to be
  }
  if (sizeDeserialized != sizeSeserializedAgain) {
    std::cout << "Test abstract vector serialization failed. Serialization "
                 "size differs from serialization of deserialized vector."
              << std::endl;
    free(array2);
    return false;  // not equal as they are supposed to be
  }

  for (int i = 0; i < sizeSeserialized; i++) {  // for all vector elements
    if (array1[i] != array2[i]) {
      if (verbose)
        std::cout << "Test abstract vector serialization failed. Serialized "
                     "bytes on position "
                  << i << " do not match." << std::endl;
      free(array2);
      return false;
    }
  }

  free(array2);
  return true;  // serialized arrayRank of bytes is equal to deserialized and
                // serialized again.
}

template <typename T>
bool test_abstract_vector_serialization(std::vector<T> to_test,
                                        unsigned char *array1,
                                        bool verbose = false) {
  if (verbose)
    std::cout << "--------------------- Testing serialization and "
                 "deserialization ---------------------"
              << std::endl;

  int sizeSeserialized = 0;
  serialize_abstract_vector<T>(to_test, array1, sizeSeserialized, verbose);
  std::vector<T> deserialized_vector;
  int sizeDeserialized = 0;
  deserialize_abstract_vector<T>(deserialized_vector, array1, sizeDeserialized,
                                 verbose);
  unsigned char *array2 = (unsigned char *)malloc(sizeSeserialized);
  if (!array2) {
    std::cout << "Memory allocation failed." << std::endl;
    exit(1);
  }
  int sizeSeserializedAgain = 0;
  serialize_abstract_vector<T>(deserialized_vector, array2,
                               sizeSeserializedAgain, verbose);

  if (verbose)
    std::cout << "Tseserialized vector size = " << sizeSeserialized << " bytes."
              << std::endl;
  if (verbose)
    std::cout << "Tdeseserialized vector size = " << sizeDeserialized
              << " bytes." << std::endl;
  if (sizeSeserialized != sizeDeserialized) {
    std::cout
        << "Test abstract vector serialization failed. Sizes do not match."
        << std::endl;
    free(array2);
    return false;  // not equal as they are supposed to be
  }
  if (sizeDeserialized != sizeSeserializedAgain) {
    std::cout << "Test abstract vector serialization failed. Serialization "
                 "size differs from serialization of deserialized vector."
              << std::endl;
    free(array2);
    return false;  // not equal as they are supposed to be
  }

  for (int i = 0; i < sizeSeserialized; i++) {  // for all vector elements
    if (array1[i] != array2[i]) {
      if (verbose)
        std::cout << "Test abstract vector serialization failed. Serialized "
                     "bytes on position "
                  << i << " do not match." << std::endl;
      free(array2);
      return false;
    }
  }

  free(array2);
  return true;  // serialized arrayRank of bytes is equal to deserialized and
                // serialized again.
}
#endif
