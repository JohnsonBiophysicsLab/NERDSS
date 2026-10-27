#pragma once

#include <math.h>
#ifdef mpi_
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "classes/mpi_functions.hpp"

#define min2(a, b) ((a) < (b) ? (a) : (b))

struct Membrane;
struct SimulVolume;
struct Molecule;
struct Complex;

// The members below are guarded by mpi_, so this struct has one layout in a
// serial build and another in an MPI build: 176 bytes versus 272, with xOffset
// at byte 56 versus 136.  Mixing objects from the two modes in one link used to
// succeed and then misread every domain bound at run time.
//
// The inline namespace puts the build mode into the mangled name of everything
// that mentions MpiContext, so such a link now fails with an undefined symbol
// naming the mode it wanted.  It is inline, so `MpiContext` still resolves
// unqualified everywhere it already appears -- no call site changes.  The
// Makefile keeps the two modes in separate object trees so this should never
// fire; it is here because the failure it catches is silent and numerical.
#ifdef mpi_
inline namespace nerdss_build_mpi {
#else
inline namespace nerdss_build_serial {
#endif

typedef struct structMpiContext {  // holds MPI related data
  structMpiContext(int nRanks, int neigborRankBufferSize) {
    xOffset = 0;
    startCell = 0;
    endCell = 0;
    startGhosted = 0;
    endGhosted = 0;
    nprocs = nRanks;

    #ifdef mpi_
    // Reserve memory for sending/receiving requests regarding complexes:
    requestsSend = (MPI_Request*)malloc(sizeof(MPI_Request) * nRanks);
    requestsRecv = (MPI_Request*)malloc(sizeof(MPI_Request) * nRanks);
    #endif

    // Reserve memory for communication with neighboring ranks:
    sendBufferSize = neigborRankBufferSize;
    recvBufferSize = neigborRankBufferSize;
    MPIArrayFromLeft = (unsigned char*)malloc(recvBufferSize);
    MPIArrayToLeft = (unsigned char*)malloc(sendBufferSize);
    MPIArrayFromRight = (unsigned char*)malloc(recvBufferSize);
    MPIArrayToRight = (unsigned char*)malloc(sendBufferSize);
  }

  // The smallest x-column count a rank between two others can own.  A
  // step acts on each rank's columns in two halves, and in the first half
  // the right neighbour binds into this rank's last owned column (bin n-2)
  // through its own left ghost.  This rank's first half must therefore not
  // reach bin n-2: an initiator in the left region pairs with partners one
  // column further, so the left region has to end two columns short of it.
  // With the split of region_mid_index() a middle rank owning 2 columns
  // pairs bin 1 with bin 2 == bin n-2, and one owning 1 column acts on bin
  // n-2 itself; both ranks then bind the same molecule in the same half-step
  // and the exchange keeps one of the two bonds (see deserialize.cpp,
  // report_exchange_conflicts()).  Rank 0 can always avoid it by doing all
  // its work in the second half, and the last rank has no right neighbour,
  // so only the middle ranks are bound by this.
  static constexpr int minMiddleColumns = 3;

  // The most ranks that xTotal columns can be split between under the rule
  // above.
  static int max_ranks_for_columns(int xTotal) {
    if (xTotal < 2) return xTotal;
    // ends own 1 column each at the least, the middles minMiddleColumns
    int most = (xTotal - 2) / minMiddleColumns + 2;
    return most < 2 ? 2 : most;
  }

  // How many columns rank tempRank owns, and the first of them.  As even a
  // split as possible -- 13 columns onto 4 ranks is 4, 3, 3, 3 -- unless
  // that leaves a middle rank short of minMiddleColumns, in which case every
  // middle rank gets exactly that many and the two ends share the rest, the
  // first end taking the odd one.  The caller has checked that nprocs is at
  // most max_ranks_for_columns(xTotal).
  void column_share(int xTotal, int tempRank, int& first, int& count) const {
    int divCells = xTotal / nprocs;
    int modCells = xTotal % nprocs;
    if (nprocs < 3 || divCells >= minMiddleColumns) {
      first = divCells * tempRank + min2(tempRank, modCells);
      count = divCells + (tempRank < modCells ? 1 : 0);
      return;
    }
    const int ends = xTotal - minMiddleColumns * (nprocs - 2);
    const int firstEnd = ends - ends / 2;
    if (tempRank == 0) {
      first = 0;
      count = firstEnd;
    } else if (tempRank == nprocs - 1) {
      first = xTotal - ends / 2;
      count = ends / 2;
    } else {
      first = firstEnd + minMiddleColumns * (tempRank - 1);
      count = minMiddleColumns;
    }
  }

  void init_x_domain_and_offset(int xTotal, int& tempRank) {
    // Determine the domain of cells that tempRank is responsible for; see
    // column_share() for how the columns are divided.
    int count = 0;
    column_share(xTotal, tempRank, startCell, count);
    endCell = startCell + count - 1;
    // Include ghosted zones:
    startGhosted = startCell;
    endGhosted = endCell;
    if (startGhosted > 0) startGhosted--;
    if (endGhosted < xTotal - 1) endGhosted++;
    // xOffset represents x bin coordinate of most left cell belonging to
    // current rank; It is used for calculating currBin in update_memberMolLists
    // that requires substracting this offset.
    xOffset = startGhosted;
    // printf("== Rank %2d cell-domain:   startCell: %3d   endCell: %3d, diff+1:
    // %3d ==\n", tempRank, startCell, endCell, endCell - startCell + 1);
    //  determine the boundary of the rank in x axis
    int leftNonGhosted = 0;
    if (tempRank > 0) leftNonGhosted = xOffset + 1;
    xLeft = 1.0 * leftNonGhosted / xTotal;
    int rightNonGhosted = xOffset + (endCell - startCell + 1);
    if (tempRank > 0) rightNonGhosted = xOffset + 1 + (endCell - startCell + 1);
    xRight = 1.0 * rightNonGhosted / xTotal;
    if (tempRank == nprocs - 1) xRight = 1.0;
  }

  inline void check_buffer_size(int serializedSize) {
    if (serializedSize > sendBufferSize * 0.8) {
      std::cout << " --- WARNING ---" << std::endl;
      std::cout << "Serialized buffer size = " << serializedSize << std::endl;
      std::cout << "Sending buffer size = " << sendBufferSize << std::endl;
    }
    if (serializedSize > sendBufferSize) {
      std::cout << "Serialized buffer size = " << serializedSize << std::endl;
      std::cout << "Sending buffer size = " << sendBufferSize << std::endl;
      std::cout << "Please increase NEIGHBOR_BUFFER_SIZE" << std::endl;
      std::cout << "Exiting..." << std::endl;
      exit(1);
    }
  }
  void increase_recv_buffers() {
    if ((nMPIArrayFromLeft > 0.8 * recvBufferSize) ||
        (nMPIArrayFromRight > 0.8 * recvBufferSize)) {
      free(MPIArrayFromLeft);
      free(MPIArrayFromRight);
      recvBufferSize *= 1.2;
      MPIArrayFromLeft = (unsigned char*)malloc(recvBufferSize);
      MPIArrayFromRight = (unsigned char*)malloc(recvBufferSize);
    }
  }
  void increase_send_buffers() {
    if ((nMPIArrayToLeft > 0.8 * sendBufferSize) ||
        (nMPIArrayToRight > 0.8 * sendBufferSize)) {
      free(MPIArrayToLeft);
      free(MPIArrayToRight);
      sendBufferSize *= 1.2;
      MPIArrayToLeft = (unsigned char*)malloc(sendBufferSize);
      MPIArrayToRight = (unsigned char*)malloc(sendBufferSize);
    }
  }
  ~structMpiContext() {
    #ifdef mpi_
    free(requestsSend);
    free(requestsRecv);
    #endif
    free(MPIArrayFromLeft);
    free(MPIArrayToLeft);
    free(MPIArrayFromRight);
    free(MPIArrayToRight);
  }
  void print_spaces() {
    if (rank) {
      printf("%*c", rank, ' ');
      printf("%*c", rank, ' ');
    }
  }

  void print_xBins() {
    printf(
        "-----xBins of rank %d: startCell = %d, endCell = %d, startGhosted = "
        "%d, endGhosted = %d, xOffset = %d\n",
        rank, startCell, endCell, startGhosted, endGhosted, xOffset);
  }

  /*
  Function serialize serializes structMpiSynchronizationData struct
  into arrayRank of bytes.
  */
  void serialize(unsigned char* arrayRank, int& nArrayRank) {
    // Determines the domain of cells that tempRank is responsible for in gready
    // manner, e.g. divide 12 cells onto 5 ranks as: 3, 3, 2, 2, 2:
    PUSH(xOffset);
    PUSH(startCell);
    PUSH(endCell);
    PUSH(startGhosted);
    PUSH(endGhosted);
    PUSH(xLeft);
    PUSH(xRight);
  }
  /*
  Function deserialize deserializes structMpiSynchronizationData struct
  from arrayRank of bytes.
  */
  void deserialize(unsigned char* arrayRank, int& nArrayRank) {
    POP(xOffset);
    POP(startCell);
    POP(endCell);
    POP(startGhosted);
    POP(endGhosted);
    POP(xLeft);
    POP(xRight);
  }

  #ifdef mpi_
  // For non-blocking sending of actions:
  MPI_Request requestSendToLeft, requestSendToRight, requestRecvFromLeft,
      requestRecvFromRight;
  MPI_Status statusRecvFromLeft,
      statusRecvFromRight;  // for counting number of doubles received
  #endif

  int sendBufferSize, recvBufferSize;
  unsigned char* MPIArrayToRight;
  unsigned char* MPIArrayFromRight;
  unsigned char* MPIArrayToLeft;
  unsigned char* MPIArrayFromLeft;
  int nMPIArrayToRight;
  int nMPIArrayFromRight;
  int nMPIArrayToLeft;
  int nMPIArrayFromLeft;

  int xOffset;
  // Determine the domain of cells that tempRank is responsible for in gready
  // manner, e.g. divide 12 cells onto 5 ranks as: 3, 3, 2, 2, 2:
  int startCell, endCell;
  int startGhosted, endGhosted;
  double xLeft, xRight;

  int nprocs;

  int rank,
      size;  // MPI current rank and number of processing elements per width

  int simItr;
  bool checkUnimoleculeReactionPopulation;
  bool hasRankCommunicationForLargeComplex;

  #ifdef mpi_
  // MPI_Requests for sending requests to other ranks regarding complexes:
  MPI_Request *requestsSend, *requestsRecv;
  #endif

  // Vector of vectors will hold information one rank has to send to other
  // ranks:
  std::vector<std::vector<unsigned char> > toRank;

  Membrane* membraneObject;
  SimulVolume* simulVolume;
  std::vector<Molecule>* moleculeList;
  std::vector<Complex>* complexList;
} MpiContext;

#ifdef mpi_
}  // inline namespace nerdss_build_mpi
#else
}  // inline namespace nerdss_build_serial
#endif
