/*
      Fully Concurrent (FC) Communication Model prototype Code
      FC_Communications.cpp

Basic simulation parralelization steps:
1) scatter xMin-xMax: moleculs, cells, etc. + border-zones;
2...) iterate including border-zones;
3...) discard invalidated border zones and receive updated
4) collect by rank 0
5) restart files: collect->fprintf || fread->scatter

Scattering over 3 RANKS the XMIN-XMAX domain:
    0  1  2 Width
  |---------|
  | 0  1  2 |
  -----------

Scattering over both x and y axes might lead to faster execution for bigger
problems: Width 0  1  2 Height |---------| 0 | 0  5  0 | 1 | 6  0  8 | 2 | 1  7
0 |
       -----------

TODO: send only affected molecules

*/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <vector>

#define STRUCT_MOL_SIZE 17
// int NMOLS = 8;
int NMOLS = 1;

typedef struct struct_mpi_mol {
  /*
  Function serialize serializes changable portion of struct_mpi_mol
  into array of bytes.
  */
  void serialize(unsigned char *array) {
    array[0] = name;  // 0th byte contains only one character (name)
    *((double *)&(array[1])) = x;   // starting from the next (1st byte),
                                    // x coordinate is stored in next 8 bytes.
    *((double *)&(array[9])) = vx;  // starting from the next (9st byte),
                                    // speed of molecule in x direction
                                    // coordinate is stored in next 8 bytes.

    // update STRUCT_MOL_SIZE if needed
  }
  void deserialize(unsigned char *array) {
    name = array[0];  // 0th byte contains only one character (name)
    x = *((double *)&(array[1]));   // starting from the next (1st byte),
                                    // x coordinate is stored in next 8 bytes.
    vx = *((double *)&(array[9]));  // starting from the next (9st byte),
                                    // speed of molecule in x direction
                                    // coordinate is stored in next 8 bytes.
    processed = false;
    // update STRUCT_MOL_SIZE if needed
  }
  void print() {
    printf("%c:%4f, v=%4f, processed=%d\n", name, x, vx, processed);
  }

  int indexOrigin;
  char name;
  double x;
  double vx;
  bool forbidBind;
  bool processed;  // if processing a molecule in one vector resulted in
                   // shifting to another vector,
                   // when processing another vector,
                   // the molecule should not be processed again.
} MPI_MOL;

// Exchanging pieces and processing nation:
// All sending and receiving is done asynchronously.
// z part is always processed before x part of the following rank.

// Legend:
// A zone corresponds to a region in the x direction of width Rmax
//
// Assume rank # has all molecules with x positions between Xb and Xe (Xb<x<Xe)
// (We won't worry about using < and =< for this discussion.)

// molecules with x between Xb       and  Xb+Rmax  are within the "x"  zone.
// molecules with x between Xb+Rmax  and  Xb+2Rmax are within the "x*" zone.
// molecules with x between Xe-2Rmax and  Xe-Rmax  are within the "z*" zone.
// molecules with x between Xe-Rmax  and  Xe       are within the "z"  zone.

// molecules with x between  Xb+2Rmax and  (Xe+Xb)/2 are within the YL  zone.
// molecules with x between (Xe+Xb)/2 and   Xe-2Rmax are within the YR  zone.

//        zone x  zone x*  Y(left) Y(right)  zone z*  zone z
//      |---x---|---x*---|---YL---|---YR---|---z*---|---z---|

// Each rank has zones for their MolList which they perform simulations on.
// E.G.
//         rank 0              rank 1              rank 2

//    |x|x*|YL|YR|z*|z|   |x|x*|YL|YR|z*|z|   |x|x*|YL|YR|z*|z|
//     xb=0.0, xe=1.5      xb=1.5, xe=3.0      xb=3.0, xe=4.5

// With this distance zoning, one can now assume that a molecule
// in zone x will only interact with a molecule in x*, or the
// molecule z of the adjacent rank (rank to the Left, if exists). Likewise,
// molecules in zone z will only interact with molecules in zone z*
// and zone x of the adjacent rank (rank to the Right, if exists).
// Hence, information from an adjacent (Left) rank's z molecules must
// be made available (and exchanged at times) to determine the
// reactivity of x. Likewise, molecules in zone z must have
// information about molecules in the Right rank's x zone.
// Special storage, called ghost zones, are made available
// to hold that information, and is presented as a zone labeled
// with the zone name and a "g" suffix. Ghost storage is positioned
// next to the zone it will be used in evaluating reactivity,
// as shown below:

//           zone x  zone x*  Y(left) Y(right)  zone z*  zone z
//         |---x---|---x*---|---YL---|---YR---|---z*---|---z---|
//                               ||
//                               vv
//  zone zg  zone x  zone x*  Y(left) Y(right)  zone z*  zone z zone xg
// |---zg--|---x---|---x*---|---YL---|---YR---|---z*---|---z---|--xg---|
//   ghost                                                       ghost

// short form:
//               |zg|x|x*| YL | YR |z*|z|xg|

// Most of the computational work occurs in YL and YR, and is made
// independent of the computations on x and z for each iteration through
// the use of the x* an z* zone setup. This allows overlap of computation
// in YL and YR and communication between x and z zones and ghost zones.

// When discussing exchange, we should use MPI language:
// sending to a rank (send), receiving (recv) from a rank,
// or waiting for a rank (wait) with an additional value indicating
// the rank that the execution occurs "on".
// The usual communication pattern is a send/recv to/from ghost cells:

//   send(z  to 2 on 1 )   recv(zg from 1 on 2 )
//   send(xg to 2 on 1 )   recv(x  from 1 on 2 ),

// or if the communication is aggregated (packed in a single call), we use

//   send(z,xg to 2 on 1 )   recv(zg,x from 1 on 2).

// When discussing exchange, we should use MPI language:
// sending to a rank (send), receiving (recv) from a rank,
// or waiting for a rank (wait) with an additional value indicating
// the rank that the execution occurs "on".
// The usual communication pattern is a send/recv to/from ghost cells:

// Similarly, waiting for sending z and xg to the right to finish
// will be reffered to as:

//     wait(z,xg to 2 on 1)

// and waiting for receiving zg and x from the left to finish
// will be reffered to as:

//     wait(zg,x from 0 on 1)

// For sequential locations assigned to sequential ranks,
// as occurs here, one can use a send Left and send Right notation,
// when discussing "generic" (without ranks) communication:

//   send(z,xg to R )   recv(zg,x from L)
//   wait(z,xg to R)   wait(zg,x from L)

// For NERDSS, the communication will always occur in these pairs.
// That is, ghosts always receive from non-ghosts and vice versa. Hence,
// for recv(zg,x from 1 on 2) there will be a send(z,xg to 2 from 1).

typedef struct structMpiSynchronizationData {  // holds MPI related data
  // TODO: reallocate in runtime to accomodate needed memory
  // Allocate memory for arrays of bytes, each storing number of molecules and
  // serialized molecules for MPI (used for exchanging data using MPI)
  unsigned char *MPIArrayToRight;
  int nMPIArrayToRight;
  unsigned char *MPIArrayFromRight;
  int nMPIArrayFromRight;
  unsigned char *MPIArrayToLeft;
  int nMPIArrayToLeft;
  unsigned char *MPIArrayFromLeft;
  int nMPIArrayFromLeft;
} MpiSynchronizationData;

typedef struct structMpiContext {  // holds MPI related data
  void init() {
    MPI_Init(NULL,
             NULL);  // initializes MPI; all MPI programs must include this
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // receives the rank
    MPI_Comm_size(MPI_COMM_WORLD,
                  &size);  // and total number of processing elements

    MINX = 0;  // TODO: define x domain based on problem size
    MAXX = 4;
    RMAX = 0.1;

    DX = (MAXX - MINX) / static_cast<double>(size);

    // Calculating the domain in terms of 2D coordinates that current rank is
    // responsible for:
    minX = MINX + rank * DX;
    // if(minX > MINX)
    //     minX -= RMAX;
    maxX = MINX + (rank + 1) * DX;
    // if(maxX < MAXX)
    //     maxX += RMAX;
    // printf("cpu %d, iY:rank = %d:%d, minY:maxY = %f:%f, minX:maxX = %f:%f\n",
    //     mpiContext.rank, iY, rank, minY, maxY, minX, maxX);
  }
  void print_spaces() {
    if (rank) {
      printf("%*c", rank, ' ');
      printf("%*c", rank, ' ');
    }
  }

  inline int left() {  // finds rank of processing element responsible for
                       // portion of data to the left
    return (rank > 0) ? rank - 1 : -1;
  }

  inline int right() {  // finds rank of processing element responsible for
                        // portion of data to the right
    return (rank < size - 1) ? rank + 1 : -1;
  }

  void deserialize_molecules_from_0(
      unsigned char *arrayChar) {  // unpacks vector
    int nArrayRank = 0;
    int size = *((int *)&(arrayChar[nArrayRank]));
    nArrayRank += sizeof(int);
    for (int i = 0; i < size; i++) {
      MPI_MOL tempMol;
      tempMol.deserialize(&arrayChar[nArrayRank]);

      // |zg|x|x*| YL | YR |z*|z|xg|
      if (tempMol.x < (minX + maxX) / 2) {  // Left part
        if (tempMol.x < minX + RMAX) {      // zg or x
          if (tempMol.x < minX)             // zg
            vZg.push_back(tempMol);
          else
            vX.push_back(tempMol);
        } else {                            // x* or YL
          if (tempMol.x < minX + 2 * RMAX)  // x*
            vXa.push_back(tempMol);
          else
            vYL.push_back(tempMol);
        }
      } else {                              // Right part
        if (tempMol.x < maxX - RMAX) {      // YR or z*
          if (tempMol.x < maxX - 2 * RMAX)  // Yr
            vYR.push_back(tempMol);
          else
            vZa.push_back(tempMol);
        } else {                 // z or xg
          if (tempMol.x < maxX)  // z
            vZ.push_back(tempMol);
          else
            vXg.push_back(tempMol);
        }
      }
      nArrayRank += STRUCT_MOL_SIZE;
    }
  }

  void reset_processed() {  // No mol should be processed twice per iteration
    for (auto &m : vZg) m.processed = false;
    for (auto &m : vX) m.processed = false;
    for (auto &m : vXa) m.processed = false;
    for (auto &m : vYL) m.processed = false;
    for (auto &m : vYR) m.processed = false;
    for (auto &m : vZa) m.processed = false;
    for (auto &m : vZ) m.processed = false;
    for (auto &m : vXg) m.processed = false;
  }

  MpiSynchronizationData m;

  double MINX, MAXX;  // 2D domain (or 3D) that should be processed in parallel
  double DX;
  double RMAX;

  int rank,
      size;  // MPI current rank and number of processing elements per width

  // Storing the domain in terms of coordinates that current rank is responsible
  // for
  double minX, maxX;

  std::vector<MPI_MOL>
      vZg;  // includes molecules that belong to the z ghost zone
  std::vector<MPI_MOL> vX;   // includes molecules that belong to the x zone
  std::vector<MPI_MOL> vXa;  // includes molecules that belong to the x* zone
  std::vector<MPI_MOL> vYL;  // includes molecules that belong to the YL zone
  std::vector<MPI_MOL> vYR;  // includes molecules that belong to the YR zone
  std::vector<MPI_MOL> vZa;  // includes molecules that belong to the z* zone
  std::vector<MPI_MOL> vZ;   // includes molecules that belong to the z zone
  std::vector<MPI_MOL>
      vXg;  // includes molecules that belong to the x ghost zone

  std::vector<int>
      ghostedMolOwnerIndex;  // index of ghosted molecule on owner rank
  std::vector<int> ghostedMolOwnerRank;  // owner rank of ghosted molecule
} mpiContext;

// initializes moleculs
void init_mols(std::vector<MPI_MOL> &vMols, mpiContext &mpiContext) {
  for (int i = 0; i < NMOLS; i++) {
    MPI_MOL tempMol;
    tempMol.name = 'A' + i;
    tempMol.x = mpiContext.MINX + 3.51 + i * 0.5;
    tempMol.vx = 0;
    tempMol.forbidBind = false;
    vMols.push_back(tempMol);
  }
  vMols[0].vx = -0.04;
}

// prints vector of molecules
void print_mols(std::vector<MPI_MOL> array) {
  for (MPI_MOL &m : array) {
    m.print();
  }
}

// prints the array
void print_array(double array[], int n) {
  for (int i = 0; i < n; i++) printf("%6f ", array[i]);
  printf("\n");
}

int serialize_molecules(std::vector<MPI_MOL> vectorMol,
                        unsigned char *arrayChar) {
  int nArrayRank = 0;
  *((int *)&(arrayChar[nArrayRank])) = vectorMol.size();
  nArrayRank += sizeof(int);
  for (int i = 0; i < vectorMol.size(); i++) {
    vectorMol[i].serialize(&arrayChar[nArrayRank]);
    nArrayRank += STRUCT_MOL_SIZE;
  }
  return nArrayRank;
}

int serialize_molecules(std::vector<MPI_MOL> vMol, unsigned char *arrayChar,
                        double tempMinX, double tempMaxX) {
  int nArrayRank = 0;
  nArrayRank += sizeof(int);
  for (int i = 0; i < vMol.size(); i++) {
    if ((vMol[i].x >= tempMinX) && (vMol[i].x < tempMaxX)) {
      vMol[i].serialize(&arrayChar[nArrayRank]);
      nArrayRank += STRUCT_MOL_SIZE;
    }
  }
  *((int *)&(arrayChar[0])) = (nArrayRank - sizeof(int)) / STRUCT_MOL_SIZE;
  return nArrayRank;
}

int deserialize_molecules(unsigned char *arrayChar,
                          std::vector<MPI_MOL> &vectorMol,
                          bool clearVectorFirst = true) {
  int nArrayRank = 0;
  int size = *((int *)&(arrayChar[nArrayRank]));
  nArrayRank += sizeof(int);
  if (clearVectorFirst) vectorMol.clear();
  for (int i = 0; i < size; i++) {
    MPI_MOL tempMol;
    tempMol.deserialize(&arrayChar[nArrayRank]);
    vectorMol.push_back(tempMol);
    nArrayRank += STRUCT_MOL_SIZE;
  }
  return nArrayRank;
}

void test() {  // TODO: optionally upgrade
  // test molecule serialization and deserialization
  MPI_MOL m;
  m.name = 'A';
  m.x = 123.4;
  m.vx = -0.123;
  m.forbidBind = false;
  m.processed = false;
  unsigned char arr[18];
  m.serialize(arr);
  // for(int tt=0;tt<18;tt++)
  //     printf("%d ", arr[tt]);
  // printf("\n");
  MPI_MOL tempMol;
  tempMol.deserialize(arr);
  // if(!tempMol.forbidBind)
  //     printf("%c:%4f, %4f\n", tempMol.name, tempMol.x, tempMol.vx);

  // test MPI_MOL vector serialization and deserialization
  std::vector<MPI_MOL> vMols;
  vMols.push_back(m);
  vMols.push_back(m);
  vMols.push_back(m);
  // Allocate memory for array of bytes storing number of molecules and
  // serialized molecules for MPI
  unsigned char *arrayChar =
      (unsigned char *)malloc(sizeof(int) + vMols.size() * STRUCT_MOL_SIZE);
  print_mols(vMols);
  int nBytes = serialize_molecules(vMols, arrayChar);
  printf("nBytes = %d\n", nBytes);
  std::vector<MPI_MOL> vTemp;
  vTemp.clear();
  deserialize_molecules(arrayChar, vTemp);
  print_mols(vTemp);
  // Free memory reserved for array of bytes storing number of molecules and
  // serialized molecules for MPI
  free(arrayChar);
}

void generate_test_vector_with_single_molelule(char name, double x,
                                               std::vector<MPI_MOL> &vTemp) {
  MPI_MOL mTemp;
  mTemp.name = name;
  mTemp.x = x;
  mTemp.vx = 5.4321;
  mTemp.forbidBind = false;

  vTemp.clear();
  vTemp.push_back(mTemp);
}

/*
Processing elements of vector v,
having in mind possible reactions
with elements on the left and right from v.
vL - vector of elements on the left
v -  vector of elements in the middle
vR - vector of elements on the right
x1 and x2 - border x coordinates between vector elements
*/
void process_vector(std::string label, std::vector<MPI_MOL> &vL,
                    std::vector<MPI_MOL> &v, std::vector<MPI_MOL> &vR,
                    double x1, double x2) {
  // std::cout << "Processing vector " << label << std::endl;
  for (int i = 0; i < v.size(); i++) {
    if (!v[i].processed) {
      // std::cout << label << ": ";
      //
      v[i].print();
      v[i].x += v[i].vx;
      v[i].processed = true;
      if (v[i].x < x1) {
        //
        printf("Move left %lf\n", v[i].x);
        vL.push_back(v[i]);
        v.erase(v.begin() +
                i);  // TODO: if order of molecules is not important,
                     // remove last after copying it into i
      } else {
        if (v[i].x > x2) {
          //
          printf("Move right %lf\n", v[i].x);
          vR.push_back(v[i]);
          v.erase(v.begin() +
                  i);  // TODO: if order of molecules is not important,
                       // remove last after copying it into i
        }
      }
    }
  }
}

void process_iteration(mpiContext &mpiContext,
                       void (*process_vector_function)(std::string label,
                                                       std::vector<MPI_MOL> &vL,
                                                       std::vector<MPI_MOL> &v,
                                                       std::vector<MPI_MOL> &vR,
                                                       double x1, double x2)) {
  MpiSynchronizationData &m = mpiContext.m;
  std::string label;

  mpiContext
      .reset_processed();  // no molecule should be procesed twice per iteration
  // Order of execution for one MPI rank (explained):
  // 1) At the beginning of each iteration,
  //    each rank owns updated x border-zone
  //    (from the end of last iteration),
  //    but not z border-zone, so he has to
  //    send zone x with zg to the previous/left rank
  //    and receive zone xg with z from the next/right rank:
  //    send(zg,x to L), recv(z,xg from R)
  // 2) As the x border-zone is already sent to the left rank
  //    (so that it can process this border-zone, and send it back),
  //    receiving updated x border-zone is necessary before processing zone x:
  //    recv(zg,x from L)
  // 3) While waiting for communication from one node to another to finish,
  //    do aproximately half of independent calculations:
  //    process(yL)
  // 4) After receiving right border-zone, update it and send it back;
  //    Note: z* has to be processed before sending z;
  //          Otherwise z* would affect z that is already sent to N.
  //    Wait for z, process z* and z, send z,xg to the following rank (right):
  //    wait(z,xg from R), process(z*,z), send(z,xg to R);
  // 5) While waiting for second communication from one node to another to
  // finish,
  //    do the rest of independent calculations (aproximately second half):
  //    process(yR)
  // 6) Wait for receiving 2) to finish, process x1,x1+, and wait for sending
  // z1,x2 to finish
  //    wait(zg,x from L), process(x*,x), wait(z,xg to R); wait(zg,x to L)

  // Order of execution for one MPI rank (depicted):
  // |zg|x|x*| YL | YR |z*|z|xg|
  //  S___
  //                       R___
  //  R___
  //          P___
  //                       W___
  //                    P___
  //                       S___
  //               P___
  //  W___
  //     P___
  //                       W___

  // Sending & receiving are all non-blocking.
  MPI_Request request1SendToLeft, request2RecvFromRight, request3RecvFromLeft,
      request4SendToRight;  // for non-blocking receiving
  MPI_Status statusRecvFromRight,
      statusRecvFromLeft;  // for counting number of doubles received
  // 1) At the beginning of each iteration,
  //    each rank owns updated x border-zone
  //    (from the end of last iteration),
  //    but not z border-zone, so he has to
  //    send zone x with zg to the previous/left rank
  //    and receive zone xg with z from the next/right rank:
  //    send(zg,x to L), recv(z,xg from R)

  // send(zg,x to L)
  if (mpiContext.rank) {  // if not first rank
    int nTemp = serialize_molecules(mpiContext.vZg, &m.MPIArrayToLeft[0]);
    m.nMPIArrayToLeft =
        serialize_molecules(mpiContext.vX, &m.MPIArrayToLeft[nTemp]);
    m.nMPIArrayToLeft += nTemp;
    // Send portion of matrix to processing element of rank number tempRank:

    // TODO: discuss

    MPI_Isend(m.MPIArrayToLeft, m.nMPIArrayToLeft, MPI_CHAR,
              mpiContext.rank - 1, 0, MPI_COMM_WORLD, &request1SendToLeft);
  }

  // recv(z,xg from R)
  if (mpiContext.rank < mpiContext.size - 1) {  // if not last
    MPI_Irecv(m.MPIArrayFromRight, 100, MPI_CHAR, mpiContext.rank + 1, 0,
              MPI_COMM_WORLD, &request2RecvFromRight);
  }

  // 2) As the x border-zone is already sent to the left rank
  //    (so that it can process this border-zone, and send it back),
  //    receiving updated x border-zone is necessary before processing zone x:
  //    recv(zg,x from L)
  if (mpiContext.rank) {  // if not first rank
    MPI_Irecv(m.MPIArrayFromLeft, 100, MPI_CHAR, mpiContext.rank - 1, 0,
              MPI_COMM_WORLD, &request3RecvFromLeft);
  }

  // 3) While waiting for communication from one node to another to finish,
  //    do aproximately half of independent calculations:
  //    process(yL)

  label = "yL";
  process_vector(label, mpiContext.vXa, mpiContext.vYL, mpiContext.vYR,
                 mpiContext.minX + 2 * mpiContext.RMAX,
                 (mpiContext.minX + mpiContext.maxX) / 2);

  // 4) After receiving right border-zone, update it and send it back;
  //    Note: z* has to be processed before sending z;
  //          Otherwise z* would affect z that is already sent to N.
  //    Wait for z, process z* and z, send z,xg to the following rank (right):
  //    wait(z,xg from R), process(z*,z), send(z,xg to R);
  // send(z,xg to R);
  if (mpiContext.rank < mpiContext.size - 1) {  // if not last rank
    MPI_Wait(&request2RecvFromRight,
             &statusRecvFromRight);  // waiting for receiving to finish
    // determine number of received bytes:
    MPI_Get_count(&statusRecvFromRight, MPI_CHAR, &m.nMPIArrayFromRight);
    // Two arrays are received. Deserialize first vector,
    // and get number of parsed bytes from received array:
    int tempRead = deserialize_molecules(m.MPIArrayFromRight, mpiContext.vZ);
    // Starting after last parsed byte, deserialize second vector:
    deserialize_molecules(m.MPIArrayFromRight + tempRead, mpiContext.vXg);
  }

  // process(z*,z)
  label = "z* zone";
  process_vector(label, mpiContext.vYR, mpiContext.vZa, mpiContext.vZ,
                 mpiContext.maxX - 2 * mpiContext.RMAX,
                 mpiContext.maxX - mpiContext.RMAX);
  label = "z zone";
  process_vector(label, mpiContext.vZa, mpiContext.vZ, mpiContext.vXg,
                 mpiContext.maxX - mpiContext.RMAX, mpiContext.maxX);

  // send(z,xg to R);
  if (mpiContext.rank < mpiContext.size - 1) {  // if not last rank
    int nTemp = serialize_molecules(mpiContext.vZ, &m.MPIArrayToRight[0]);
    m.nMPIArrayToRight =
        serialize_molecules(mpiContext.vXg, &m.MPIArrayToRight[nTemp]);
    m.nMPIArrayToRight += nTemp;
    MPI_Isend(m.MPIArrayToRight, m.nMPIArrayToRight, MPI_CHAR,
              mpiContext.rank + 1, 0, MPI_COMM_WORLD, &request4SendToRight);
  }

  // 5) While waiting for second communication from one node to another to
  // finish,
  //    do the rest of independent calculations (aproximately second half):
  //    process(yR)

  label = "yR";
  process_vector(label, mpiContext.vYL, mpiContext.vYR, mpiContext.vZa,
                 (mpiContext.minX + mpiContext.maxX) / 2,
                 mpiContext.maxX - 2 * mpiContext.RMAX);

  // 6) Wait for sending 2) to finish, process x1,x1+, and wait for sending
  // z1,x2 to finish
  //    wait(zg,x from L), process(x*,x), wait(z,xg to R), wait(zg,x to L)

  // wait(zg,x from L)
  if (mpiContext.rank) {  // if not first rank
    MPI_Wait(&request3RecvFromLeft,
             &statusRecvFromLeft);  // waiting for receiving to finish
    MPI_Get_count(&statusRecvFromLeft, MPI_CHAR, &m.nMPIArrayFromLeft);
    int tempRead = deserialize_molecules(m.MPIArrayFromLeft, mpiContext.vZg);
    deserialize_molecules(m.MPIArrayFromLeft + tempRead, mpiContext.vX);
  }

  // process(x*,x)
  label = "x zone";
  process_vector(label, mpiContext.vZg, mpiContext.vX, mpiContext.vXa,
                 mpiContext.minX, mpiContext.minX + mpiContext.RMAX);
  label = "x* zone";
  process_vector(label, mpiContext.vX, mpiContext.vXa, mpiContext.vYL,
                 mpiContext.minX + mpiContext.RMAX,
                 mpiContext.minX + 2 * mpiContext.RMAX);

  // wait(z,xg to R);
  if (mpiContext.rank < mpiContext.size - 1) {  // if not last rank
    MPI_Wait(&request4SendToRight,
             MPI_STATUS_IGNORE);  // waiting for sending to finish
  }

  // wait(zg,x to L)
  if (mpiContext.rank) {  // if not last rank
    MPI_Wait(&request1SendToLeft,
             MPI_STATUS_IGNORE);  // waiting for sending to finish
  }

  //    std::cout << "Here" << std::endl;
}

int main(int argc, char **argv) {
  mpiContext mpiContext;  // storing all context related to MPI
  MpiSynchronizationData &m = mpiContext.m;
  std::string label;            // holding temporary strings
  std::vector<MPI_MOL> vOwned;  // includes molecules that belong to the rank
  std::vector<MPI_MOL>
      vExtended;  // includes molecules that belong to the rank and border zones

  // TODO: reallocate in runtime to accomodate needed memory
  // Allocate memory for arrays of bytes, each storing number of molecules and
  // serialized molecules for MPI (used for exchanging data using MPI)
  m.MPIArrayToRight = (unsigned char *)malloc(1000);
  m.nMPIArrayToRight;
  m.MPIArrayFromRight = (unsigned char *)malloc(1000);
  m.nMPIArrayFromRight;
  m.MPIArrayToLeft = (unsigned char *)malloc(1000);
  m.nMPIArrayToLeft;
  m.MPIArrayFromLeft = (unsigned char *)malloc(1000);
  m.nMPIArrayFromLeft;

  mpiContext.init();

  // if(!mpiContext.rank)
  //     test();

  // Allocate memory for array of bytes storing number of molecules and
  // serialized molecules (used for exchanging data with rank 0 using MPI)
  unsigned char *MPIArrayFrom0 = (unsigned char *)malloc(1000);
  int nMPIArrayFrom0;

  // Rank 0 scatters:
  if (!mpiContext.rank) {
    init_mols(vOwned, mpiContext);
    // print_mols(vOwned);
    //  Result:
    // A:0.250000, v=0.100000
    // B:0.750000, v=0.000000
    // C:1.250000, v=0.000000
    // D:1.750000, v=0.000000
    // E:2.250000, v=0.000000
    // F:2.750000, v=0.000000
    // G:3.250000, v=0.000000
    // H:3.750000, v=0.000000

    for (int tempRank = mpiContext.size - 1; tempRank >= 0; tempRank--) {
      double tempMinX = mpiContext.MINX + tempRank * mpiContext.DX;
      if (tempRank > 0) tempMinX -= mpiContext.RMAX;
      double tempMaxX = mpiContext.MINX + (tempRank + 1) * mpiContext.DX;
      if (tempRank < mpiContext.size - 1) tempMaxX += mpiContext.RMAX;
      // printf("tempRank = %d, domain X = %f : %f\n", tempRank, tempMinX,
      // tempMaxX);
      //  Result:
      // tempRank = 0, domain X = 0.000000 : 1.100000
      // tempRank = 1, domain X = 0.900000 : 2.100000
      // tempRank = 2, domain X = 1.900000 : 3.100000
      // tempRank = 3, domain X = 2.900000 : 4.000000

      // Store portion of molecules for rank number tempRank into the tempArray:
      nMPIArrayFrom0 =
          serialize_molecules(vOwned, &MPIArrayFrom0[0], tempMinX, tempMaxX);
      // Send portion of matrix to processing element of rank number tempRank:
      if (tempRank)
        MPI_Send(MPIArrayFrom0, nMPIArrayFrom0, MPI_CHAR, tempRank, 0,
                 MPI_COMM_WORLD);
    }
  }
  vOwned.clear();

  // Receiving pieces on each rank from rank 0:
  if (mpiContext.rank) {
    MPI_Status status;
    MPI_Recv(MPIArrayFrom0, 1000, MPI_CHAR, 0, 0, MPI_COMM_WORLD,
             &status);  // Each rank receives piece from rank 0
    MPI_Get_count(&status, MPI_CHAR, &nMPIArrayFrom0);
  }
  deserialize_molecules(MPIArrayFrom0, vExtended);
  mpiContext.deserialize_molecules_from_0(MPIArrayFrom0);
  // Free memory reserved for array of bytes storing number of molecules and
  // serialized molecules for MPI
  free(MPIArrayFrom0);

  // if(mpiContext.rank == 2) print_mols(vExtended);
  // Result:
  // C:1.250000, v=0.000000
  // G:3.250000, v=0.000000
  // E:2.250000, v=0.000000
  // A:0.250000, v=0.100000
  // D:1.750000, v=0.000000
  // H:3.750000, v=0.000000
  // F:2.750000, v=0.000000
  // B:0.750000, v=0.000000

  void (*process_vector_function)(
      std::string label, std::vector<MPI_MOL> & vL, std::vector<MPI_MOL> & v,
      std::vector<MPI_MOL> & vR, double x1, double x2);
  process_vector_function = &process_vector;

  // Simulation iterations:
  for (int nIter = 0; nIter < 70; nIter++) {
    process_iteration(mpiContext, process_vector_function);

  }  // for nIter

  free(m.MPIArrayFromLeft);
  free(m.MPIArrayToLeft);
  free(m.MPIArrayFromRight);
  free(m.MPIArrayToRight);

  // Returning results to rank 0:
  // Allocate memory for array of bytes storing number of molecules and
  // serialized molecules (used for exchanging data with rank 0 using MPI)
  unsigned char *MPIArrayTo0 = (unsigned char *)malloc(1000);
  int nMPIArrayTo0;

  nMPIArrayTo0 = 0;
  nMPIArrayTo0 +=
      serialize_molecules(mpiContext.vX, &MPIArrayTo0[nMPIArrayTo0]);
  nMPIArrayTo0 +=
      serialize_molecules(mpiContext.vXa, &MPIArrayTo0[nMPIArrayTo0]);
  nMPIArrayTo0 +=
      serialize_molecules(mpiContext.vYL, &MPIArrayTo0[nMPIArrayTo0]);
  nMPIArrayTo0 +=
      serialize_molecules(mpiContext.vYR, &MPIArrayTo0[nMPIArrayTo0]);
  nMPIArrayTo0 +=
      serialize_molecules(mpiContext.vZa, &MPIArrayTo0[nMPIArrayTo0]);
  nMPIArrayTo0 +=
      serialize_molecules(mpiContext.vZ, &MPIArrayTo0[nMPIArrayTo0]);
  if (mpiContext.rank)  // don't send if rank 0
    MPI_Send(MPIArrayTo0, nMPIArrayTo0, MPI_CHAR, 0, 0,
             MPI_COMM_WORLD);  // sending all molecules to rank 0

  // Gathering results:
  if (mpiContext.rank == 0) {  // only rank 0 receives all pieces of matrix from
                               // processing elements

    for (int tempRank = 0; tempRank < mpiContext.size; tempRank++) {
      // Receive data from processing element tempRank
      // (possible non-blocking receiving and waiting for finishing instead of
      // blocking receiving):
      MPI_Status status;
      if (tempRank) {
        MPI_Recv(MPIArrayTo0, 10000, MPI_CHAR, tempRank, 0, MPI_COMM_WORLD,
                 &status);
        MPI_Get_count(&status, MPI_CHAR, &nMPIArrayTo0);
      }

      int tempRead = 0;
      tempRead += deserialize_molecules(
          MPIArrayTo0 + tempRead, vOwned,
          false);  // false => append to vOwned, don't clear it
      tempRead += deserialize_molecules(
          MPIArrayTo0 + tempRead, vOwned,
          false);  // false => append to vOwned, don't clear it
      tempRead += deserialize_molecules(
          MPIArrayTo0 + tempRead, vOwned,
          false);  // false => append to vOwned, don't clear it
      tempRead += deserialize_molecules(
          MPIArrayTo0 + tempRead, vOwned,
          false);  // false => append to vOwned, don't clear it
      tempRead += deserialize_molecules(
          MPIArrayTo0 + tempRead, vOwned,
          false);  // false => append to vOwned, don't clear it
      tempRead += deserialize_molecules(
          MPIArrayTo0 + tempRead, vOwned,
          false);  // false => append to vOwned, don't clear it
    }
    // print_mols(vOwned);
    //  Result:
  }

  free(MPIArrayTo0);

  MPI_Finalize();  // finalize MPI; all MPI programs must include this
  return 0;
}
