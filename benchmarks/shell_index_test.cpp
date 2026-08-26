/*! \file shell_index_test.cpp
 * \brief Brute-force completeness check for ShellIndex.
 *
 * The lat-long stencil is only correct if, for every pair of surface points
 * within the cutoff, one of them holds the other in its own cell or in one of
 * its forward neighbours -- and if no pair is held twice, since a pair offered
 * twice would be tested twice.  Both are checked here against an O(N^2) sweep,
 * over radius-and-cutoff combinations from mild curvature (a 1000 nm sphere
 * against a 30 nm cutoff) to extreme (a 30 nm sphere against 19.4).
 *
 * Build and run:
 *   g++ -O2 -std=c++0x -Iinclude $(gsl-config --cflags) \
 *       src/classes/class_ShellIndex.cpp src/classes/class_Vec3D.cpp \
 *       benchmarks/shell_index_test.cpp -o shell_index_test $(gsl-config --libs)
 *   ./shell_index_test
 *
 * Exits non-zero if any pair is missed or double counted.
 */
#include "classes/class_ShellIndex.hpp"
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>

static double frand() { return double(rand()) / RAND_MAX; }

static int check(double R, double h, int nPts, unsigned seed)
{
    ShellIndex idx;
    idx.build(R, h, true);
    if (!idx.active) { printf("  R=%-7g h=%-6g INACTIVE\n", R, h); return 0; }

    srand(seed);
    std::vector<Vec3D> pts;
    std::vector<int> cell;
    for (int i = 0; i < nPts; ++i) {
        // uniform on the sphere
        double z = 2 * frand() - 1, phi = 2 * M_PI * frand(), r = sqrt(1 - z * z);
        const double rf = idx.radiusFloor; // worst admitted radius
        Vec3D p; p.x = rf * r * cos(phi); p.y = rf * r * sin(phi); p.z = rf * z;
        pts.push_back(p);
        cell.push_back(idx.cell_of(p));
    }

    // what the stencil offers
    std::set<std::pair<int,int>> offered;
    long long offeredCount = 0;
    std::map<int, std::vector<int>> members;
    for (int i = 0; i < nPts; ++i) members[cell[i]].push_back(i);
    for (std::map<int, std::vector<int>>::iterator it = members.begin(); it != members.end(); ++it) {
        const std::vector<int>& m = it->second;
        for (size_t a = 0; a < m.size(); ++a) {
            for (size_t b = a + 1; b < m.size(); ++b) {
                offered.insert(std::make_pair(std::min(m[a],m[b]), std::max(m[a],m[b]))); ++offeredCount;
            }
            for (size_t n = 0; n < idx.neighborList[it->first].size(); ++n) {
                int nb = idx.neighborList[it->first][n];
                std::map<int, std::vector<int>>::iterator jt = members.find(nb);
                if (jt == members.end()) continue;
                for (size_t b = 0; b < jt->second.size(); ++b) {
                    int q = jt->second[b];
                    offered.insert(std::make_pair(std::min(m[a],q), std::max(m[a],q))); ++offeredCount;
                }
            }
        }
    }

    // Ground truth is the property the simulation actually relies on: a pair
    // whose straight-line separation is within the cutoff must be offered.
    // Points are placed at the radius floor, which is the worst case -- the
    // same chord subtends a wider angle the closer in it sits.
    long long inRange = 0, missed = 0, angleMissed = 0;
    for (int i = 0; i < nPts; ++i)
        for (int j = i + 1; j < nPts; ++j) {
            double dx = pts[i].x - pts[j].x, dy = pts[i].y - pts[j].y, dz = pts[i].z - pts[j].z;
            if (sqrt(dx*dx + dy*dy + dz*dz) > h) continue;
            ++inRange;
            if (!offered.count(std::make_pair(i, j))) ++missed;
        }
    // and separately that the stencil covers its own declared angle
    for (int i = 0; i < nPts; ++i)
        for (int j = i + 1; j < nPts; ++j) {
            double d = (pts[i].x*pts[j].x + pts[i].y*pts[j].y + pts[i].z*pts[j].z);
            double n = sqrt((pts[i].x*pts[i].x+pts[i].y*pts[i].y+pts[i].z*pts[i].z)
                          * (pts[j].x*pts[j].x+pts[j].y*pts[j].y+pts[j].z*pts[j].z));
            d /= n; if (d > 1) d = 1; if (d < -1) d = -1;
            if (acos(d) > idx.gammaCut) continue;
            if (!offered.count(std::make_pair(i, j))) ++angleMissed;
        }

    long long dup = offeredCount - (long long)offered.size();
    printf("  R=%-7g h=%-6g bands=%-3d cells=%-6d inCutoff=%-7lld missed=%-3lld angleMissed=%-3lld doubleCounted=%lld\n",
           R, h, idx.nBands, idx.totalCells, inRange, missed, angleMissed, dup);
    return (missed != 0 || angleMissed != 0 || dup != 0) ? 1 : 0;
}

int main()
{
    int bad = 0;
    printf("ShellIndex stencil check (points uniform on the shell)\n");
    bad += check(70.0,  19.4476, 4000, 1);   // gagsphere
    bad += check(100.0,  9.7704, 4000, 2);   // sphere
    bad += check(1000.0, 30.0,   6000, 3);   // large sphere
    bad += check(50.0,  19.4476, 4000, 4);   // strong curvature
    bad += check(30.0,  19.4476, 3000, 5);   // cutoff a large fraction of R
    bad += check(500.0,  5.0,    6000, 6);   // many small cells
    bad += check(70.0,   2.0,    6000, 7);   // very fine
    printf("%s\n", bad ? "FAILURES PRESENT" : "all cases: complete and single-counted");
    return bad;
}
