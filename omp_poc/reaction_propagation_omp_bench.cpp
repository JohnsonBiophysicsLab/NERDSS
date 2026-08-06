#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct ReactionPair {
  double x1, y1, z1, x2, y2, z2;
  double d1x, d1y, d1z, d2x, d2y, d2z;
  double dr1z, dr2z, mag1, mag2, bindRadius;
};

struct ReactionResult {
  double distance, separation, effectiveDiffusion;
  unsigned char withinRmax;
};

struct Motion {
  double cx, cy, cz, tx, ty, tz, rx, ry, rz;
};

struct Quat {
  double x, y, z, w;
};

struct Point {
  double x, y, z;
  std::uint32_t complex;
};

inline ReactionResult checkReaction(const ReactionPair& p, double dt) {
  double d = (p.d1x + p.d2x + p.d1y + p.d2y + p.d1z + p.d2z) / 3.0;
  bool flat = std::fabs(p.d1z) < 1e-10;
  d += 2 * p.mag1 * (1 - std::cos(std::sqrt((flat ? 2.0 : 4.0) * p.dr1z * dt)))
      / ((flat ? 4.0 : 6.0) * dt);
  flat = std::fabs(p.d2z) < 1e-10;
  d += 2 * p.mag2 * (1 - std::cos(std::sqrt((flat ? 2.0 : 4.0) * p.dr2z * dt)))
      / ((flat ? 4.0 : 6.0) * dt);
  const double rmax = 3 * std::sqrt(6 * d * dt) + p.bindRadius;
  const double dx = p.x1 - p.x2;
  const double dy = p.y1 - p.y2;
  const double dz = p.z1 - p.z2;
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  return {r, r - p.bindRadius, d, static_cast<unsigned char>(r < rmax)};
}

inline Quat makeQuat(const Motion& m) {
  const double cz = std::cos(m.rz * .5), sz = std::sin(m.rz * .5);
  const double cy = std::cos(m.ry * .5), sy = std::sin(m.ry * .5);
  const double cx = std::cos(m.rx * .5), sx = std::sin(m.rx * .5);
  Quat q{sx * cy * cz - cx * sy * sz,
         cx * sy * cz + sx * cy * sz,
         cx * cy * sz - sx * sy * cz,
         cx * cy * cz + sx * sy * sz};
  const double scale = 1 / std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
  q.x *= scale;
  q.y *= scale;
  q.z *= scale;
  q.w *= scale;
  return q;
}

inline Point propagate(const Point& p, const Motion& m, const Quat& q) {
  const double x = p.x - m.cx;
  const double y = p.y - m.cy;
  const double z = p.z - m.cz;
  const double tx = 2 * (q.y * z - q.z * y);
  const double ty = 2 * (q.z * x - q.x * z);
  const double tz = 2 * (q.x * y - q.y * x);
  return {x + q.w * tx + q.y * tz - q.z * ty + m.cx + m.tx,
          y + q.w * ty + q.z * tx - q.x * tz + m.cy + m.ty,
          z + q.w * tz + q.x * ty - q.y * tx + m.cz + m.tz,
          p.complex};
}

void reactionSerial(const std::vector<ReactionPair>& input,
                    std::vector<ReactionResult>& output, double dt) {
  for (std::size_t k = 0; k < input.size(); ++k) {
    output[k] = checkReaction(input[k], dt);
  }
}

void reactionOpenMp(const std::vector<ReactionPair>& input,
                    std::vector<ReactionResult>& output, double dt,
                    int threads) {
  const std::int64_t count = static_cast<std::int64_t>(input.size());
#pragma omp parallel for schedule(static) num_threads(threads)
  for (std::int64_t k = 0; k < count; ++k) {
    output[static_cast<std::size_t>(k)] =
        checkReaction(input[static_cast<std::size_t>(k)], dt);
  }
}

void propagationSerial(const std::vector<Motion>& motions,
                       const std::vector<Point>& input,
                       std::vector<Quat>& quaternions,
                       std::vector<Point>& output) {
  for (std::size_t k = 0; k < motions.size(); ++k) {
    quaternions[k] = makeQuat(motions[k]);
  }
  for (std::size_t k = 0; k < input.size(); ++k) {
    output[k] = propagate(input[k], motions[input[k].complex],
                          quaternions[input[k].complex]);
  }
}

void propagationOpenMp(const std::vector<Motion>& motions,
                       const std::vector<Point>& input,
                       std::vector<Quat>& quaternions,
                       std::vector<Point>& output, int threads) {
  const std::int64_t motionCount = static_cast<std::int64_t>(motions.size());
  const std::int64_t pointCount = static_cast<std::int64_t>(input.size());
#pragma omp parallel num_threads(threads)
  {
#pragma omp for schedule(static)
    for (std::int64_t k = 0; k < motionCount; ++k) {
      quaternions[static_cast<std::size_t>(k)] =
          makeQuat(motions[static_cast<std::size_t>(k)]);
    }
#pragma omp for schedule(static)
    for (std::int64_t k = 0; k < pointCount; ++k) {
      const auto index = static_cast<std::size_t>(k);
      output[index] = propagate(input[index], motions[input[index].complex],
                                quaternions[input[index].complex]);
    }
  }
}

template <class Function>
double medianMilliseconds(Function&& function, int repetitions, int trials) {
  std::vector<double> samples;
  samples.reserve(static_cast<std::size_t>(trials));
  for (int trial = 0; trial < trials; ++trial) {
    const auto begin = std::chrono::steady_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition) {
      function();
      asm volatile("" ::: "memory");
    }
    const auto elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    samples.push_back(elapsed / repetitions);
  }
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2];
}

int repetitions(std::size_t n) {
  if (n <= 4096) return 300;
  if (n <= 65536) return 80;
  if (n <= 1048576) return 15;
  return 4;
}

double relativeError(double actual, double expected) {
  return std::fabs(actual - expected) / std::max(1.0, std::fabs(expected));
}

std::vector<ReactionPair> makePairs(std::size_t n) {
  std::mt19937_64 generator(0x4e4552445353ULL);
  std::uniform_real_distribution<double> coordinate(-500, 500), offset(-12, 12);
  std::uniform_real_distribution<double> diffusion(.001, 3), rotation(1e-5, .05),
      radius(.5, 5);
  std::vector<ReactionPair> values(n);
  for (std::size_t k = 0; k < n; ++k) {
    auto& p = values[k];
    p.x1 = coordinate(generator);
    p.y1 = coordinate(generator);
    p.z1 = coordinate(generator);
    p.x2 = p.x1 + offset(generator);
    p.y2 = p.y1 + offset(generator);
    p.z2 = p.z1 + offset(generator);
    p.d1x = diffusion(generator);
    p.d1y = diffusion(generator);
    p.d1z = diffusion(generator);
    p.d2x = diffusion(generator);
    p.d2y = diffusion(generator);
    p.d2z = diffusion(generator);
    if ((k & 7) == 0) p.d1z = 0;
    if ((k & 15) == 0) p.d2z = 0;
    p.dr1z = rotation(generator);
    p.dr2z = rotation(generator);
    const double a = radius(generator), b = radius(generator);
    p.mag1 = a * a;
    p.mag2 = b * b;
    p.bindRadius = radius(generator);
  }
  return values;
}

void makePropagation(std::size_t n, std::vector<Motion>& motions,
                     std::vector<Point>& points) {
  constexpr std::size_t pointsPerComplex = 4;
  std::mt19937_64 generator(0x475055504f43ULL);
  std::uniform_real_distribution<double> coordinate(-500, 500);
  std::normal_distribution<double> movement(0, .2), rotation(0, .03);
  motions.resize((n + pointsPerComplex - 1) / pointsPerComplex);
  points.resize(n);
  for (auto& motion : motions) {
    motion = {coordinate(generator), coordinate(generator), coordinate(generator),
              movement(generator), movement(generator), movement(generator),
              rotation(generator), rotation(generator), rotation(generator)};
  }
  for (std::size_t k = 0; k < n; ++k) {
    const auto complex = static_cast<std::uint32_t>(k / pointsPerComplex);
    points[k] = {motions[complex].cx + 10 * movement(generator),
                 motions[complex].cy + 10 * movement(generator),
                 motions[complex].cz + 10 * movement(generator), complex};
  }
}

std::vector<int> threadCounts() {
  const int maximum = omp_get_num_procs();
  if (const char* requested = std::getenv("OMP_BENCH_THREADS")) {
    std::vector<int> counts;
    std::stringstream stream(requested);
    std::string token;
    while (std::getline(stream, token, ',')) {
      const int count = std::stoi(token);
      if (count < 1 || count > maximum) {
        std::cerr << "OMP_BENCH_THREADS value " << count
                  << " is outside [1," << maximum << "]\n";
        std::exit(2);
      }
      if (counts.empty() || counts.back() != count) counts.push_back(count);
    }
    if (counts.empty()) {
      std::cerr << "OMP_BENCH_THREADS did not contain a thread count\n";
      std::exit(2);
    }
    return counts;
  }
  std::vector<int> counts;
  for (const int candidate : {1, 2, 4, 8, 16, 32, 64}) {
    if (candidate <= maximum) counts.push_back(candidate);
  }
  if (counts.empty() || counts.back() != maximum) counts.push_back(maximum);
  return counts;
}

void printRow(const char* stage, std::size_t n, int threads, int reps,
              int trials, double serialMs, double openMpMs, double maxError,
              std::size_t mismatches) {
  std::cout << stage << ',' << n << ',' << threads << ',' << reps << ','
            << trials << ',' << serialMs << ',' << openMpMs << ','
            << serialMs / openMpMs << ','
            << serialMs / openMpMs / threads << ',' << maxError << ','
            << mismatches << '\n';
}

void benchmarkReaction(std::size_t n, const std::vector<int>& threads,
                       int trials) {
  auto input = makePairs(n);
  std::vector<ReactionResult> serial(n), parallel(n);
  const int reps = repetitions(n);
  reactionSerial(input, serial, .1);
  const double serialMs = medianMilliseconds(
      [&] { reactionSerial(input, serial, .1); }, reps, trials);

  for (const int count : threads) {
    reactionOpenMp(input, parallel, .1, count);
    const double openMpMs = medianMilliseconds(
        [&] { reactionOpenMp(input, parallel, .1, count); }, reps, trials);
    double maxError = 0;
    std::size_t mismatches = 0;
    for (std::size_t k = 0; k < n; ++k) {
      maxError = std::max({maxError,
          relativeError(parallel[k].distance, serial[k].distance),
          relativeError(parallel[k].separation, serial[k].separation),
          relativeError(parallel[k].effectiveDiffusion,
                        serial[k].effectiveDiffusion)});
      mismatches += parallel[k].withinRmax != serial[k].withinRmax;
    }
    printRow("reaction", n, count, reps, trials, serialMs, openMpMs,
             maxError, mismatches);
    if (maxError > 1e-12 || mismatches != 0) std::exit(3);
  }
}

void benchmarkPropagation(std::size_t n, const std::vector<int>& threads,
                          int trials) {
  std::vector<Motion> motions;
  std::vector<Point> input;
  makePropagation(n, motions, input);
  std::vector<Quat> serialQuaternions(motions.size());
  std::vector<Quat> parallelQuaternions(motions.size());
  std::vector<Point> serial(n), parallel(n);
  const int reps = repetitions(n);
  propagationSerial(motions, input, serialQuaternions, serial);
  const double serialMs = medianMilliseconds(
      [&] { propagationSerial(motions, input, serialQuaternions, serial); },
      reps, trials);

  for (const int count : threads) {
    propagationOpenMp(motions, input, parallelQuaternions, parallel, count);
    const double openMpMs = medianMilliseconds(
        [&] { propagationOpenMp(motions, input, parallelQuaternions, parallel,
                                count); },
        reps, trials);
    double maxError = 0;
    std::size_t mismatches = 0;
    for (std::size_t k = 0; k < motions.size(); ++k) {
      maxError = std::max({maxError,
          relativeError(parallelQuaternions[k].x, serialQuaternions[k].x),
          relativeError(parallelQuaternions[k].y, serialQuaternions[k].y),
          relativeError(parallelQuaternions[k].z, serialQuaternions[k].z),
          relativeError(parallelQuaternions[k].w, serialQuaternions[k].w)});
    }
    for (std::size_t k = 0; k < n; ++k) {
      maxError = std::max({maxError,
          relativeError(parallel[k].x, serial[k].x),
          relativeError(parallel[k].y, serial[k].y),
          relativeError(parallel[k].z, serial[k].z)});
      mismatches += parallel[k].complex != serial[k].complex;
    }
    printRow("propagation", n, count, reps, trials, serialMs, openMpMs,
             maxError, mismatches);
    if (maxError > 1e-12 || mismatches != 0) std::exit(3);
  }
}

}  // namespace

int main(int argc, char** argv) {
  omp_set_dynamic(0);
  constexpr int trials = 7;
  const auto threads = threadCounts();
  std::vector<std::size_t> sizes{1024, 4096, 16384, 65536, 262144,
                                 1048576, 4194304};
  if (argc > 1) {
    sizes.clear();
    for (int i = 1; i < argc; ++i) sizes.push_back(std::stoull(argv[i]));
  }

  std::cout << std::setprecision(9)
            << "# omp_num_procs=" << omp_get_num_procs()
            << ",omp_max_threads=" << omp_get_max_threads()
            << ",schedule=static,trials=" << trials
            << ",precision=double,points_per_complex=4\n"
            << "stage,n,threads,repetitions,trials,serial_ms,openmp_ms,"
               "speedup,efficiency,max_relative_error,id_or_flag_mismatches\n";
  for (const auto n : sizes) benchmarkReaction(n, threads, trials);
  for (const auto n : sizes) benchmarkPropagation(n, threads, trials);
}
