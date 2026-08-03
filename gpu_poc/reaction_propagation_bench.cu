#include <cuda_runtime.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#define CUDA_CHECK(x) do { cudaError_t e=(x); if(e!=cudaSuccess){ \
  std::cerr<<"CUDA error: "<<cudaGetErrorString(e)<<" at "<<__LINE__<<"\n"; \
  std::exit(2); } } while(0)

namespace {

struct ReactionPair {
  double x1,y1,z1,x2,y2,z2;
  double d1x,d1y,d1z,d2x,d2y,d2z;
  double dr1z,dr2z,mag1,mag2,bindRadius;
};
struct ReactionResult {
  double distance,separation,effectiveDiffusion;
  unsigned char withinRmax;
};
struct Motion { double cx,cy,cz,tx,ty,tz,rx,ry,rz; };
struct Quat { double x,y,z,w; };
struct Point { double x,y,z; uint32_t complex; };

__host__ __device__
ReactionResult checkReaction(const ReactionPair& p,double dt) {
  double d=(p.d1x+p.d2x+p.d1y+p.d2y+p.d1z+p.d2z)/3.0;
  bool flat=fabs(p.d1z)<1e-10;
  d+=2*p.mag1*(1-cos(sqrt((flat?2.0:4.0)*p.dr1z*dt)))
      /((flat?4.0:6.0)*dt);
  flat=fabs(p.d2z)<1e-10;
  d+=2*p.mag2*(1-cos(sqrt((flat?2.0:4.0)*p.dr2z*dt)))
      /((flat?4.0:6.0)*dt);
  double rmax=3*sqrt(6*d*dt)+p.bindRadius;
  double dx=p.x1-p.x2,dy=p.y1-p.y2,dz=p.z1-p.z2;
  double r=sqrt(dx*dx+dy*dy+dz*dz);
  return {r,r-p.bindRadius,d,(unsigned char)(r<rmax)};
}

__host__ __device__
Quat makeQuat(const Motion& m) {
  double cz=cos(m.rz*.5),sz=sin(m.rz*.5);
  double cy=cos(m.ry*.5),sy=sin(m.ry*.5);
  double cx=cos(m.rx*.5),sx=sin(m.rx*.5);
  Quat q{sx*cy*cz-cx*sy*sz,cx*sy*cz+sx*cy*sz,
         cx*cy*sz-sx*sy*cz,cx*cy*cz+sx*sy*sz};
  double s=1/sqrt(q.x*q.x+q.y*q.y+q.z*q.z+q.w*q.w);
  q.x*=s;q.y*=s;q.z*=s;q.w*=s;
  return q;
}

__host__ __device__
Point propagate(const Point& p,const Motion& m,const Quat& q) {
  double x=p.x-m.cx,y=p.y-m.cy,z=p.z-m.cz;
  double tx=2*(q.y*z-q.z*y),ty=2*(q.z*x-q.x*z),tz=2*(q.x*y-q.y*x);
  return {
    x+q.w*tx+q.y*tz-q.z*ty+m.cx+m.tx,
    y+q.w*ty+q.z*tx-q.x*tz+m.cy+m.ty,
    z+q.w*tz+q.x*ty-q.y*tx+m.cz+m.tz,p.complex
  };
}

__global__ void reactionKernel(const ReactionPair* i,ReactionResult* o,
                               size_t n,double dt) {
  size_t k=blockIdx.x*blockDim.x+threadIdx.x;
  if(k<n)o[k]=checkReaction(i[k],dt);
}
__global__ void quatKernel(const Motion* i,Quat* o,size_t n) {
  size_t k=blockIdx.x*blockDim.x+threadIdx.x;
  if(k<n)o[k]=makeQuat(i[k]);
}
__global__ void propagationKernel(const Point* i,const Motion* m,
                                  const Quat* q,Point* o,size_t n) {
  size_t k=blockIdx.x*blockDim.x+threadIdx.x;
  if(k<n)o[k]=propagate(i[k],m[i[k].complex],q[i[k].complex]);
}

template<class F> double cpuMs(F&& f,int reps) {
  auto begin=std::chrono::steady_clock::now();
  for(int i=0;i<reps;++i) {
    f();
    asm volatile("" ::: "memory");
  }
  return std::chrono::duration<double,std::milli>(
    std::chrono::steady_clock::now()-begin).count()/reps;
}
template<class F> double gpuMs(F&& f,int reps) {
  cudaEvent_t begin,end;
  CUDA_CHECK(cudaEventCreate(&begin));
  CUDA_CHECK(cudaEventCreate(&end));
  CUDA_CHECK(cudaEventRecord(begin));
  for(int i=0;i<reps;++i)f();
  CUDA_CHECK(cudaEventRecord(end));
  CUDA_CHECK(cudaEventSynchronize(end));
  float ms=0;
  CUDA_CHECK(cudaEventElapsedTime(&ms,begin,end));
  CUDA_CHECK(cudaEventDestroy(begin));
  CUDA_CHECK(cudaEventDestroy(end));
  return ms/reps;
}
int repetitions(size_t n) {
  return n<=4096?200:n<=65536?100:n<=1048576?30:10;
}
double relativeError(double a,double b) {
  return fabs(a-b)/std::max(1.0,fabs(b));
}

std::vector<ReactionPair> makePairs(size_t n) {
  std::mt19937_64 g(0x4e4552445353ULL);
  std::uniform_real_distribution<double> c(-500,500),dx(-12,12);
  std::uniform_real_distribution<double> d(.001,3),dr(1e-5,.05),r(.5,5);
  std::vector<ReactionPair> v(n);
  for(size_t k=0;k<n;++k) {
    auto& p=v[k];
    p.x1=c(g);p.y1=c(g);p.z1=c(g);
    p.x2=p.x1+dx(g);p.y2=p.y1+dx(g);p.z2=p.z1+dx(g);
    p.d1x=d(g);p.d1y=d(g);p.d1z=d(g);
    p.d2x=d(g);p.d2y=d(g);p.d2z=d(g);
    if((k&7)==0)p.d1z=0;
    if((k&15)==0)p.d2z=0;
    p.dr1z=dr(g);p.dr2z=dr(g);
    double a=r(g),b=r(g);
    p.mag1=a*a;p.mag2=b*b;p.bindRadius=r(g);
  }
  return v;
}

void makePropagation(size_t n,std::vector<Motion>& m,std::vector<Point>& p) {
  constexpr size_t pointsPerComplex=4;
  std::mt19937_64 g(0x475055504f43ULL);
  std::uniform_real_distribution<double> c(-500,500);
  std::normal_distribution<double> move(0,.2),rot(0,.03);
  m.resize((n+pointsPerComplex-1)/pointsPerComplex);
  p.resize(n);
  for(auto& v:m)
    v={c(g),c(g),c(g),move(g),move(g),move(g),rot(g),rot(g),rot(g)};
  for(size_t k=0;k<n;++k) {
    uint32_t z=k/pointsPerComplex;
    p[k]={m[z].cx+10*move(g),m[z].cy+10*move(g),
          m[z].cz+10*move(g),z};
  }
}

void benchReaction(size_t n) {
  auto input=makePairs(n);
  std::vector<ReactionResult> cpu(n),gpu(n);
  ReactionPair* di;
  ReactionResult* dout;
  CUDA_CHECK(cudaMalloc(&di,n*sizeof(*di)));
  CUDA_CHECK(cudaMalloc(&dout,n*sizeof(*dout)));
  CUDA_CHECK(cudaMemcpy(di,input.data(),n*sizeof(*di),cudaMemcpyHostToDevice));
  int blocks=(n+255)/256,reps=repetitions(n);
  reactionKernel<<<blocks,256>>>(di,dout,n,.1);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  double c=cpuMs([&] {
    for(size_t k=0;k<n;++k)cpu[k]=checkReaction(input[k],.1);
  },reps);
  double resident=gpuMs([&] {
    reactionKernel<<<blocks,256>>>(di,dout,n,.1);
    CUDA_CHECK(cudaGetLastError());
  },reps);
  double endToEnd=gpuMs([&] {
    CUDA_CHECK(cudaMemcpyAsync(di,input.data(),n*sizeof(*di),
                               cudaMemcpyHostToDevice));
    reactionKernel<<<blocks,256>>>(di,dout,n,.1);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(gpu.data(),dout,n*sizeof(*dout),
                               cudaMemcpyDeviceToHost));
  },reps);
  CUDA_CHECK(cudaMemcpy(gpu.data(),dout,n*sizeof(*dout),
                        cudaMemcpyDeviceToHost));

  double error=0;size_t mismatch=0;
  for(size_t k=0;k<n;++k) {
    error=std::max({error,
      relativeError(gpu[k].distance,cpu[k].distance),
      relativeError(gpu[k].separation,cpu[k].separation),
      relativeError(gpu[k].effectiveDiffusion,cpu[k].effectiveDiffusion)});
    mismatch+=gpu[k].withinRmax!=cpu[k].withinRmax;
  }
  if(error>1e-12 || mismatch) {
    std::cerr<<"reaction correctness gate failed"<<std::endl;
    std::exit(3);
  }
  std::cout<<"reaction,"<<n<<','<<reps<<','<<c<<','<<resident<<','
           <<endToEnd<<','<<c/resident<<','<<c/endToEnd<<','
           <<error<<','<<mismatch<<'\n';
  CUDA_CHECK(cudaFree(di));CUDA_CHECK(cudaFree(dout));
}

void benchPropagation(size_t n) {
  std::vector<Motion> motions;
  std::vector<Point> input;
  makePropagation(n,motions,input);
  std::vector<Quat> cpuQuat(motions.size());
  std::vector<Point> cpu(n),gpu(n);
  Motion* dm;Quat* dq;Point *di,*dout;
  CUDA_CHECK(cudaMalloc(&dm,motions.size()*sizeof(*dm)));
  CUDA_CHECK(cudaMalloc(&dq,motions.size()*sizeof(*dq)));
  CUDA_CHECK(cudaMalloc(&di,n*sizeof(*di)));
  CUDA_CHECK(cudaMalloc(&dout,n*sizeof(*dout)));
  CUDA_CHECK(cudaMemcpy(dm,motions.data(),motions.size()*sizeof(*dm),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(di,input.data(),n*sizeof(*di),cudaMemcpyHostToDevice));
  int mb=(motions.size()+255)/256,pb=(n+255)/256,reps=repetitions(n);
  quatKernel<<<mb,256>>>(dm,dq,motions.size());
  CUDA_CHECK(cudaGetLastError());
  propagationKernel<<<pb,256>>>(di,dm,dq,dout,n);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  double c=cpuMs([&] {
    for(size_t k=0;k<motions.size();++k)cpuQuat[k]=makeQuat(motions[k]);
    for(size_t k=0;k<n;++k)
      cpu[k]=propagate(input[k],motions[input[k].complex],
                       cpuQuat[input[k].complex]);
  },reps);
  double resident=gpuMs([&] {
    quatKernel<<<mb,256>>>(dm,dq,motions.size());
    CUDA_CHECK(cudaGetLastError());
    propagationKernel<<<pb,256>>>(di,dm,dq,dout,n);
    CUDA_CHECK(cudaGetLastError());
  },reps);
  double endToEnd=gpuMs([&] {
    CUDA_CHECK(cudaMemcpyAsync(dm,motions.data(),motions.size()*sizeof(*dm),
                               cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyAsync(di,input.data(),n*sizeof(*di),
                               cudaMemcpyHostToDevice));
    quatKernel<<<mb,256>>>(dm,dq,motions.size());
    CUDA_CHECK(cudaGetLastError());
    propagationKernel<<<pb,256>>>(di,dm,dq,dout,n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(gpu.data(),dout,n*sizeof(*dout),
                               cudaMemcpyDeviceToHost));
  },reps);
  CUDA_CHECK(cudaMemcpy(gpu.data(),dout,n*sizeof(*dout),
                        cudaMemcpyDeviceToHost));

  double error=0;size_t mismatch=0;
  for(size_t k=0;k<n;++k) {
    error=std::max({error,relativeError(gpu[k].x,cpu[k].x),
                    relativeError(gpu[k].y,cpu[k].y),
                    relativeError(gpu[k].z,cpu[k].z)});
    mismatch+=gpu[k].complex!=cpu[k].complex;
  }
  if(error>1e-12 || mismatch) {
    std::cerr<<"propagation correctness gate failed"<<std::endl;
    std::exit(3);
  }
  std::cout<<"propagation,"<<n<<','<<reps<<','<<c<<','<<resident<<','
           <<endToEnd<<','<<c/resident<<','<<c/endToEnd<<','
           <<error<<','<<mismatch<<'\n';
  CUDA_CHECK(cudaFree(dm));CUDA_CHECK(cudaFree(dq));
  CUDA_CHECK(cudaFree(di));CUDA_CHECK(cudaFree(dout));
}

} // namespace

int main(int argc,char** argv) {
  CUDA_CHECK(cudaSetDevice(0));
  cudaDeviceProp p{};
  CUDA_CHECK(cudaGetDeviceProperties(&p,0));
  std::cout<<"# gpu="<<p.name<<",cc="<<p.major<<'.'<<p.minor
           <<",memory_bytes="<<p.totalGlobalMem
           <<"\n# warmups=1,repetitions=adaptive,precision=double,"
             "points_per_complex=4\n"
           <<"stage,n,repetitions,cpu_ms,gpu_resident_ms,"
             "gpu_end_to_end_ms,resident_speedup,end_to_end_speedup,"
             "max_relative_error,flag_or_id_mismatches\n";
  std::vector<size_t> sizes{
    1024,4096,16384,65536,262144,1048576,4194304
  };
  if(argc>1) {
    sizes.clear();
    for(int i=1;i<argc;++i)sizes.push_back(std::stoull(argv[i]));
  }
  for(auto n:sizes)benchReaction(n);
  for(auto n:sizes)benchPropagation(n);
}
