// Optional FP64 matmat only. No BLAS, solver, sampling, or quadrature changes.
// Distributed under the Boost Software License, Version 1.0.
// See LICENSE in the Rodin repository root.
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef BEM_CUDA_SOURCE_SHA256
#error "Build with scripts/build_bem_cuda.py to bind the source digest"
#endif

namespace {
thread_local std::string error_message;
void check(cudaError_t status, const char* stage) {
  if (status != cudaSuccess)
    throw std::runtime_error(std::string(stage) + ": " + cudaGetErrorString(status));
}
struct Operator {
  int n = 0, groups = 0, device = 0;
  double epsilon2 = 0, coefficient = 0;
  double *points = nullptr, *sqrt_area = nullptr, *input = nullptr;
  double *packed = nullptr, *output = nullptr;
  int* rotations = nullptr;
  std::string metadata;
  ~Operator() {
    cudaSetDevice(device);
    cudaFree(points); cudaFree(sqrt_area); cudaFree(input);
    cudaFree(packed); cudaFree(output); cudaFree(rotations);
  }
};

template<int RHS>
__global__ void pack_input(int n, const double* input, const double* sqrt_area,
                           double* packed) {
  int entry = blockIdx.x * blockDim.x + threadIdx.x;
  if (entry < n * 9 * RHS) {
    int j = entry / (9 * RHS), component = entry % (9 * RHS);
    packed[component * n + j] = sqrt_area[j] * input[entry];
  }
}

template<int RHS>
__global__ void apply_operator(int n, int groups, const double* points,
                              const double* sqrt_area, const int* rotations,
                              double epsilon2, double coefficient,
                              const double* packed, double* output) {
  constexpr int SIZE = 9 * RHS;
  constexpr int WARPS = 4; // launch exactly 128 threads per output panel
  __shared__ double partial[WARPS][SIZE];
  const int i = blockIdx.x, lane = threadIdx.x % 32, warp = threadIdx.x / 32;
  double sum[SIZE] = {};
  const double xi0 = points[i], xi1 = points[n+i], xi2 = points[2*n+i];
  for (int group = 0; group < groups; ++group) {
    const int* g = rotations + 6 * group;
    for (int j = threadIdx.x; j < n; j += blockDim.x) {
      const double r0 = xi0 - g[3] * points[g[0]*n+j];
      const double r1 = xi1 - g[4] * points[g[1]*n+j];
      const double r2 = xi2 - g[5] * points[g[2]*n+j];
      const double distance2 = r0*r0 + r1*r1 + r2*r2;
      const double denominator = distance2 + epsilon2;
      // FP64 sqrt/division, without --use_fast_math or reduced precision.
      const double beta = coefficient / (denominator * sqrt(denominator));
      const double alpha = (distance2 + 2.0 * epsilon2) * beta;
      #pragma unroll
      for (int load = 0; load < 3; ++load) {
        #pragma unroll
        for (int rhs = 0; rhs < RHS; ++rhs) {
          const double f0 = (g[3] * g[3+load]) * packed[(g[0]*3*RHS+g[load]*RHS+rhs)*n+j];
          const double f1 = (g[4] * g[3+load]) * packed[(g[1]*3*RHS+g[load]*RHS+rhs)*n+j];
          const double f2 = (g[5] * g[3+load]) * packed[(g[2]*3*RHS+g[load]*RHS+rhs)*n+j];
          const double dot = beta * (r0*f0 + r1*f1 + r2*f2);
          sum[load*RHS+rhs] += alpha*f0 + r0*dot;
          sum[3*RHS+load*RHS+rhs] += alpha*f1 + r1*dot;
          sum[6*RHS+load*RHS+rhs] += alpha*f2 + r2*dot;
        }
      }
    }
  }
  #pragma unroll
  for (int component = 0; component < SIZE; ++component) {
    double value = sum[component];
    for (int offset = 16; offset; offset /= 2)
      value += __shfl_down_sync(0xffffffff, value, offset);
    if (lane == 0) partial[warp][component] = value;
  }
  __syncthreads();
  if (warp == 0) {
    #pragma unroll
    for (int component = 0; component < SIZE; ++component) {
      double value = lane < WARPS ? partial[lane][component] : 0.0;
      for (int offset = 16; offset; offset /= 2)
        value += __shfl_down_sync(0xffffffff, value, offset);
      if (lane == 0) output[i*SIZE+component] = sqrt_area[i] * value;
    }
  }
}

template<int RHS> void launch(Operator* op) {
  pack_input<RHS><<<(op->n*9*RHS+255)/256,256>>>(op->n, op->input, op->sqrt_area, op->packed);
  check(cudaGetLastError(), "pack-input launch");
  apply_operator<RHS><<<op->n,128>>>(op->n, op->groups, op->points, op->sqrt_area,
    op->rotations, op->epsilon2, op->coefficient, op->packed, op->output);
  check(cudaGetLastError(), "FP64 quotient matmat launch");
}
}

extern "C" const char* bem_cuda_error() { return error_message.c_str(); }
extern "C" const char* bem_cuda_source_sha256() { return BEM_CUDA_SOURCE_SHA256; }
extern "C" int bem_cuda_abi_version() { return 1; }

extern "C" void* bem_cuda_create(int n, int groups, const double* points,
    const double* sqrt_area, const int* matrices, double epsilon, double viscosity,
    int device) {
  Operator* op = nullptr;
  try {
    if (n <= 0 || groups <= 0 || !std::isfinite(epsilon) || epsilon <= 0 ||
        !std::isfinite(viscosity) || viscosity <= 0)
      throw std::runtime_error("invalid CUDA operator dimensions or physical parameters");
    op = new Operator;
    op->n = n; op->groups = groups; op->device = device;
    op->epsilon2 = epsilon*epsilon;
    op->coefficient = 1.0 / (8.0 * std::acos(-1.0) * viscosity);
    check(cudaSetDevice(device), "select CUDA device");
    cudaDeviceProp properties{};
    check(cudaGetDeviceProperties(&properties, device), "query CUDA device");
    if (properties.major < 7 || (properties.major == 7 && properties.minor < 5))
      throw std::runtime_error("CUDA backend requires compute capability >=7.5");
    std::vector<double> transposed(3*n);
    for (int j = 0; j < n; ++j) {
      if (!std::isfinite(sqrt_area[j]) || sqrt_area[j] <= 0)
        throw std::runtime_error("CUDA panel weights must be finite and positive");
      for (int component = 0; component < 3; ++component) {
        if (!std::isfinite(points[3*j+component])) throw std::runtime_error("nonfinite CUDA panel centroid");
        transposed[component*n+j] = points[3*j+component];
      }
    }
    std::vector<int> permutations(6*groups);
    for (int group = 0; group < groups; ++group) {
      int columns[3] = {};
      for (int row = 0; row < 3; ++row) {
        int count = 0;
        for (int col = 0; col < 3; ++col) {
          int value = matrices[9*group+3*row+col];
          if (value) {
            if (value != 1 && value != -1) throw std::runtime_error("rotation is not a signed permutation");
            permutations[6*group+row] = col; permutations[6*group+3+row] = value;
            ++count; ++columns[col];
          }
        }
        if (count != 1) throw std::runtime_error("rotation row is not a signed permutation");
      }
      if (columns[0] != 1 || columns[1] != 1 || columns[2] != 1)
        throw std::runtime_error("rotation columns are not a permutation");
    }
    check(cudaMalloc(&op->points, 3*n*sizeof(double)), "allocate centroids");
    check(cudaMalloc(&op->sqrt_area, n*sizeof(double)), "allocate panel weights");
    check(cudaMalloc(&op->rotations, 6*groups*sizeof(int)), "allocate rotations");
    check(cudaMalloc(&op->input, 18*n*sizeof(double)), "allocate input");
    check(cudaMalloc(&op->packed, 18*n*sizeof(double)), "allocate packed response");
    check(cudaMalloc(&op->output, 18*n*sizeof(double)), "allocate output");
    check(cudaMemcpy(op->points, transposed.data(), 3*n*sizeof(double), cudaMemcpyHostToDevice), "copy centroids");
    check(cudaMemcpy(op->sqrt_area, sqrt_area, n*sizeof(double), cudaMemcpyHostToDevice), "copy panel weights");
    check(cudaMemcpy(op->rotations, permutations.data(), 6*groups*sizeof(int), cudaMemcpyHostToDevice), "copy rotations");
    int runtime = 0, driver = 0;
    check(cudaRuntimeGetVersion(&runtime), "query CUDA runtime");
    check(cudaDriverGetVersion(&driver), "query CUDA driver");
    char uuid[33] = {}, data[1024] = {};
    for (int i = 0; i < 16; ++i) std::snprintf(uuid+2*i, 3, "%02x", static_cast<unsigned char>(properties.uuid.bytes[i]));
    std::snprintf(data, sizeof(data),
      "{\"device\":\"cuda\",\"dtype\":\"float64\",\"device_index\":%d,"
      "\"device_name\":\"%s\",\"uuid_hex\":\"%s\",\"compute_capability\":\"%d.%d\","
      "\"runtime_version\":%d,\"driver_version\":%d,\"device_total_memory_bytes\":%llu,"
      "\"operator_device_allocation_bytes\":%llu,\"fast_math\":false,\"fmad\":false}",
      device, properties.name, uuid, properties.major, properties.minor, runtime, driver,
      static_cast<unsigned long long>(properties.totalGlobalMem),
      static_cast<unsigned long long>(58*n*sizeof(double)+6*groups*sizeof(int)));
    op->metadata = data;
    return op;
  } catch (const std::exception& error) {
    error_message = error.what(); delete op; return nullptr;
  }
}

extern "C" const char* bem_cuda_metadata(void* pointer) {
  return static_cast<Operator*>(pointer)->metadata.c_str();
}

extern "C" int bem_cuda_apply(void* pointer, const double* input, int rhs, double* output) {
  try {
    auto* op = static_cast<Operator*>(pointer);
    if (!op || (rhs != 1 && rhs != 2)) throw std::runtime_error("CUDA matmat requires one or two RHS");
    check(cudaSetDevice(op->device), "restore CUDA device");
    size_t bytes = size_t(op->n)*9*rhs*sizeof(double);
    check(cudaMemcpy(op->input, input, bytes, cudaMemcpyHostToDevice), "copy matmat input");
    if (rhs == 1) launch<1>(op); else launch<2>(op);
    check(cudaDeviceSynchronize(), "synchronize FP64 matmat");
    check(cudaMemcpy(output, op->output, bytes, cudaMemcpyDeviceToHost), "copy matmat output");
    return 0;
  } catch (const std::exception& error) { error_message = error.what(); return 1; }
}

extern "C" void bem_cuda_destroy(void* pointer) { delete static_cast<Operator*>(pointer); }
