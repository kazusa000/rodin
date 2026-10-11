# Distributed under the Boost Software License, Version 1.0.
# See LICENSE in the Rodin repository root.
"""Explicit optional CUDA FP64 matmat; CPU remains the solver/reference."""
from __future__ import annotations

import ctypes
import hashlib
import json
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parents[2] / "cuda/bem_quotient_cuda.cu"
DOUBLE = ctypes.POINTER(ctypes.c_double)
INTEGER = ctypes.POINTER(ctypes.c_int)


class CudaMatmat:
    def __init__(self, *, points, sqrt_areas, rotations, epsilon, viscosity, library, device_index=0):
        if not library:
            raise ValueError("CUDA selection requires an explicit --cuda-library; no CPU fallback")
        path = Path(library).resolve(strict=True)
        self.library = ctypes.CDLL(str(path))
        lib = self.library
        lib.bem_cuda_abi_version.restype = ctypes.c_int
        lib.bem_cuda_source_sha256.restype = ctypes.c_char_p
        lib.bem_cuda_error.restype = ctypes.c_char_p
        lib.bem_cuda_create.restype = ctypes.c_void_p
        lib.bem_cuda_create.argtypes = [ctypes.c_int, ctypes.c_int, DOUBLE, DOUBLE, INTEGER,
                                       ctypes.c_double, ctypes.c_double, ctypes.c_int]
        lib.bem_cuda_apply.restype = ctypes.c_int
        lib.bem_cuda_apply.argtypes = [ctypes.c_void_p, DOUBLE, ctypes.c_int, DOUBLE]
        lib.bem_cuda_metadata.restype = ctypes.c_char_p
        lib.bem_cuda_metadata.argtypes = [ctypes.c_void_p]
        lib.bem_cuda_destroy.argtypes = [ctypes.c_void_p]
        lib.bem_cuda_destroy.restype = None
        if lib.bem_cuda_abi_version() != 1:
            raise RuntimeError("incompatible CUDA operator ABI")
        source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
        if lib.bem_cuda_source_sha256().decode("ascii") != source_hash:
            raise RuntimeError("CUDA library source hash differs from this checkout; rebuild fixed source")
        points = np.ascontiguousarray(points, dtype=np.float64)
        areas = np.ascontiguousarray(sqrt_areas, dtype=np.float64)
        original_rotations = np.asarray(rotations)
        matrices = np.ascontiguousarray(rotations, dtype=np.int32)
        if (points.ndim != 2 or points.shape[1] != 3 or not len(points)
                or areas.shape != (len(points),) or matrices.ndim != 3
                or matrices.shape[1:] != (3, 3) or not len(matrices)):
            raise ValueError("invalid CUDA panel/rotation array shapes")
        if not np.array_equal(original_rotations, matrices):
            raise ValueError("CUDA rotations must be exact signed permutations")
        self.n = len(points)
        self.handle = lib.bem_cuda_create(self.n, len(matrices), points.ctypes.data_as(DOUBLE),
            areas.ctypes.data_as(DOUBLE), matrices.ctypes.data_as(INTEGER), epsilon, viscosity, device_index)
        if not self.handle:
            raise RuntimeError("CUDA operator initialization failed: " + lib.bem_cuda_error().decode())
        self.metadata = json.loads(lib.bem_cuda_metadata(self.handle).decode())
        self.metadata.update({"source_sha256": source_hash, "library": str(path),
                              "library_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})

    def matmat(self, values):
        if values.ndim != 2 or values.shape[0] != 9*self.n or values.shape[1] not in (1, 2):
            raise ValueError("CUDA quotient matmat expects (9*n, 1|2)")
        if values.dtype != np.float64 or not np.all(np.isfinite(values)):
            raise ValueError("CUDA quotient matmat requires finite float64 inputs")
        array = np.ascontiguousarray(values)
        output = np.empty_like(array)
        if self.library.bem_cuda_apply(self.handle, array.ctypes.data_as(DOUBLE), array.shape[1],
                                       output.ctypes.data_as(DOUBLE)):
            raise RuntimeError("CUDA quotient matmat failed: " + self.library.bem_cuda_error().decode())
        if not np.all(np.isfinite(output)):
            raise RuntimeError("CUDA quotient matmat returned nonfinite values")
        return output

    def close(self):
        if getattr(self, "handle", None):
            self.library.bem_cuda_destroy(self.handle)
            self.handle = None

    def __del__(self):
        self.close()


def make_cuda_operator(cpu_operator, *, library, device_index=0):
    """Delegate preconditioning unchanged; override only the operator action."""
    class CudaOperator:
        def __init__(self):
            self.cuda = CudaMatmat(points=cpu_operator.points, sqrt_areas=cpu_operator.sqrt_areas,
                                  rotations=cpu_operator.rotations, epsilon=cpu_operator.epsilon,
                                  viscosity=cpu_operator.viscosity, library=library, device_index=device_index)
            self.n = cpu_operator.n
            self.row_chunk = cpu_operator.row_chunk
            self.applications = 0

        def matmat(self, values):
            self.applications += 1
            return self.cuda.matmat(values)

        def make_preconditioner(self, panel_block_size):
            return cpu_operator.make_preconditioner(panel_block_size)

    return CudaOperator()
