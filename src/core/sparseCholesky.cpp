#include "core/sparseCholesky.h"

#if PP_USE_ACCELERATE
#include <Accelerate/Accelerate.h>
#include <algorithm>
#include <cstdlib>
#include <vector>
#elif PP_USE_CHOLMOD
#include <cholmod.h>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>
#else
#include <Eigen/SparseCholesky>
#endif

namespace PP
{

#if PP_USE_ACCELERATE || PP_USE_CHOLMOD
namespace
{

// The BLAS under the factorization (Accelerate's own, or the one CHOLMOD's
// supernodal factorization calls: Accelerate on macOS, OpenBLAS or MKL on
// Linux) runs its own threads, controlled only by environment variables read
// when the BLAS initializes. Unless they are already set, they are set to 1
// before main() (static initialization of the core library), so that
// factorizations are reproducible; PP_MULTITHREADED_FACTORIZATION=1 in the
// environment keeps the BLAS default (all cores).
struct _SingleThreadedBlas
{
    _SingleThreadedBlas()
    {
        const char* multi = getenv("PP_MULTITHREADED_FACTORIZATION");
        if (not (multi and multi[0] == '1')) {
            setenv("VECLIB_MAXIMUM_THREADS", "1", 0);
            setenv("OPENBLAS_NUM_THREADS", "1", 0);
            setenv("MKL_NUM_THREADS", "1", 0);
        }
    }
};
const _SingleThreadedBlas _singleThreadedBlas;

} // namespace
#endif

#if PP_USE_ACCELERATE

// The lower triangle of A in Accelerate's block CSC arrays (Eigen's default
// storage is CSC too): the block columns, with the block rows >= the column
// and each block in column-major order, and the positions of the entries in
// A's values.
struct SparseCholesky::Impl
{
    int                 _requestedBlockSize = 1;
    int                 _blockSize = 1;  // of the current pattern
    std::vector<long>   _columnStarts;
    std::vector<int>    _rowIndices;
    std::vector<double> _values;
    std::vector<long>   _source;  // index in A.valuePtr() of each entry
    SparseOpaqueSymbolicFactorization _symbolic = {};
    SparseOpaqueFactorization_Double  _numeric = {};
    bool _hasSymbolic = false;
    bool _hasNumeric = false;

    ~Impl() { Clear(); }

    void Clear()
    {
        if (_hasNumeric) {
            SparseCleanup(_numeric);
            _hasNumeric = false;
        }
        if (_hasSymbolic) {
            SparseCleanup(_symbolic);
            _hasSymbolic = false;
        }
    }

    SparseMatrix_Double GetMatrix()
    {
        SparseMatrix_Double M;
        M.structure.rowCount = int(_columnStarts.size()) - 1;  // blocks
        M.structure.columnCount = M.structure.rowCount;
        M.structure.columnStarts = _columnStarts.data();
        M.structure.rowIndices = _rowIndices.data();
        M.structure.attributes = SparseAttributes_t();
        M.structure.attributes.kind = SparseSymmetric;
        M.structure.attributes.triangle = SparseLowerTriangle;
        M.structure.blockSize = uint8_t(_blockSize);
        M.data = _values.data();
        return M;
    }
};

namespace
{

// Errors are reported through info(): no messages on stderr.
void
_IgnoreError(const char*)
{
}

} // namespace

SparseCholesky::SparseCholesky(int blockSize) : _impl(new Impl)
{
    SetBlockSize(blockSize);
}

void
SparseCholesky::SetBlockSize(int blockSize)
{
    _impl->_requestedBlockSize = std::max(1, std::min(blockSize, 255));
}


SparseCholesky::~SparseCholesky() = default;

const char*
SparseCholesky::GetBackend()
{
    return "Accelerate";
}

void
SparseCholesky::analyzePattern(const SparseMatrix& A)
{
    Impl& m = *_impl;
    const int bs = A.cols() % m._requestedBlockSize == 0 ? m._requestedBlockSize : 1;
    const Eigen::Index nb = A.cols() / bs;
    std::vector<long> columnStarts(size_t(nb) + 1, 0);
    std::vector<int> rowIndices;
    std::vector<long> source;
    rowIndices.reserve(m._rowIndices.size());
    source.reserve(m._source.size());
    const int* outer = A.outerIndexPtr();
    const int* inner = A.innerIndexPtr();
    const int* nnz = A.innerNonZeroPtr();  // null when compressed
    auto end = [&](Eigen::Index j) { return nnz ? outer[j] + nnz[j] : outer[j + 1]; };
    // Block rows >= the block column, from the first column of each block
    // column; the other columns must have the same block rows, as dense
    // blocks (otherwise fall back to scalar entries).
    bool blocked = bs > 1;
    for (Eigen::Index a = 0; a < nb and blocked; ++a) {
        const Eigen::Index j0 = a * bs;
        const int count = end(j0) - outer[j0];
        for (int c = 1; c < bs and blocked; ++c) {
            blocked = end(j0 + c) - outer[j0 + c] == count;
        }
        for (int k = outer[j0]; k < end(j0) and blocked; k += bs) {
            const int row = inner[k];
            blocked = row % bs == 0 and k + bs <= end(j0) and inner[k + bs - 1] == row + bs - 1;
            for (int c = 1; c < bs and blocked; ++c) {
                blocked = inner[outer[j0 + c] + (k - outer[j0])] == row;
            }
        }
    }
    const int b = blocked ? bs : 1;
    if (b != bs) {
        columnStarts.assign(size_t(A.cols()) + 1, 0);
    }
    const Eigen::Index nc = A.cols() / b;
    for (Eigen::Index a = 0; a < nc; ++a) {
        const Eigen::Index j0 = a * b;
        for (int k = outer[j0]; k < end(j0); k += b) {
            if (inner[k] / b < a) {
                continue;
            }
            rowIndices.push_back(inner[k] / b);
            // The block, column-major: column c at outer[j0 + c] + offset.
            for (int c = 0; c < b; ++c) {
                for (int r = 0; r < b; ++r) {
                    source.push_back(outer[j0 + c] + (k - outer[j0]) + r);
                }
            }
        }
        columnStarts[size_t(a) + 1] = long(rowIndices.size());
    }

    // The same lower-triangle pattern (e.g. the Newton iterations once the
    // diagram topology settles): keep the symbolic factorization.
    if (m._hasSymbolic and m._symbolic.status == SparseStatusOK and b == m._blockSize and
        columnStarts == m._columnStarts and rowIndices == m._rowIndices) {
        m._source = std::move(source);
        _info = Eigen::Success;
        return;
    }
    m.Clear();
    _info = Eigen::InvalidInput;
    m._blockSize = b;
    m._columnStarts = std::move(columnStarts);
    m._rowIndices = std::move(rowIndices);
    m._source = std::move(source);
    m._values.resize(m._source.size());  // b x b values per block

    SparseSymbolicFactorOptions options = {};
    options.control = SparseDefaultControl;
    options.orderMethod = SparseOrderDefault;
    options.malloc = malloc;
    options.free = free;
    options.reportError = _IgnoreError;
    m._symbolic = SparseFactor(SparseFactorizationCholesky, m.GetMatrix().structure, options);
    m._hasSymbolic = true;
    _info = m._symbolic.status == SparseStatusOK ? Eigen::Success : Eigen::InvalidInput;
}

void
SparseCholesky::factorize(const SparseMatrix& A)
{
    Impl& m = *_impl;
    if (not m._hasSymbolic or m._symbolic.status != SparseStatusOK) {
        _info = Eigen::InvalidInput;
        return;
    }
    const double* values = A.valuePtr();
    for (size_t k = 0; k < m._source.size(); ++k) {
        m._values[k] = values[m._source[k]];
    }
    if (m._hasNumeric) {
        SparseRefactor(m.GetMatrix(), &m._numeric);
    } else {
        m._numeric = SparseFactor(m._symbolic, m.GetMatrix());
        m._hasNumeric = true;
    }
    _info = m._numeric.status == SparseStatusOK ? Eigen::Success : Eigen::NumericalIssue;
}

SparseCholesky::Vector
SparseCholesky::solve(const Vector& b) const
{
    Vector x = b;
    DenseVector_Double xb;
    xb.count = int(x.size());
    xb.data = x.data();
    SparseSolve(_impl->_numeric, xb);
    return x;
}

#elif PP_USE_CHOLMOD

// CHOLMOD reads A in place (Eigen's compressed CSC with int indices), as a
// symmetric matrix of which only the lower triangle is used. The ordering
// (AMD, or METIS when CHOLMOD is built with it) and the choice between the
// simplicial and the supernodal factorization are CHOLMOD's defaults.
struct SparseCholesky::Impl
{
    cholmod_common   _common;
    cholmod_factor*  _factor = nullptr;
    std::vector<int> _outer;  // pattern of the current analysis
    std::vector<int> _inner;

    Impl()
    {
        cholmod_start(&_common);
        _common.print = 0;  // errors are reported through info()
        _common.error_handler = nullptr;
        _common.quick_return_if_not_posdef = 1;
    }

    ~Impl()
    {
        if (_factor) {
            cholmod_free_factor(&_factor, &_common);
        }
        cholmod_finish(&_common);
    }

    // A view of A, which must be compressed.
    static cholmod_sparse View(const SparseMatrix& A)
    {
        cholmod_sparse M = {};
        M.nrow = size_t(A.rows());
        M.ncol = size_t(A.cols());
        M.nzmax = size_t(A.nonZeros());
        M.p = const_cast<int*>(A.outerIndexPtr());
        M.i = const_cast<int*>(A.innerIndexPtr());
        M.x = const_cast<double*>(A.valuePtr());
        M.stype = -1;  // symmetric, lower triangle
        M.itype = CHOLMOD_INT;
        M.xtype = CHOLMOD_REAL;
        M.dtype = CHOLMOD_DOUBLE;
        M.sorted = 1;
        M.packed = 1;
        return M;
    }
};

SparseCholesky::SparseCholesky(int) : _impl(new Impl) {}

void
SparseCholesky::SetBlockSize(int)
{
}

SparseCholesky::~SparseCholesky() = default;

const char*
SparseCholesky::GetBackend()
{
    return "CHOLMOD";
}

void
SparseCholesky::analyzePattern(const SparseMatrix& A)
{
    Impl& m = *_impl;
    if (not A.isCompressed()) {
        SparseMatrix C = A;
        C.makeCompressed();
        analyzePattern(C);
        return;
    }
    const int* outer = A.outerIndexPtr();
    const int* inner = A.innerIndexPtr();
    const size_t nc = size_t(A.cols()) + 1;
    const size_t nz = size_t(A.nonZeros());
    // The same pattern (e.g. the Newton iterations once the diagram topology
    // settles): keep the symbolic factorization.
    if (m._factor and m._outer.size() == nc and m._inner.size() == nz and
        std::memcmp(m._outer.data(), outer, nc * sizeof(int)) == 0 and
        std::memcmp(m._inner.data(), inner, nz * sizeof(int)) == 0) {
        _info = Eigen::Success;
        return;
    }
    if (m._factor) {
        cholmod_free_factor(&m._factor, &m._common);
    }
    m._outer.assign(outer, outer + nc);
    m._inner.assign(inner, inner + nz);
    cholmod_sparse M = Impl::View(A);
    m._factor = cholmod_analyze(&M, &m._common);
    _info = m._factor ? Eigen::Success : Eigen::InvalidInput;
}

void
SparseCholesky::factorize(const SparseMatrix& A)
{
    Impl& m = *_impl;
    if (not m._factor) {
        _info = Eigen::InvalidInput;
        return;
    }
    if (not A.isCompressed()) {
        SparseMatrix C = A;
        C.makeCompressed();
        factorize(C);
        return;
    }
    cholmod_sparse M = Impl::View(A);
    const int ok = cholmod_factorize(&M, m._factor, &m._common);
    _info = ok and m._factor->minor == m._factor->n ? Eigen::Success : Eigen::NumericalIssue;
}

SparseCholesky::Vector
SparseCholesky::solve(const Vector& b) const
{
    Impl& m = *_impl;
    cholmod_dense B = {};
    B.nrow = size_t(b.size());
    B.ncol = 1;
    B.nzmax = B.nrow;
    B.d = B.nrow;
    B.x = const_cast<double*>(b.data());
    B.xtype = CHOLMOD_REAL;
    B.dtype = CHOLMOD_DOUBLE;
    cholmod_dense* X = cholmod_solve(CHOLMOD_A, m._factor, &B, &m._common);
    Vector x(b.size());
    if (X) {
        std::memcpy(x.data(), X->x, size_t(b.size()) * sizeof(double));
        cholmod_free_dense(&X, &m._common);
    } else {
        x.setConstant(std::numeric_limits<double>::quiet_NaN());
    }
    return x;
}

#else // Eigen

struct SparseCholesky::Impl
{
    Eigen::SimplicialLLT<SparseMatrix, Eigen::Lower> _llt;
};

SparseCholesky::SparseCholesky(int) : _impl(new Impl) {}

void
SparseCholesky::SetBlockSize(int)
{
}


SparseCholesky::~SparseCholesky() = default;

const char*
SparseCholesky::GetBackend()
{
    return "Eigen";
}

void
SparseCholesky::analyzePattern(const SparseMatrix& A)
{
    _impl->_llt.analyzePattern(A);
    _info = _impl->_llt.info();
}

void
SparseCholesky::factorize(const SparseMatrix& A)
{
    _impl->_llt.factorize(A);
    _info = _impl->_llt.info();
}

SparseCholesky::Vector
SparseCholesky::solve(const Vector& b) const
{
    return _impl->_llt.solve(b);
}

#endif

} // namespace PP
