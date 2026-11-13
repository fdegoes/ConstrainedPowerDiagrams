#ifndef PP_CORE_SPARSECHOLESKY_H
#define PP_CORE_SPARSECHOLESKY_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <memory>

namespace PP
{

// Sparse Cholesky factorization of a symmetric positive definite matrix,
// with the interface of Eigen's SimplicialLLT: one symbolic analysis per
// sparsity pattern, then numeric factorizations of matrices with that
// pattern. Only the lower triangle of the matrix is read.
//
// The backend is chosen at build time (CMake option PP_SPARSE_SOLVER): Apple
// Accelerate's supernodal Cholesky (macOS), SuiteSparse's CHOLMOD (macOS and
// Linux), or Eigen's SimplicialLLT (portable, slowest). All are reproducible
// run to run: the BLAS of Accelerate and CHOLMOD is limited to one thread
// (multithreaded factorizations change the round-off from run to run), unless
// PP_MULTITHREADED_FACTORIZATION=1 is set in the environment.
class SparseCholesky
{
public:
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using Vector = Eigen::VectorXd;

    // blockSize: the pattern of the matrices is made of dense blockSize x
    // blockSize blocks (e.g. the d + 1 unknowns of a site), which the
    // Accelerate backend factorizes as blocks. Must divide the matrix size.
    explicit SparseCholesky(int blockSize = 1);
    ~SparseCholesky();
    // Block size of the next analyzePattern().
    void SetBlockSize(int blockSize);

    SparseCholesky(const SparseCholesky&) = delete;
    SparseCholesky& operator=(const SparseCholesky&) = delete;

    // Symbolic analysis (ordering and elimination tree) of the pattern of A.
    // Invalidates the factorization, except that the Accelerate backend keeps
    // it when the pattern is that of the previous analysis.
    void analyzePattern(const SparseMatrix& A);
    // Numeric factorization of A, which must have the pattern given to
    // analyzePattern().
    void factorize(const SparseMatrix& A);
    // analyzePattern() and factorize().
    void compute(const SparseMatrix& A)
    {
        analyzePattern(A);
        factorize(A);
    }

    // Success after a successful factorize(), NumericalIssue when A is not
    // positive definite.
    Eigen::ComputationInfo info() const { return _info; }

    // A^-1 b, from the last successful factorization.
    Vector solve(const Vector& b) const;

    // Name of the backend ("Accelerate", "CHOLMOD" or "Eigen").
    static const char* GetBackend();

private:
    struct Impl;
    std::unique_ptr<Impl>  _impl;
    Eigen::ComputationInfo _info = Eigen::InvalidInput;
};

} // namespace PP

#endif // PP_CORE_SPARSECHOLESKY_H
