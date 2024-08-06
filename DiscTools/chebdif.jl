# Chebyshev Collocation Method: Directly ported from MATLAB code by JAC Weideman and SC Reddy
# %  The function [x, DM] =  chebdif(N,M) computes the differentiation 
# %  matrices D1, D2, ..., DM on Chebyshev nodes. 
# % 
# %  Input:
# %  N:        Size of differentiation matrix.        
# %  M:        Number of derivatives required (integer).
# %  Note:     0 < M <= N-1.
# %
# %  Output:
# %  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
# %
# %  The code implements two strategies for enhanced 
# %  accuracy suggested by W. Don and S. Solomonoff in 
# %  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
# %  The two strategies are (a) the use of trigonometric 
# %  identities to avoid the computation of differences 
# %  x(k)-x(j) and (b) the use of the "flipping trick"
# %  which is necessary since sin t can be computed to high
# %  relative precision when t is small whereas sin (pi-t) cannot.
# %  Note added May 2003:  It may, in fact, be slightly better not to
# %  implement the strategies (a) and (b).   Please consult the following
# %  paper for details:   "Spectral Differencing with a Twist", by
# %  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp. 

# %  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by 
# %  JACW, May 2003.
using LinearAlgebra:I as I
using LinearAlgebra:diagind as diagind
using LinearAlgebra:diag as diag

function chebdif(N, M)
   n1 = floor(Int64,N/2)
   n2 = ceil(Int64, N/2)

   k = collect(0:N-1)  # Compute theta vector.
   th = k .* pi / (N-1)
   x = sin.(pi .* (N-1:-2:1-N) ./ (2*(N-1)))  # Compute Chebyshev points.
   T = repeat(th/2, 1, N)
   DX = 2 .* sin.(T' .+ T) .* sin.(T' .- T)  # Trigonometric identity.
   DX = [DX[1:n1, :]; -reverse(DX[1:n2, :])]  # Flipping trick.
   DX[diagind(DX)] .= 1.0  # Put 1's on the main diagonal of DX.

   C = [(-1.0)^(i+j) for i in 0:N-1, j in 0:N-1]  # C is the matrix with entries c(k)/c(j)
   C[1, :] .= C[1, :] .* 2.0
   C[N, :] .= C[N, :] .* 2.0
   C[:, 1] .= C[:, 1] ./ 2.0
   C[:, N] .= C[:, N] ./ 2.0

   Z = 1 ./ DX  # Z contains entries 1/(x(k)-x(j))
   Z[diagind(Z)] .= 0.0  # with zeros on the diagonal.
   D = Matrix{Float64}(I,N,N)  # D contains diff. matrices.

   DM = zeros(N, N, M)  # Initialize DM matrix.

   for ell = 1:M
      D = ell .* Z .* (C .* repeat(diag(D), 1, N) .- D)  # Off-diagonals
      D[diagind(D)] .= -sum(D, dims=2)  # Correct main diagonal of D
      DM[:, :, ell] = D  # Store current D in DM
   end

   return x, DM
end

## Chebyshev Interpolation Method: Directly ported from MATLAB code by JAC Weideman and SC Reddy
# %  The function p = chebint(fk, x) computes the polynomial interpolant
# %  of the data (xk, fk), where xk are the Chebyshev nodes.  
# %  Two or more data points are assumed.
# %
# %  Input:
# %  fk:  Vector of y-coordinates of data, at Chebyshev points 
# %       x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
# %  x:   Vector of x-values where polynomial interpolant is to be evaluated.
# %
# %  Output:
# %  p:    Vector of interpolated values.
# %
# %  The code implements the barycentric formula; see page 252 in
# %  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
# %  (Note that if some fk > 1/eps, with eps the machine epsilon,
# %  the value of eps in the code may have to be reduced.)

# %  J.A.C. Weideman, S.C. Reddy 1998
function chebint(fk, x)
   N = length(fk)
   M = length(x)
   
   xk = sin.(pi .* (N-1:-2:1-N) ./ (2*(N-1)))
   
   w = ones(N) .* (-1).^(0:N-1)
   w[1] = w[1] / 2
   w[N] = w[N] / 2
   
   D = x .- xk'
   D = 1 ./ (D .+ eps().*(D .== 0))
   
   p = D * (w .* fk) ./ (D * w)
   
   return p
end