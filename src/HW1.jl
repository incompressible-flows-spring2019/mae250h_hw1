module HW1

# We need FFTW to do the circulant tridiagonal solution
using FFTW

# The next block is the documentation associated with our function. It is
# available if we type ?HW1.scale_it at a Julia prompt.
"""
    scale_it(f::Function,x::Real,c::Real) -> Real

Takes as input a function `f` and two real-valued arguments `x` and `c` and returns
value `c*f(x)`.

For example,
```
julia> g(x) = x^2

julia> HW1.scale_it(g,1,2)
2
```
"""
function scale_it(f::Function,x::Real,c::Real)
  return c*f(x)
end

"""
    numerical_integrate(f::Function,a::Real,b::Real,scheme::String;N::Int=100) -> Real

Numerically integrate a function `f` from `a` to `b`, using quadrature of the
specified by `scheme`. This can be either `"trapezoidal"` or `"simpson"`.
The interval [a,b] is subdivided into `N` uniform intervals; `N` defaults to 100. The
result is returned as a real value.
"""
function numerical_integrate(f::Function,a::Real,b::Real,scheme::String;N::Int = 100)

  dx = (b-a)/N
  if scheme == "trapezoidal"
    res = 0.5*(f(a)+f(b))
    x = a
    for i = 1:N-1
      x += dx
      res += f(x)
    end
    return res*dx

  elseif scheme == "simpson"
    # This is a way to check that N is even and throw an assertion error
    # if it is not
    @assert mod(N,2) == 0 "N must be even for Simpson's rule"

    res = f(a) + f(b)
    x = a
    for i = 1:N/2-1
      x += dx
      res += 4*f(x)
      x += dx
      res += 2*f(x)
    end
    x += dx
    res += 4*f(x)

    return res*dx/3

  else
    error("No scheme by this name")
  end

end

"""
    trisolve(a,b,c,f::Vector{<:Real},sys_type::String) -> Vector{Real}

Solve a tridiagonal system, in which the matrix consists of elements `a`, `b`, and `c`,
on the subdiagonal, diagonal, and superdiagonal, respectively, and the right-hand
side of the system is given by vector `f`. The `sys_type` of system can be either `"circulant"`
or `"regular"`. The solution is returned as a vector the same size as `f`.

For `"circulant"` systems, `a`, `b` and `c` can only be scalars, since the diagonals
must each be uniform. For `"regular"` systems, these can each be vectors, but `b` must
be the same length as `f` and `a` and `c` should be one element shorter.
"""
function trisolve(a,b,c,f::Vector{<:Real},sys_type::String)
  # Note that `f` is set up to only accept vectors whose elements are of type Real
  # or some subtype of Real (e.g. Float64). That is what the <: does.

  # Return immediately with zeros if f is all zeros
  if all(f .== 0)
    return f
  end

  if sys_type == "circulant"

    # For circulant matrices, use this special form
    return circulant(a,b,c,f)

  elseif sys_type == "regular"

    # For regular matrices, use this special form
    return thomas(a,b,c,f)

  end

end

# ======  Circulant tridiagonal matrices ==============#

function circulant(a::Real,b::Real,c::Real,f::Vector{<:Real})
  # For circulant matrices with only allow a, b, c to be scalars

  M = length(f)

  # create a solution vector that is just a copy of f. We will operate on this
  # vector "in place". We are also converting it to complex type for use in the
  # fft
  x = ComplexF64.(f)

  m = 0:M-1
  # Get the eigenvalues of the matrix
  # In Julia, you can get special symbols like π by typing backslash, then pi,
  # then [TAB].
  lam = b .+ (a+c)*cos.(2π*m/M) .- im*(a-c)*sin.(2π*m/M)

  # x now temporarily holds the transform of the right-hand side
  fft!(x)

  # Scale x by M to make it agree with what we usually call the discrete
  # Fourier transform
  x ./= M

  # now divide the transform of f by the eigenvalues of the matrix
  # This is a shorthand take each element of x, divide it by the corresponding
  # entry in lam, and put the result in x.
  x ./= lam

  # Now inverse transform the result and scale appropriately and return just
  # the real part of each element
  ifft!(x)

  x .*= M # multiply by the vector length to get the "correct" form of IFFT

  return real.(x)

end

# ======  Regular tridiagonal matrices ==============#

# For scalar-valued a, b, c, call the vector-valued form. This ensures
# that we have only one function to write and debug
function thomas(a::Real,b::Real,c::Real,f::Vector{<:Real})
  n = length(f)
  ftype = eltype(f) # this is to ensure that a, b, c elements have same type as f
  return thomas(a*ones(ftype,n-1),
                b*ones(ftype,n),
                c*ones(ftype,n-1),f)
end

function thomas(a::Vector{<:Real},b::Vector{<:Real},c::Vector{<:Real},f::Vector{<:Real})
  # This is the main algorithm for regular circulant matrices

  # We will develop the solution x in place
  x = copy(f)

  M = length(f)

  if b[1] == 0
    error("First pivot is zero")
  end

  # Modify the first row of coefficients
  x[1] /= b[1]

  if M == 1
    # There is only one entry, so we are done!
    return x
  end

  c[1] /= b[1]
  tmp1 = c[1]
  tmp2 = x[1]
  for i = 2:M-1
    tmp = b[i] - a[i-1]*tmp1
    tmp1 = c[i]/tmp
    tmp2 = (x[i]-a[i-1]*tmp2)/tmp
    c[i] = tmp1
    x[i] = tmp2
  end
  x[M] = (x[M]-a[M-1]*x[M-1])/(b[M]-a[M-1]*c[M-1])

  # Back substitute
  tmp = x[M]
  for i = M-1:-1:1
    tmp = x[i] - c[i]*tmp
    x[i] = tmp
  end

  return x

end


end # module
