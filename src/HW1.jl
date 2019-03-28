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

  M = length(f)


  # Return immediately with zeros if f is all zeros
  if all(f .== 0)
    return f
  end


  if sys_type == "circulant"
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

    x .*= M

    return real.(x)


  elseif sys_type == "regular"
    println("That is your job!")
    return nothing
  end

end



end # module
