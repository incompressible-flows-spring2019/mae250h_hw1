module HW1

"""
    numerical_integration(fcn,a,b;N=100) -> Number

Numerically integrate a function `fcn` using trapezoidal quadrature. The
function `fcn` should take a single argument \$x\$. The function is integrated from
\$x = a\$ to \$x = b\$, and return the result. The optional argument `N` sets the
number of quadrature steps, and defaults to 100.
"""
function numerical_integration(fcn,a::Real,b::Real;N=100)

  dx = (b-a)/N
  x = a
  result = 0.5*(fcn(a)+fcn(b))
  for i = 1:N-1
    x += dx
    result += fcn(x)
  end

  return result*dx
end



end # module
