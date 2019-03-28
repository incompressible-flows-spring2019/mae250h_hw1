using Test
using HW1


fcn(x) = 2*x

@test HW1.numerical_integration(fcn,0,1) ≈ 1.0 atol = 1e-14

fcn(x) = x^3

@test HW1.numerical_integration(fcn,-1,1) ≈ 0.0 atol = 1e-14
