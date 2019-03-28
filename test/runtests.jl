# A package of useful testing commands
using Test

# We need to tell this test to load the module we are testing
using HW1

# The Random package is useful for making tests as arbitrary as possible
using Random


#=
This test checks that scale_it will fail if we try to pass something other than
a function into the first argument. There should be no form of scale_it that
works for such a case. We check that it "throws" an error
=#
@test_throws MethodError HW1.scale_it(1,1,1)

# This test checks that the desired result is obtained
fcn(x) = x^2
@test HW1.scale_it(fcn,1,2) == 2

#=
This check is similar to the previous one, but floating-point tests can
be challenging because floating-point 0 is never quite zero on a computer.
Here, we use an approximate check, to within a tolerance that would be
reasonable for 64-bit floating-point numbers (Float64)
=#
fcn(x) = sin(x)
@test HW1.scale_it(fcn,pi,2) ≈ 0 atol=1e-14


# Setting up a separate test set is useful to compartmentalize the testing
# of different parts of the code
@testset "trisolve" begin

  # Our routine should not accept arguments that are inappropriate:
  @test_throws MethodError HW1.trisolve(1,1,1,1,1)

  #=
  Our routine should return a vector of 0 if the rhs is 0, regardless
  of the coefficients in the matrix. (Though it would be useful to check
  for singular matrices
  =#
  x = HW1.trisolve(rand(),rand(),rand(),zeros(10),"circulant")
  @test all(x.== 0)

  #=
  Now a test of a circulant matrix with non-trivial result. We will test a
  4x4 system with random choices for coefficients and right-hand side vector
  and compare our answer with A\rhs
  =#
  a = rand(); b = rand(); c = rand()

  A = [b c 0 a; a b c 0; 0 a b c; c 0 a b]
  rhs = rand(4)

  x = HW1.trisolve(a,b,c,rhs,"circulant") # our solution
  x2 = A\rhs  # Julia's solution with native solver

  @test x ≈ x2

end
