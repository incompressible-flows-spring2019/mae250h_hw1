function [ x ] = trisolve(a, b, c, f, type )
%TRISOLVE Solve a simple tridiagonal system of equations
%   Given a matrix A with a tridiagonal structure specified by the
%   elements a, b, and c, and the right-hand side column array f, solve the
%   system to find x. The size of the system is determined from the size of
%   f.
%   
%     x = trisolve(a,b,c,f,'reg') would set A as a regular tridiagonal
%       matrix. Here, a, b and c could be scalars (in which it is assumed
%       that each diagonal is constant) or vectors (in which b is the main
%       diagonal and should have the same length as f, and a and c are the
%       sub- and super-diagonal and have length one smaller than b).
%
%     x = trisolve(a,b,c,f,'circ') would set A as a circulant tridiagonal
%       matrix. Here, a, b and c are restricted to be scalars.
%
%   MAE 259A
%   J. D. Eldredge
%   1/5/2014

M = length(f);

if (strcmp(type,'circ')),
    % For a circulant matrix, get the eigenvalues
    m = 0:M-1;
    lam = b + (a+c)*cos(2*pi*m/M) - 1i*(a-c)*sin(2*pi*m/M);
    
    % Transform the rhs (fhat = V^-1*f).
    fhat = fft(f)/M;
    
    % Get xhat = Lam^-1*fhat. Note the use of .' rather than ' to compute
    % the non-conjugate transpose of the eigenvalue array.
    xhat = fhat./lam.';
    
    % Get x = V*xhat by transforming back.
    x = real(M*ifft(xhat));
elseif (strcmp(type,'reg'))
    
    % You do this part.
 
    % Should have separate treatments for scalar or vector a,b,c
        

   
    
    
end


end

