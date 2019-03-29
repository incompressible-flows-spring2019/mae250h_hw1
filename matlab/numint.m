function [ integ ] = numint( func, llim, ulim, intops, params )
%NUMINT Use a specified numerical integration scheme to evaluate integral
%   integ = NUMINT(func,llim,ulim,intops,params) numerically integrates a
%      1-d function (specified by the function pointer 'func') between llim
%      and ulim. The integrand function must have the form
%          y = func(x,params)
%      where params is a structure of parameters used by the function.
%      Various numerical integration parameters are specified
%      by the structure 'intops':
%          intops.N = number of integration steps
%          intops.scheme = integration scheme
%                          1 = Trapezoidal
%                          2 = Simpson's rule
%                          3 = Adaptive Simpson's rule, using Matlab quad
%
%   MAE 259A
%   J. D. Eldredge
%   1/3/2014


%% Unpack the integration parameters
scheme = intops.scheme;
N = intops.N;

%% Perform the integration
integ = 0;
dx = (ulim-llim)/N;
x = llim;
switch scheme
    case 1,
        % Trapezoidal integration
        integ = 0.5*(func(llim,params) + func(ulim,params));
        for ii = 1:N-1
            x = x + dx;
            integ = integ + func(x,params);
        end
        integ = integ*dx;
    case 2,
        % Simpson's rule with trapezoidal end
        if (mod(N,2)~=0),
            disp('ERROR. N must be even for Simpson''s rule.');
            return
        end
        integ = func(llim,params) + func(ulim,params);
        for ii = 1:N/2-1
            x = x + dx;
            integ = integ + 4*func(x,params);
            x = x + dx;
            integ = integ + 2*func(x,params);
        end
        x = x + dx;
        integ = integ + 4*func(x,params);
        integ = integ*dx/3;
    case 3
        % Adaptive Simpson, using built-in Matlab 'quad' function
        integ = quad(@(x) func(x,params),llim,ulim,intops.tol,1);
    otherwise
        disp('Invalid choice of integration scheme.');
        return
end


end

