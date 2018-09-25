function [Tout,Xout,FuncEval]=ExplicitRungeKutta(fun,tspan,x0,h,s,AT,b,c,varargin)
% Parameters related to constant step size
hAT = h*AT;
hb = h*b;
hc = h*c;
%hd = h*d;

% Size parameters
x = x0;
t = tspan(1); % Initial time
tf = tspan(end); % Final time
N = length(tspan) - 1; %round((tf-t)/h); % Number of steps
nx = length(x0); % System size (dim(x))

% Allocate space
T = zeros(1,s); % Stage T
X = zeros(nx,s); % Stage X
F = zeros(nx,s); % Stage F

Tout = zeros(N+1,1); % Time for output
Xout = zeros(N+1,nx); % States for output

Xout(1,:) = x;
FuncEval = 0;
for n=1:N
    % Stage 1
    T(1) = t;
    X(:,1) = x;
    F(:,1) = feval(fun,T(1),X(:,1),varargin{:});
    % Stage 2,3,...,s
    T(2:s) = t + hc(2:s);
    for i=2:s
        X(:,i) = x + F(:,1:i-1)*hAT(1:i-1,i);
        F(:,i) = feval(fun,T(i),X(:,i),varargin{:});
        FuncEval = FuncEval +1;
    end
    % Next step
    t = tspan(n+1);
    x = x + F*hb;
    %e = F*hd;

    % Save output
    Tout(n+1) = t;
    Xout(n+1,:) = x';
end
end


