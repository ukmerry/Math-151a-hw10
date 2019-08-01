tf% qSplineTest
%
% A program to test the accuracy of a quadratic
% spline approximation based upon data values at n+1 equispaced 
% knots oveçr the interval [xMin,xMax].
%              
%
% Math 151A, Winter 2018 (03/08/2018)
%


F     = @(x)sin(x);  % Target test function 
dF    = @(x)cos(x);

xMin = 0.0;
xMax = 2.0;

nStart  = 10;       % Starting number of panels 

nRefine = 5;        % Number of refinements. 0 refinements 
                    % leads to a single construction with 
                    % nStart panels. 

for kStep = 1:nRefine+1

n     = nStart*2^(kStep-1);  % number of panels
h     = (xMax-xMin)/n;       % panel size 


% Create equispaced data samples

f = zeros(n+1,1);

for i = 0:n
  x_i    = xMin + i*h;
  f(i+1) = F(x_i);
end


% Data allocation for spline coefficients using a spline 
% representation over the ith panal of 
%
% S_i(x) = a(i) + b(i)*(x - x_(i-1)) + c(i)*(x - x_(i-1))^2
%

a     = zeros(n,1);  
b     = zeros(n,1); 
c     = zeros(n,1); 

% 
% Data allocation for n+1 equations for the q(i) values where 
% q(i) = S'(x_(i-1)) i = 1 .. n + 1 
%

q     = zeros(3*n,1);     
A     = zeros(3*n,3*n);  % Matrix determining q values
qRHS  = zeros(3*n,1);  % Right hand side of equations for q


% Set up A and qRHS, the equations for the q(i) values  

g=1; 
for i=1:n
    A(i,g:g+2) = [1 h h^2];
    qRHS(i)    = f(i);
    g          = g+3;
end
p=1;j=1;
for i=n+1:2*n
    A(i,p:p+2) = [1 2*h 4*h^2];
    qRHS(i)    = f(j+1);
    p          = p+3;
    j          = j+1;
end
e=2;
for i=(2*n+1):(3*n-1)
    A(i,e:e+3) = [1 2*h 0 -1];
    qRHS(i)    = 0;
    e = e+3;
end
A(3*n,2)       = 1;
qRHS(3*n)      = dF(xMin)
% Solve the system 

q = A\qRHS;

% Determine the values of the coefficients
% a(i), b(i) and c(i) for i = 1 .. n from the q(i) values.
%

j=1;
for i=1:n
    a(i)=q(j);
    j = j+3;
end
j=2;
for i=1:n
    b(i)=q(j)
    j = j+3;
end
j=3;
for i=3:3:3*n
    c(i)=q(j);
    j = j+3;
end

% Evaluate the spline at 200 equispaced points 

nSample = 200;
x       = xMin:(xMax-xMin)/nSample:xMax;

s       = qSplineEval(x,a,b,c,xMin,xMax); % spline evaluation

exactF = F(x);                            % exact  evaluation 


% Evaluate the max of the error over sampling points 
% and store the error in the vector splineErrorNorm 

splineError            = s - exactF;
splineErrorNorm(kStep) = norm(splineError,inf);

fprintf('Panels : %-3d  Error Max:  %-15.10e \n',n,splineErrorNorm(kStep));


plot(x,exactF,x,s);
axis([xMin,xMax,-0.5,1.5]);
legend('y = sin(x)','y = spline approximation');

% Pause the computation after each refinement step so that 
% one can view the plots 

fprintf ('press any key to advance to next refinement \n');
pause



end % End of loop over refinements 

%
% End of computation loop 
%

% Process the error 

Rate = zeros(nRefine+1,1);
for i=3:(nRefine+1)
    Rate(i) = log(splineErrorNorm(i-1)/splineErrorNorm(i))/log(splineErrorNorm(i-2)/splineErrorNorm(i-1));
end

fs = [' # panels            Error                        Convergence',sprintf('\n')];

for kStep = 1:nRefine+1
  n  = nStart*2^(kStep-1);
  fs = [fs,sprintf('  %-3d              %-15.10e              %-15.10e \n',n,splineErrorNorm(kStep), Rate(kStep))];
end
fs