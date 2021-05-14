clear
close all

%% Set up the physical parameters
x=chebfun('x'); % variable x in [-1,1]
PDE.a=cosh(x); % coefficient a, the 0*x is important (so it's treated as a function)
PDE.b=sin(pi*x)+2; % coefficient b
RHO=tanh(x)+2;
PDE.rho_inv=1/(RHO); % represents 1/rho
PDE.nu=0.8; % the fractional order for time derivative
PDE.BC=3; % type of boundary condition: 1 is CC, 2 is SS, 3 is CS, 4 is SC (C=clamped, SS=supported)
PDE.FRAC_TYPE='C'; % RL is Riemann-Liouville, C is Caputo

u0=sin(2*pi*x)*(1+x)*(1-x)^2; % initial condition
u1=x*0; % initial time derivative
OMEGA=20; % hack to tell code whether to apply Cauchy to poles of not
PDE.F=@(t) cos(20*t)*sin(pi*x); % forcing
fhat=@(z) z./(400+z.^2)*sin(pi*x); % Laplace transform of the forcing

SING.POLES=[20*1i,-20*1i]; % poles of fhat
SING.RES={@(z) sin(pi*x)*(20*1i)/(20*1i+20*1i),@(z) sin(pi*x)*(-20*1i)/(-20*1i-20*1i)}; % WARNING - code currently only deals with simple poles
% SING.CAUCHY=1; % set to 1 if poles to right of contour, otherwise set to 0

NUM.t1=10; % final time we are interested in (must be finite!)
NUM.t0=1; % smallest time we are interested in (must be >0)

%% Set up the numerical parameters
NUM.DiscMin=10;         % min discretisation size for linear system solvers
NUM.DiscMax=1000;      % max discretisation size for linear system solvers
NUM.tol=10^(-16);       % tolerance of computation (in norm of vector of Chebyshev coeffs)
NUM.Parallel="on";      % set to "on" if you have access to multiple processors and wish to solve linear systems in parallel
NUM.Adaptive="on";      % set to "off" if you don't want to solve the linear systems adaptively
NUM.N=1000;             % use this to force size of discretisations when p.Adaptive="off"
NUM.Prog="off";          % set to "on" if you want to display progress of solving linear systems
NUM.CONT=1;             % type of contour used, 1 is hyperbolic, 2 is hyperbolic with optimised parameters, 3 is parabolic with optimal parameters

%% Compute suitable contour parameters - leave these lines unchanged!
% hyperbolic
NUM.beta=2; % max real part in complex exponential (determines the shape of hyperbola)
C=2*sqrt(max(PDE.a/PDE.b));
THETA1=pi*0.5000001:0.000001:min([pi/(2-PDE.nu),0.9999*pi,pi/(2*abs(1-PDE.nu))]);
RR1=(C*sqrt(-cos(THETA1)./(cos(THETA1*(PDE.nu-1))))./abs(sin(THETA1)-cos(THETA1).*(tan(THETA1*(PDE.nu-1))))).^(2/PDE.nu);
THETA2=-THETA1;
RR2=(C*sqrt(-cos(THETA2)./(cos(THETA2*(PDE.nu-1))))./abs(sin(THETA2)-cos(THETA2).*(tan(THETA2*(PDE.nu-1))))).^(2/PDE.nu);
NUM.shift=2/NUM.t1;
ANG=angle(RR1.*cos(THETA1)-NUM.shift+RR1.*sin(THETA1)*1i);
NUM.delta=pi-min(ANG);
% parabolic
if PDE.nu<1
    XL=log(10^(-16))/NUM.t0;
    XLL=abs(XL-RR1.*cos(THETA1));
    I=find(XLL==min(XLL));
    NUM.bb=abs(cos(THETA1(I(1)))/(RR1(I(1))*sin(THETA1(I(1))).^2));
else 
    NUM.bb=1/(C^2);
end

%% Perform the quadrature
NUM.Nquad=200; % number of quadrature points
[zi,vi]=frac_beam_quad(PDE,fhat,u0,u1,NUM,SING,OMEGA);

%% Compute and plot the solution - as set up, should agree with middle portio of figure 6 of paper
figure
plot(u0,'--k'); hold on
TVEC=[1,5,10]; % change this to any value in [NUM.t0,NUM.t1]
for j=1:length(TVEC)
    t=TVEC(j);
    U=chebfun(vi*exp(t*zi(:)),'coeffs');
    plot(real(U))
end



    
    