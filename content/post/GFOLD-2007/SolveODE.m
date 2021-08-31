syms Gamma lambda0(t) lambda3(t) m(t) H r1(t) r2(t) r3(t) rdot1(t) rdot2(t) rdot3(t) l11(t) l12(t) l13(t) l21(t) l22(t) l23(t) tfinal
r(t) = [r1(t), r2(t), r3(t)];
rdot(t) = [rdot1(t), rdot2(t), rdot3(t)];
lambda1(t) = [l11(t), l12(t), l13(t)];
lambda2(t) = [l21(t), l22(t), l23(t)];

g = [0 0 -3.7114];
cosphi = cos(deg2rad(27));
alpha = 1 /(225 * 9.807 * cosphi);

T_max = 3100;
n = 6;
T1 = 0.3 * T_max;
T2 = 0.8 * T_max;
rho1 = n * T1 * cosphi;
rho2 = n * T2 * cosphi;

mwet = 1905;

R12 = lambda0 - alpha * lambda3 + norm(lambda2 / m);
Gamma = rho2 ;%rho1 + heaviside(R12) * (rho2 - rho1);
T_c = lambda2 / norm(lambda2) * Gamma;

eq1 = diff(r, t) == rdot;
eq2 = diff(rdot, t) == g + T_c / m;
eq3 = diff(m, t) == - alpha * Gamma;

cost = Gamma;

H = [lambda0 lambda1 lambda2 lambda3] * [cost, rdot, g + T_c / m, - alpha * Gamma].';

eq4 = diff(lambda1, t) == - [0 0 0];
eq5 = diff(lambda2, t) == - lambda1;
eq6 = diff(lambda3, t) == - norm(lambda2) * Gamma / m^2;

%sol = dsolve([eq1, eq2, eq3, eq4, eq5, eq6], [r(0) == [1500 0 2000], rdot(0) == [-75 0 100], m(0) == mwet, r(tfinal) == [0 0 0], rdot(tfinal) == [0 0 0], lambda3(tfinal) == 0])

r0 = [2000, 0, 1500]; % x y z, unlike in paper BE CAREFUL
rdot0 = [100, 0, -75];

rfinal = [0, 0, 0];
rdotfinal = [0, 0, 0];

lambda0 = -1;
lambda_guess = randn(1, 7);

t0 = 0;
Deltat_guess = 72;

%jac(0, [r0, rdot0, mwet, lambda_guess], Deltat_guess)

solinit = bvpinit(linspace(0,1),[r0, rdot0, mwet, lambda_guess], Deltat_guess);
opts = bvpset('FJacobian',@(s, u, Deltat) jac(s, u, Deltat, j_fun, lambda0), 'Stats','on');
sol = bvp4c(@(s, u, Deltat)ode(s,u,Deltat, lambda0, g, alpha, rho1, rho2, t0), @(ya, yb, Deltat)bc(ya,yb,Deltat, r0, rdot0, rfinal, rdotfinal, mwet, g, alpha, rho1, rho2, lambda0), solinit, opts);
y = sol.y;
time = sol.parameters*sol.x;
ut = -y(4,:);

% -----------------------------------------------------------
% ODE’s of augmented states
function du = ode(s,u,Deltat, lambda0, g, alpha, rho1, rho2, t0)
  r = u(1:3);
  rdot = u(4:6);
  m = u(7);
  lambda1 = u(8:10);
  lambda2 = u(11:13);
  lambda3 = u(14);

  real_t = t0 + Deltat * s;

  R12 = lambda0 - alpha * lambda3 + norm(lambda2 / m);
  if (R12 > 0) 
    Gamma = rho2;
  else
    Gamma = rho1;
  end
  
  T_c = lambda2 / norm(lambda2) * Gamma;

  du = zeros(14, 1);
  du(1:3) = rdot;
  du(4:6) = g.' + T_c / m;
  du(7) = -alpha * Gamma;
  du(8:10) = [0; 0; 0];
  du(11:13) = -lambda1;
  du(14) = dot(lambda2, T_c) / m^2;

  du = Deltat * du;
end
% -----------------------------------------------------------
% boundary conditions: x1(0)=1;x2(0)=2, x1(tf)=3, p2(tf)=0;
% p1(tf)*x2(tf)-0.5*p2(2)^2
function residual = bc(ya,yb,Deltat, r0, rdot0, rfinal, rdotfinal, mwet, g, alpha, rho1, rho2, lambda0)
  residual = zeros(15, 1);
  residual(1:3) = ya(1:3) - r0.';
  residual(4:6) = ya(4:6) - rdot0.';
  residual(7) = ya(7) - mwet;
  residual(8:10) = yb(1:3) - rfinal.';
  residual(11:13) = yb(4:6) - rdotfinal.';
  residual(14) = yb(14) - 0;
  if (Deltat > 0) 
    residual(15) = 0; 
  else
    residual(15) = Deltat;
  end
  residual(15) = 0; 
%   rdot = yb(4:6);
%   m = yb(7);
%   lambda1 = yb(8:10);
%   lambda2 = yb(11:13);
%   lambda3 = yb(14);
%   
%   R12 = lambda0 - alpha * lambda3 + norm(lambda2 / m);
%   if (R12 > 0) 
%     Gamma = rho2;
%   else
%     Gamma = rho1;
%   end
%   T_c = lambda2 / norm(lambda2) * Gamma;
%   
%   H = dot([lambda0; lambda1; lambda2; lambda3], [Gamma; rdot; g.' + T_c / m; - alpha * Gamma]);
%   residual(16) = H - 0;
end


function [dfdy, dfdpar] = jac(s, u, Deltat, j_fun, lambda0)
  r = u(1:3);
  rdot = u(4:6);
  m = u(7);
  lambda1 = u(8:10);
  lambda2 = u(11:13);
  lambda3 = u(14);
  
  dfdy = Deltat * j_fun(lambda0,lambda3,lambda2(1),lambda2(2),lambda2(3), m);
  dfdpar = ones(14, 1);
end