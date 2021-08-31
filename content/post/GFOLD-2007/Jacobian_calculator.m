syms lambda3 m lambda0
r = sym('r', [1 3]);
rdot = sym('rdot', [1 3]);
lambda1 = sym('lambda1', [1 3]);
lambda2 = sym('lambda2', [1 3]);

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

normlambda2 = ((lambda2(1))^2 + (lambda2(2))^2 + (lambda2(3))^2)^(1/2);
R12 = lambda0 - alpha * lambda3 + normlambda2 / m;
Gamma = rho1 + heaviside(R12) * (rho2 - rho1);
T_c = lambda2 / normlambda2 * Gamma;

f = [rdot, g + T_c / m, - alpha * Gamma, - [0 0 0], - lambda1, - norm(lambda2) * Gamma / m^2];
j = jacobian(f, [r, rdot, m, lambda1, lambda2, lambda3]);
j_fun = matlabFunction(j);