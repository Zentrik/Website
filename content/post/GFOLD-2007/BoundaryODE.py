import pykep

g = [0, 0, -3.7114];
cosphi = cos(deg2rad(27));
α = 1 /(225 * 9.807 * cosphi);

T_max = 3100;
n = 6;
T1 = 0.3 * T_max;
T2 = 0.8 * T_max;
rho1 = n * T1 * cosphi;
rho2 = n * T2 * cosphi;

mwet = 1905;

r0 = [2000, 0, 1500] # x y z, unlike in paper BE CAREFUL
rdot0 = [100, 0, -75]

rfinal = [0, 0, 0]
rdotfinal = [0, 0, 0]

x0 = [r0, rdot0, mwet]

l = pykep.pontryagin.leg(pykep.epoch(), x0, l0, tf, xf)