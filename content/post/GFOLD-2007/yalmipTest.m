g0    = 9.80665;         % standard gravity [m/s**2]
m_dry = 1505;            % dry mass kg
m_fuel= 400;             % fuel in kg
T_max = 3100 * 6 * cos(deg2rad(27));           % thrust max
throt = [0.3, 0.8];       % throttle ability
Isp   = 225;             % fuel efficiency (specific impulse)
V_max = 90;              % velocity max
y_gs  = deg2rad(4); % glide slope cone, must be 1 < Degrees < 91
p_cs  = deg2rad(45);  % thrust pointing constraint
alpha = 1/(Isp*g0);      % fuel consumption parameter
m_wet = (m_dry+m_fuel);  % wet mass kg
r1    = throt(1)*T_max;  % lower thrust bound
r2    = throt(2)*T_max;  % upper thrust bound

g = [-3.7114; 0; 0];                   % gravity
nh= [1; 0; 0];                     % thrust vector reference direction

r0 = [1500; 0; 2000];              % initial position
v0 = [-75;  0;   100];                 % initial velocity
x0 = [r0; v0];

rf = [0; 0; 0];                      % final position target
vf = [0; 0; 0];                       % final velocity target

c = [1; 0; 0] / tan(y_gs);


tf = 81;
dt = 0.5;
N = floor(tf / dt + 1);

% It's good practice to start by clearing YALMIPs internal database 
% Every time you call sdpvar etc, an internal database grows larger
yalmip('clear')

% Define variables
x = sdpvar(6, N); % state vector (3position,3velocity)
u = sdpvar(3, N); % u = Tc/mass because Tc[:,n)/m[n) is not allowed by DCP
z = sdpvar(1, N);  % z = ln(mass)
s = sdpvar(1, N); % thrust slack parameter

% Define constraints 
Constraints = [x(:, 1) == x0];
Constraints = [Constraints, x(4:6, N) == vf]; % don't forget to slow down, buddy!

Constraints = [Constraints, s(N) == 0]; % thrust at the end must be zero. WE NEED BOUND ON FINAL THRUST.
%u(:, 1) == s(1) * [1; 0; 0]; % thrust direction starts straight
%u(:, N) == s(N-1) * [1; 0; 0]; % and ends straight
Constraints = [Constraints, z(1) == log(m_wet)]; % convexified (7)
%Constraints = [Constraints, log(m_dry) <= z(2:end) <= log(m_wet)];
Constraints = [Constraints, z(N) >= log(m_dry)];

for n = 1:N-1 % any t in [0,tf] maps to any n in [1,N-1)

    % Leapfrog Integration Method
    %    accurate +/- sqrt( (dt*df/dr)**2 + 1)
    %    https://goo.gl/jssWkB
    %    https://en.wikipedia.org/wiki/Leapfrog_integration

    % Dynamics --> v = A(w)*x + B*(g + u)

    Constraints = [Constraints, x(4:6, n+1) == x(4:6, n) + (dt/2)*((u(:,n)+g) + (u(:,n+1)+g))];
    Constraints = [Constraints, x(1:3, n+1) == x(1:3, n) + (dt/2)*(x(4:6,n+1)+x(4:6,n)) + dt^2/12 * (u(:,n+1) - u(:,n))];

    %Constraints = [Constraints, norm( (x(2:3,n) - rf(2:3)) ) - c(1)*(x(1,n)-rf(1)) <= 0]; % specific, but faster

    Constraints = [Constraints, z(n+1) == z(n) - (alpha*dt/2)*(s(n) + s(n+1))]; % mass decreases
    Constraints = [Constraints, norm(u(:,n)) <= s(n)]; % limit thrust magnitude & also therefore, mass

    z0_term = m_wet - alpha * r2 * (n - 1) * dt;  % see ref [2], eq 34,35,36
    z1_term = m_wet - alpha * r1 * (n - 1) * dt;
    z0 = log( z0_term );
    z1 = log( z1_term );
    mu_1 = r1 * exp(-z0);
    mu_2 = r2 * exp(-z0);
    delta = (z(n) - z0);

    Constraints = [Constraints, mu_1 * (1 - delta + delta^2 / 2) <= s(n) <= mu_2 * (1 - delta)]; % thrust bound

    Constraints = [Constraints, z0 <= z(n) <= z1];

end
%Constraints = [Constraints, x(1, :) >= 0]; % no, this is not the Boring Company!

% Define an objective
Objective = -z(N);
Constraints = [Constraints, x(1:3, N) == rf];

% Objective = norm(x(2:3,N)-rf(2:3));
% Constraints = [Constraints, x(1, N) == 0]; % land on ground

% Set some options for YALMIP and solver
options = sdpsettings('savesolveroutput',1);

% Solve the problem
sol = optimize(Constraints, Objective, options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 xval = value(x);
 T = 0:dt:tf;
 plot(T, xval(1:3, :))
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end