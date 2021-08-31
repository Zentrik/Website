function value = Yalmip(tf)
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

  mass = m_dry;
  Thrust = r1;

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

  Constraints = [Constraints, x(1, :) >= 0]; % no, this is not the Boring Company!

  Constraints = [Constraints, x(4:6, N) == vf]; % don't forget to slow down, buddy!

  for n = 1:N-1 % any t in [0,tf] maps to any n in [1,N-1)

      % Leapfrog Integration Method
      %    accurate +/- sqrt( (dt*df/dr)**2 + 1)
      %    https://goo.gl/jssWkB
      %    https://en.wikipedia.org/wiki/Leapfrog_integration

      % Dynamics --> v = A(w)*x + B*(g + u)

      Constraints = [Constraints, x(4:6, n+1) == x(4:6, n) + (dt/2)*((u(:,n)+g) + (u(:,n+1)+g))];
      Constraints = [Constraints, x(1:3, n+1) == x(1:3, n) + (dt/2)*(x(4:6,n+1)+x(4:6,n)) + dt^2/12 * (u(:,n+1) - u(:,n))];
  end
  for n = 1:N
      Constraints = [Constraints, norm( (x(2:3,n) - rf(2:3)) ) - c(1)*(x(1,n)-rf(1)) <= 0]; % specific, but faster

      Constraints = [Constraints, norm(u(:,n)) <= Thrust / mass]; % limit thrust magnitude & also therefore, mass
  end

  % Define an objective
  % Objective = -z(N);
  % Constraints = [Constraints, x(1:3, N) == rf];

  Objective = norm(x(2:3,N)-rf(2:3));
  Constraints = [Constraints, x(1, N) == 0]; % land on ground

  % Set some options for YALMIP and solver
  options = sdpsettings('savesolveroutput',1);

  % Solve the problem
  sol = optimize(Constraints, Objective, options);

  % Analyze error flags
  if sol.problem == 0
   % Extract and display value
   %sol = value(x);
   %T = 0:dt:tf;
   %plot(T, sol(1:3, :))
   
   value = sol;
  else
   display('Hmm, something went wrong!');
   sol.info
   yalmiperror(sol.problem)
   value = 'Fail';
  end
end