solinit = bvpinit(linspace(0,1),[2;3;1;1],2);
sol = bvp4c(@ode, @bc, solinit);
y = sol.y;
time = sol.parameters*sol.x;
ut = -y(4,:);

% -----------------------------------------------------------
% ODE’s of augmented states
function dydt = ode(t,y,T)
dydt = T*[ y(2);-y(4);0;-y(3) ];
end
% -----------------------------------------------------------
% boundary conditions: x1(0)=1;x2(0)=2, x1(tf)=3, p2(tf)=0;
% p1(tf)*x2(tf)-0.5*p2(2)^2
function res = bc(ya,yb,T)
res = [ ya(1) - 1; ya(2) - 2; yb(1) - 3; yb(4);
yb(3)*yb(2)-0.5*yb(4)^2];
end

%http://solmaz.eng.uci.edu/Teaching/MAE274/SolvingOptContProb_MATLAB.pdf