global Iext
global Phi

Iext = 0;
Phi = 0.02;

% Plotting nullclines
nc1 = @(V)(-2.2*(1+tanh((V+1.2)/18)).*(V-120) - 2*(V+60))./(8*(V+84));
nc2 = @(V)0.5*(1+tanh((V-2)/30));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;
plot(V,100*nc2(V));

% Quiver
[V,m] = meshgrid(-100:5:100,0:0.1:1);
dV = -0.11*(1+tanh((V+1.2)/18)).*(V-120) - 0.4*m.*(V+84) - 0.1*(V+60);
dm = (0.5*(1+tanh((V-2)/30)) - m)*0.02.*cosh((V-2)/60);
q = quiver(V,100*m,dV,100*dm); 
% Finding equilibrium points
y_init = [-60; 0];
equiPoint = fsolve(@mle,y_init);
plot(equiPoint(1),100*equiPoint(2),'r*')
disp(equiPoint)

% Jacobi and Eigen values at equilibrium points 
syms v m
J = jacobian([-0.11*(1+tanh((v+1.2)/18))*(v-120) - 0.4*m.*(v+84) - 0.1*(v+60), (0.5*(1+tanh((v-2)/30)) - m)*0.02.*cosh((v-2)/60)],[v,m]);
A = double(subs(J,[v,m],equiPoint'));
disp(A);
[V,D] = eig(A)

deltaT = [0, 500];
equiPoint1 = equiPoint

Iext = 4
% Generating action potential
[t, y] = ode45(@mle4ode,deltaT,equiPoint);
plot(y(:,1),100*y(:,2))

Iext = 0
% Depolarizing current pulses
equiPoint1(1) = equiPoint(1) + 50;

Phi= 0.01
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
%plot(y(:,1),100*y(:,2))

Phi = 0.02
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
%plot(y(:,1),100*y(:,2))

Phi = 0.04
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
%plot(y(:,1),100*y(:,2))

legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point')

ylim([0 100])
xlim([-100 100])

function dYdt = mle4ode(t,y)
    global Iext
    global Phi
    dYdt1 = (Iext -0.11 * (1+tanh((y(1)+1.2)/18))*(y(1)-120) - 0.4*y(2)*(y(1)+84) - 0.1*(y(1)+60));
    dYdt2 = ((0.5*(1+tanh((y(1)-2)/30)) - y(2))*Phi*cosh((y(1)-2)/60));
    dYdt = [dYdt1;dYdt2]; 
end

function dYdt = mle(y)
    dYdt1 = -0.11*(1+tanh((y(1)+1.2)/18))*(y(1)-120) - 0.4*y(2)*(y(1)+84) - 0.1*(y(1)+60);
    dYdt2 = (0.5*(1+tanh((y(1)-2)/30)) - y(2))*0.02*cosh((y(1)-2)/60);
    dYdt = [dYdt1;dYdt2]; 
end
