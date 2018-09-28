global Iext
global Phi
global Gca
global Gk
global Gl
global Vca
global Vk
global Vl
global V1
global V2
global V3
global V4
global C

Iext = 0;
Phi = 0.02;
Gca = 4.4;
Gk = 8.0;
Gl = 2;
Vca = 120;
Vk = -84;
Vl = -60;
V1 = -1.2;
V2 = 18;
V3 = 2;
V4 = 30;
C = 20;

%%% Question 2

figure(1)
% Plotting nullclines
nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));

% Quiver
[V,w] = meshgrid(-100:5:100,0:0.1:1);
dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
q = quiver(V,100*w,dV,100*dw); 

% Finding equilibrium points
y_init = [-60; 0];
equiPoint = fsolve(@mle,y_init);
plot(equiPoint(1),100*equiPoint(2),'r*')
disp(equiPoint)

%%% Question 3

% Jacobi and Eigen values at equilibrium points 
syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint'));
disp(A);
[V,D] = eig(A)


deltaT = [0, 500];
equiPoint1 = equiPoint

%%% Question 5

% Action potential 

Iext = 0;
equiPoint1(1) = equiPoint(1) + 1000/C;

Phi= 0.01;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

Phi = 0.02;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

Phi = 0.04;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

%%% Question 6

Iext = 0;
Phi = 0.02;

equiPoint1(1) = equiPoint(1) + 300/C;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1(1) = equiPoint(1) + 500/C;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1(1) = equiPoint(1) + 2000/C;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1(1) = equiPoint(1) + 3000/C;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))


legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point','Phi=0.01','Phi=0.02','Phi=0.04','Excitation 1','Excitation 2','Excitation 3','Excitation 4');

%%% Question 7

figure(2)

Iext = 86;

nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));

% Quiver
[V,w] = meshgrid(-100:5:100,0:0.1:1);
dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
q = quiver(V,100*w,dV,100*dw); 

% Finding equilibrium points
y_init = [-60; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint2)

% Shows limit cycle
equiPoint1(1) = equiPoint(1);
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1 = equiPoint2;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

% Inward spiral
[t, y] = ode45(@mle4ode,deltaT,[-27.9;0.17]);
plot(y(:,1),100*y(:,2))

%%% Question 8

equiPoint1 = equiPoint2 ;
[t, y] = ode45(@mle4odeRev,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

%%% Question 9


function dYdt = mle4ode(t,y)
    global Iext
    global Phi
    global Gca
    global Gk
    global Gl
    global Vca
    global Vk
    global Vl
    global V1
    global V2
    global V3
    global V4
    global C
    V = y(1);
    w = y(2);
    dYdt1 = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
    dYdt2 = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
    dYdt = [dYdt1;dYdt2]; 
end

function dYdt = mle4odeRev(t,y)
    global Iext
    global Phi
    global Gca
    global Gk
    global Gl
    global Vca
    global Vk
    global Vl
    global V1
    global V2
    global V3
    global V4
    global C
    V = y(1);
    w = y(2);
    dYdt1 = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
    dYdt2 = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
    dYdt = [-1*dYdt1;-1*dYdt2]; 
end


function dYdt = mle(y)
    global Iext
    global Phi
    global Gca
    global Gk
    global Gl
    global Vca
    global Vk
    global Vl
    global V1
    global V2
    global V3
    global V4
    global C
    V = y(1);
    w = y(2);
    dYdt1 = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
    dYdt2 = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
    dYdt = [dYdt1;dYdt2]; 
end
