global Iext
global Phi
global Gca
global Gk Gna
global Gl
global Vca Vna
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
[V,D] = eig(A)


deltaT = [0, 500];
equiPoint1 = equiPoint;
legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point')
%%% Question 5
figure(2)

V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));

[V,w] = meshgrid(-100:5:100,0:0.1:1);
dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
q = quiver(V,100*w,dV,100*dw); 

plot(equiPoint(1),100*equiPoint(2),'r*')

% Action potential 

Iext = 0;
equiPoint1(1) = equiPoint(1) + 1*1000/C;

Phi= 0.01;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

Phi = 0.02;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

Phi = 0.04;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point','Phi=0.01','Phi=0.02','Phi=0.04');

%%% Question 6
figure(3)


V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));

Iext = 0;
Phi = 0.02;

[V,w] = meshgrid(-100:5:100,0:0.1:1);
dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
q = quiver(V,100*w,dV,100*dw); 

plot(equiPoint(1),100*equiPoint(2),'r*')

legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point');


offsets = 0:1:200;
Peaks = zeros(size(offsets));
for i = 1:length(offsets) 
    equiPoint1(1) = equiPoint(1) + offsets(i);
    [t, y] = ode45(@mle4ode,deltaT,equiPoint1);
    plot(y(:,1),100*y(:,2));
    Peaks(i) = max(y(:,1));
end
figure(4)
plot(equiPoint(1)+offsets,Peaks);


legend('V vs Peak Amplitude of Action Potential');

%%% Question 7

figure(5)

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

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)


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

legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point','From equilibrium point of Iext = 0','On equilibrium Point','Start point [-27.9;0.17]');

%%% Question 8

equiPoint1(1) = equiPoint2(1) +1 ;
[t, y] = ode45(@mle4odeRev,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

%%% Question 9

figure(6)

%%% Iext = 80

for i = 80:1:100

    Iext = i;

    nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
    nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
    V = -100:1:100;
    plot(V,100*nc1(V)); hold on;

    ylim([0 100])
    xlim([-100 100])

    plot(V,100*nc2(V));
    
    y_init = [-60; 0];
    equiPoint2 = fsolve(@mle,y_init);
    plot(equiPoint2(1),100*equiPoint2(2),'r*')
    disp(equiPoint2)

end
% Quiver
% [V,w] = meshgrid(-100:5:100,0:0.1:1);
% dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
% dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
% q = quiver(V,100*w,dV,100*dw); 

% Finding equilibrium points
y_init = [-60; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint2)

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)


%%% Iext = 86

Iext = 86;

nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));


% Quiver
% [V,w] = meshgrid(-100:5:100,0:0.1:1);
% dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
% dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
% q = quiver(V,100*w,dV,100*dw); 

% Finding equilibrium points
y_init = [-60; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint2)

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)



%%% Iext = 90


Iext = 90;

nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));


% Quiver
% [V,w] = meshgrid(-100:5:100,0:0.1:1);
% dV = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
% dw = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
% q = quiver(V,100*w,dV,100*dw); 

% Finding equilibrium points
y_init = [-60; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint2)

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)


%%% Question 10 



Iext = 30;
Phi = 0.0667;
Gca = 4;
Gk = 8.0;
Gl = 2;
Vca = 120;
Vk = -84;
Vl = -60;
V1 = -1.2;
V2 = 18;
V3 = 12;
V4 = 17.4;
C = 20;


figure(7)
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

y_init = [-20; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint)

y_init = [10; 0];
equiPoint3 = fsolve(@mle,y_init);
plot(equiPoint3(1),100*equiPoint3(2),'r*')
disp(equiPoint)

% Jacobi and Eigen values at equilibrium points 
syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint'));
[V,D] = eig(A)


deltaT = [0, 2000];
equiPoint1 = equiPoint;
equiPoint1(1) = equiPoint(1)+10;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1 = equiPoint2;
equiPoint1(1) =  equiPoint2(1) + 1;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1 = equiPoint3;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))


%%% Question 11

%%% Question 12
% figure(8)

Iext = 0;
Gk = 36;
Gl = 0.3;
Gna = 120;
Vna = 55;
Vk = -72;
Vl = 10;
C = 1;

% nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
% nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
% V = -100:1:100;
% plot(V,100*nc1(V)); hold on;

% ylim([0 100])
% xlim([-100 100])

% plot(V,100*nc2(V));

% Quiver
% [V,m,n] = meshgrid(-100:5:100,0:0.1:1,0:0.1:1);
% alphan = -0.01*(V+50)./(exp(-(V+50)/10)-1+10^-12);
% betan = 0.125*exp(-(V+60)/80);
% alpham = -0.1*(V+35)./(exp(-(V+35)/10)-1 + 10^-12);
% betam = 4*exp(-(V+60)/18);
% dVdt = (Iext -Gna.*m.^3.*(V-Vna) - Gk*n.^4.*(V-Vk) - Gl*(V-Vl))/C;
% dmdt = alpham.*(1-m) - betam.*m;
% dndt = alphan.*(1-n) - betan.*n;
% q = quiver3(V,100*m,100*n,dVdt,100*dmdt,100*dndt); 


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


function dYdt = hh4ode(t,y)
    global Iext
    global Gk
    global Gl Gna
    global Vna
    global Vk
    global Vl
    global C
    V = y(1);
    n = y(2);
    m = y(3);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12);
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12);
    betam = 4*exp(-(V+60)/18);
    dVdt = (Iext -Gna*m^3*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dYdt = [dVdt;dndt;dmdt]; 
end

