%%
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
global fni
global hinit

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


global equiP
options = odeset('Events',@myEvent);
optionsfsolve = optimset('Display','off');
%%% Question 2
%%
figure(1)
disp("Question 2")
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
xlabel('Membane Potential (in mV)') 
ylabel('% of K+ channels open') 
legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point')
title('Phase Plane plot')
%%
%%% Question 3
disp("Question 3")
% Jacobi and Eigen values at equilibrium points 
syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);

A = double(subs(J,[V,w],equiPoint'));
disp(A)
[V,D] = eig(A)


deltaT = [0, 500];
equiPoint1 = equiPoint;

%%
%%% Question 5
figure(2)
disp("Question 5")

y_init = [-60; 0];
equiPoint = fsolve(@mle,y_init);

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
deltaT = [0 500]
Iext = 0;
equiPoint1(1) = equiPoint(1) + 1*1000/C;
equiP = equiPoint;
Phi= 0.01;
[t, y,te,ye,ze] = ode45(@mle4ode,deltaT,equiPoint1,options);
plot(y(:,1),100*y(:,2))
display(te)
Phi = 0.02;
[t, y,te,ye,ze] = ode45(@mle4ode,deltaT,equiPoint1,options);
plot(y(:,1),100*y(:,2))
display(te)

Phi = 0.04;
[t, y,te,ye,ze] = ode45(@mle4ode,deltaT,equiPoint1,options);
plot(y(:,1),100*y(:,2))
display(te)
xlabel('Membane Potential (in mV)') 
ylabel('% of K+ channels open') 
title('Phase plane plot with varying Phi')
legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point','Phi=0.01','Phi=0.02','Phi=0.04');

%%
%%% Question 6
figure(3)
disp("Question 6")

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
xlabel('Membane Potential (in mV)') 
ylabel('% of K+ channels open') 
title('Phase plane plot with varying Phi')



offsets = 0:1:150;
Peaks = zeros(size(offsets));
thresould = -1000;
ev = double(equiPoint1(1));
for i = 1:length(offsets) 
    equiPoint1(1) = equiPoint(1) + offsets(i);
    [t, y] = ode45(@mle4ode,deltaT,equiPoint1);
    plot(y(:,1),100*y(:,2));
    Peaks(i) = max(y(:,1));
    if (thresould == -1000) && ((Iext -Gca*0.5*(1+tanh((ev-V1)/V2))*(ev-Vca) - Gk*w*(ev-Vk) - Gl*(ev-Vl))/C > 0)
        thresould = equiPoint1(1);
    end
end
display(thresould)
figure(4)
plot(equiPoint(1)+offsets,Peaks);

xlabel('starting voltage') 
ylabel('Peak voltage') 
title('V vs Peak Amplitude of Action Potential')

%%% Question 7
%%
disp("Question 7")
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

equiPoint1 = equiPoint;
plot(equiPoint(1),100*equiPoint(2),'g*')
% Shows limit cycle
equiPoint1(1) = equiPoint(1);
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1 = equiPoint2;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

% Inward spiral
plot(-27.9, 17,'b*')
[t, y] = ode45(@mle4ode,deltaT,[-27.9;0.17]);
plot(y(:,1),100*y(:,2))

xlabel('Membane Potential (in mV)') 
ylabel('% of K+ channels open') 
title('Phase plane plot')
legend('V Nullcline', 'm nullcline', 'Gradients' ,'Stabe Equilibrium point','Old Equilibrium point','Going to limit cycle','point inside UPO', 'stable Sprial');

%%
disp("Question 8")
figure(6)
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
J = jacobian([ -1*(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, -1*(0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)


%%% Question 8
deltaT = [0, 1000];
equiPoint1 = equiPoint2;
equiPoint1(1) = equiPoint2(1) +4 ;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

equiPoint1(1) = equiPoint2(1) +5;
[t, y] = ode45(@mle4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,2))

xlabel('Membane Potential (in mV)') 
ylabel('% of K+ channels open') 
title('Phase plane plot')
legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point','From inside limit cycle','From outside limit cycle');

%%
%%% Question 9

figure
disp("Question 9")

%%% Iext = 80

for i = 80:1:100

    Iext = i;
    
    equiP = equiPoint2;
    [t, y,te,ye,ze] = ode45(@mle4ode,deltaT,equiPoint,options);
    if size(ze) == 1
        disp("No spiking");
        disp(i);
    end        
    %plot(y(:,1),100*y(:,2))
    plot(t,y(:,1),'DisplayName',int2str(i)); hold on;
end
legend
figure

Iext = 80;

nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));

% Finding equilibrium points
y_init = [-60; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint2)

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)

Iext = 86;


nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));

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

% Finding equilibrium points
y_init = [-60; 0];
equiPoint2 = fsolve(@mle,y_init);
plot(equiPoint2(1),100*equiPoint2(2),'r*')
disp(equiPoint2)

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)

%%
%%% Question 10 

disp("Question 10")

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


figure(8)
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
plot(equiPoint2(1),100*equiPoint2(2),'g*')
disp(equiPoint2)

y_init = [10; 0];
equiPoint3 = fsolve(@mle,y_init);
plot(equiPoint3(1),100*equiPoint3(2),'b*')
disp(equiPoint3)

% Jacobi and Eigen values at equilibrium points 
syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);
A = double(subs(J,[V,w],equiPoint'));
[V,D] = eig(A)

syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);

A = double(subs(J,[V,w],equiPoint2'));
[V,D] = eig(A)
syms V w
J = jacobian([ (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Vca) - Gk*w*(V-Vk) - Gl*(V-Vl))/C, (0.5*(1+tanh((V-V3)/V4)) - w)*Phi*cosh((V-V3)/(V4*2))],[V,w]);

A = double(subs(J,[V,w],equiPoint3'));
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

xlabel('Membane Potential (in mV)') 
ylabel('% of K+ channels open') 
title('Phase plane plot')
legend('V Nullcline', 'm nullcline', 'Gradients' ,'Equilibrium point 1', 'Equilibrium point 2', 'Equilibrium point 3', 'deviation for equ1','deviation for equ2','deviation for equ3');



%%

%%% Question 11

figure
for i = 30:1:50

    Iext = i;
    
    nc1 = @(V)(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gl*(V-Vl))./(Gk*(V-Vk));
    nc2 = @(V)0.5*(1+tanh((V-V3)/V4));
    V = -100:1:100;
    plot(V,100*nc1(V)); hold on;

    ylim([0 100])
    xlim([-100 100])

    plot(V,100*nc2(V));
end

figure
for i = 30:1:50
    Iext = i;

    equiPoint1 = equiPoint2;
    equiPoint1(1) = equiPoint2(1) + 1;
    [t, y] = ode45(@mle4ode,deltaT,equiPoint1);
    %plot(y(:,1),100*y(:,2))
    plot(t,y(:,1),'DisplayName',int2str(i)); hold on;
end
legend


%%



%%% Question 12
Iext = 0;
Gk = 36;
Gl = 0.3;
Gna = 120;
Vna = 55;
Vk = -72;
C = 1;



%%
%%% Question 13

V = -60;
alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12);
betan = 0.125*exp(-(V+60)/80);
alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12);
betam = 4*exp(-(V+60)/18);
alphah = 0.07*(exp(-(V+60)/20));
betah = 1/(exp(-(V+30)/10)+1);

n = alphan/(alphan + betan);
m = alpham/(alpham + betam);
h = alphah/(alphah + betah);

Vl = V - (Iext -Gna*m^3*h*(V-Vna) - Gk*n^4*(V-Vk))/Gl;
display(Vl)
%%
%%% Question 14
disp("Question 14")
figure
Iext = 0;
y_init = [-60; 0;0;0];
equiPoint3 = fsolve(@hh,y_init,optionsfsolve);
plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)

nc1 = @(V)((Iext - Gk*equiPoint3(2)^4*(V-Vk) - Gl*(V-Vl))./(Gna*equiPoint3(4)*(V-Vna))).^1/3;
nc2 = @(V) (-0.1*(V+35)./(exp(-(V+35)/10)-1 + 10^-12))./((-0.1*(V+35)./(exp(-(V+35)/10)-1 + 10^-12))+ (4*exp(-(V+60)/18)));

V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));
equiPoint1 = equiPoint3;

for i = 0:1:20
    equiPoint1(1) = equiPoint3(1) + i ;
    deltaT = [0,1000];
    [t, y] = ode45(@hh4ode,deltaT,equiPoint1);
    plot(y(:,1),100*y(:,3))
end

figure
for i = 0:1:20
    equiPoint1(1) = equiPoint3(1) + i ;
    deltaT = [0,1000];
    [t, y] = ode45(@hh4ode,deltaT,equiPoint1);
    plot(t,y(:,1)); hold on;
end
xlim([0 10])
syms V m
J = jacobian([ (Iext -(Gna*m^3*(equiPoint3(4))*(V-Vna)) - Gk*(equiPoint3(2))^4*(V-Vk) - Gl*(V-Vl))/C, (-0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12))*(1-m) - (4*exp(-(V+60)/18))*m],[V,m]);
equi = zeros(1,2);
equi(1) = equiPoint3(1);
equi(2) = equiPoint3(3);
A = double(subs(J,[V,m],equi));
[V,D] = eig(A)

%% Question 15
disp("Question 15")
figure
for i = 8:1:12
    Iext = i;
    y_init = [-60; 0;0;0];
    equiPoint3 = fsolve(@hh,y_init,optionsfsolve);
    disp(equiPoint3)

    syms V m
    J = jacobian([ (Iext -(Gna*m^3*((0.07*(exp(-(V+60)/20)))/((0.07*(exp(-(V+60)/20)))+ (1/(exp(-(V+30)/10)+1))))*(V-Vna)) - Gk*((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))/((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))+ (0.125*exp(-(V+60)/80))))^4*(V-Vk) - Gl*(V-Vl))/C, (-0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12))*(1-m) - (4*exp(-(V+60)/18))*m],[V,m]);
    equi = zeros(1,2);
    equi(1) = equiPoint3(1);
    equi(2) = equiPoint3(3);
    A = double(subs(J,[V,m],equi));
    [V,D] = eig(A)
    if D(1,1) < 0 && D(2,2) < 0
        display(i)
        display("Stable")
    else
        display(i)
        display("Unstable")
    end
    deltaT = [0,1000];
    [t, y] = ode45(@hh4ode,deltaT,equiPoint1);
    plot(t,y(:,1)); hold on;
end    

xlim([0 100])
%%
%%% Question 16
figure(11)
fni = 0
Iext = 0;
y_init = [-60; 0;0;0];
equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3)); hold on;

fni = 0.1

equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3))

fni = 0.17

equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3))

fni = 0.2
equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3))


%% Question 17 

figure

Iext = 0;

y_init = [-60; 0;0;0];
equiPoint3 = fsolve(@hh,y_init,optionsfsolve);
plot(equiPoint3(1),100*equiPoint3(2),'b*')
disp(equiPoint3)
hinit = equiPoint3(4);
equiPoint1 = zeros(2,1);
equiPoint1(1) =  equiPoint3(1);
equiPoint1(2) =  equiPoint3(2);
for i = 0:1:20
    equiPoint1(1) = equiPoint3(1) + i ;
    deltaT = [0,1000];
    [t, y] = ode45(@hhReduced4ode,deltaT,equiPoint1);
    plot(t, y(:,1));hold on;
end
ylim([-100 100])
xlim([0 10])

%% Question 18
figure

am = @(V) -0.1 * (35+V) ./ (exp(-0.1*(35+V)) - 1);
bm = @(V) 4 * exp(-(60+V)/18);

an = @(V) 0.01 * (-(50+V)) ./ (exp(-0.1*(50+V)) - 1);
bn = @(V) 0.125 * exp(-(60+V)/80);


nc1 = @(V)((Iext - Gk*equiPoint3(4)*((am(V))./((am(V))+(bm(V)))).^3.*(V-Vna) - Gl*(V-Vl))./(Gk*(V-Vk))).^1/4;
nc2 = @(V) an(V)./(an(V)+bn(V));

V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));



for i = 0.2:0.1:0.4
    fni = i;
    equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
    %plot(equiPoint3(1),100*equiPoint3(3),'b*')
    disp(equiPoint3)
    hinit = equiPoint3(4);
    equiPoint1 = zeros(2,1);
    equiPoint1(1) =  equiPoint3(1);
    equiPoint1(2) =  equiPoint3(2);
    equiPoint1(1) = equiPoint3(1) + 10 ;
    deltaT = [0,1000];
    [t, y] = ode45(@hhReducedPara4ode,deltaT,equiPoint1);
    plot(y(:,1),100*y(:,2)); hold on;
end


%% Question 19
Iext = 0;
Gk = 36;
Gl = 0.3;
Gna = 120;
Vna = 55;
Vk = -72;
C = 1;

am = @(V) -0.1 * (35+V) ./ (exp(-0.1*(35+V)) - 1);
bm = @(V) 4 * exp(-(60+V)/18);

ah = @(V) 0.07 * exp(-(60+V)/20);
bh = @(V) 1 ./ (exp(-(30+V)/10) + 1);
 

an = @(V) 0.01 * (-(50+V)) ./ (exp(-0.1*(50+V)) - 1);
bn = @(V) 0.125 * exp(-(60+V)/80);


Iext = 0;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               am(y(1)).*(1-y(2))-bm(y(1)).*y(2); ...
               an(y(1)).*(1-y(3))-bn(y(1)).*y(3); ...
               ah(y(1)).*(1-y(4))-bh(y(1)).*y(4) ];
[T1,Y1]=ode15s(f, [0 100], [-60 0.052 0.317 0.596]);
Iext = -3;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               am(y(1)).*(1-y(2))-bm(y(1)).*y(2); ...
               an(y(1)).*(1-y(3))-bn(y(1)).*y(3); ...
               ah(y(1)).*(1-y(4))-bh(y(1)).*y(4) ];
[T2,Y2]=ode15s(f, [T1(length(T1)) 20 + T1(length(T1))], [Y1(length(Y1),1),Y1(length(Y1),2),Y1(length(Y1),3),Y1(length(Y1),4)]);
Iext = 0;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               am(y(1)).*(1-y(2))-bm(y(1)).*y(2); ...
               an(y(1)).*(1-y(3))-bn(y(1)).*y(3); ...
               ah(y(1)).*(1-y(4))-bh(y(1)).*y(4) ];
[T3,Y3]=ode15s(f, [T2(length(T2)) 150 + T2(length(T2))], [Y2(length(Y2),1),Y2(length(Y2),2),Y2(length(Y2),3),Y2(length(Y2),4)]);
figure
plot([T1;T2;T3], [Y1(:,1); Y2(:,1); Y3(:,1)])
xlabel('t'); ylabel('V(t)');
title('Membrane potential vs time');

%% Question 20

n = equiPoint3(2);
h = equiPoint3(4);

myfun1 = @(V,m) (-Gk*n.^4.*(V - Vk) - Gna*m.^3.*h.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
myfun2 = @(V,m) am(V).*(1-m)-bm(V).*m;
figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
title('Phase plane for n and h fixed at rest values');

n = Y2(length(Y2),3);
h = Y2(length(Y2),4);

myfun1 = @(V,m) (-Gk*n.^4.*(V - Vk) - Gna*m.^3.*h.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
myfun2 = @(V,m) am(V).*(1-m)-bm(V).*m;
figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
title('Phase plane for n and h fixed at end of stimulus values');

%%
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


function dYdt = hhReducedPara4ode(t,y)
    global Iext
    global Gk
    global Gl Gna
    global Vna
    global Vk
    global Vl
    global C
    global hinit
    global fni
    V = y(1);
    n = y(2);
    
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    m = alpham/(alpham+betam);
    
    h = hinit;
    
    dVdt = (Iext -Gna*m^3*(V-Vna)*fni -Gna*m^3*h*(V-Vna)*(1-fni) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    
    dndt = alphan*(1-n) - betan*n;
    
    dYdt = [dVdt;dndt]; 
end

function dYdt = hhReduced4ode(t,y)
    global Iext
    global Gk
    global Gl Gna
    global Vna
    global Vk
    global Vl
    global C
    global hinit
    V = y(1);
    n = y(2);
    
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    m = alpham/(alpham+betam);
    
    h = hinit;
    
    dVdt = (Iext -Gna*m^3*h*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    
    dndt = alphan*(1-n) - betan*n;
    
    dYdt = [dVdt;dndt]; 
end


function dYdt = hh(y)
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
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*m^3*h*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
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
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*m^3*h*(V-Vna) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end


function dYdt = hhPara(y)
    global Iext
    global Gk
    global Gl Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*m^3*(V-Vna)*fni -Gna*m^3*h*(V-Vna)*(1-fni) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end



function dYdt = hhPara4ode(t,y)
    global Iext
    global Gk
    global Gl Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*m^3*(V-Vna)*fni -Gna*m^3*h*(V-Vna)*(1-fni) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end



function [value, isterminal, direction] = myEvent(t,y)
    global equiP
    if (norm(y-equiP) < 10^-2)
        value = 0;
    else
        value = 1;
    end
    isterminal = 1;
    direction = 0;
end
