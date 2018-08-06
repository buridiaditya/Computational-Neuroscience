global mu


y_init = [1; 0];
deltaT = [0, 100];
mu = 0.1; 
tic;
[t, y] = ode45(@systemDE,deltaT,y_init);
toc;
tic;
[t, y] = ode15s(@systemDE,deltaT,y_init);
toc;

subplot(2,3,1)
plot(t,y)
title('Time history u = 0.1 ode15s')

subplot(2,3,4)
plot(y(:,1),y(:,2))
title('Phase plane plot u = 0.1 ode15s')


y_init = [1; 0];
deltaT = [0, 100];
mu = 1; 
tic;
[t, y] = ode45(@systemDE,deltaT,y_init);
toc;
tic;
[t, y] = ode15s(@systemDE,deltaT,y_init);
toc;

subplot(2,3,2)
plot(t,y)
title('Time history u = 1 ode15s')

subplot(2,3,5)
plot(y(:,1),y(:,2))
title('Phase plane plot u = 1 ode15s')

y_init = [1; 0];
deltaT = [0, 1000];
mu = 100; 
tic;
[t, y] = ode45(@systemDE,deltaT,y_init);
toc;
tic;
[t, y] = ode15s(@systemDE,deltaT,y_init);
toc;

subplot(2,3,3)
plot(t,y)
title('Time history u = 100 ode15s')

subplot(2,3,6)
plot(y(:,1),y(:,2))
title('Phase plane plot u = 100 ode15s')

function dYdt = systemDE(t,y)
    global mu
    dYdt = y;
    dYdt(1) = mu*y(2);
    dYdt(2) = mu*(1-y(1).^2)*y(2) - y(1)*1/(mu);
end
