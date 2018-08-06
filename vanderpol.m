global mu
mu = 100; 

y_init = [1; 0];

deltaT = [0, 500];
tic;
[t, y] = ode23(@systemDE,deltaT,y_init);
toc;
subplot(1,2,1)
plot(t,y)
title('Time history')

subplot(1,2,2)
plot(y(:,1),y(:,2))
title('Phase plane plot')

function dYdt = systemDE(t,y)
    global mu
    dYdt = y;
    dYdt(1) = mu*y(2);
    dYdt(2) = mu*(1-y(1).^2)*y(2) - y(1)*1/(mu);
end
