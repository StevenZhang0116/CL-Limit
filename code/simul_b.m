%% Parameters Settings
close all
clc
clear

omega = 0.3;
A = 1;
blst = [0.4,0.8];
period = 2*pi/omega;
T = 1.7;

k1lst = [];
k2lst = [];

figure()
hold on

for j = 1:length(blst)
    b_1 = blst(j);

    epsilon = 8.85418782*10^-12;
    
    jfunc = @(x) A + b_1*cos(omega*x);
    xfunc = @(x) 1/2*(A*T^3/3+b_1/omega^3*((2-omega^2*T^2)*sin(omega*(x-T))+2*omega*T*cos(omega*(x-T))-2*sin(omega*x)));
    
    xslot = 0:0.001:period;
    y = xfunc(xslot); % x_a set

    % x_a = 1;
    y0 = []; % j set
    y1 = []; % V set
    
    for k = 1:length(xslot)
        % x_a = 1;
        x_a = y(k);
    
        vfunc = @(x) 1/2*(A^2*T^4/4 + 2/omega*b_1*x_a*sin(omega*(x))-A*b_1*T^3/omega*sin(omega*(x-T))...
        +2*A*b_1/(omega^4)*((T^2*omega^2-2)*cos(omega*(x-T))+2*omega*T*sin(omega*(x-T))...
        +2*cos(omega*x))-b_1^2/(8*omega^4)*((2*T^2*omega^2-1)*cos(2*omega*(x-T))...
        +2*T*omega*sin(2*omega*(x-T))+cos(2*omega*x)));
    
        y0(k) = jfunc(xslot(k));
        y1(k) = vfunc(xslot(k));
    end
    
    plot(xslot,y1,'LineWidth',2,'LineStyle','-')
  
end

hold off
legend('j','V','FontSize',12)
xlabel('t','FontSize',20)
ylabel('V','FontSize',20)
xlim([0 period])
ylim([0 4])



