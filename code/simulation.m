%% Parameters Settings
A = 1;
b_1 = 0.4;
% b_2 = 0.8;
period = 2*pi/omega;
T = 1.8;

omega = 0.3;


k1lst = [];
k2lst = [];

% for j = 1:length(omegalst)
    epsilon = 8.85418782*10^-12;
    
    jfunc = @(x) A + b_1*cos(omega*x);
    xfunc = @(x) 1/2*(A*T^3/3+b_1/omega^3*((2-omega^2*T^2)*sin(omega*(x-T))+2*omega*T*cos(omega*(x-T))-2*sin(omega*x)));
    
    xslot = 0:0.001:period;
    y = xfunc(xslot); % x_a set
%     plot(xslot,y)
    % save('xslot.mat',"xslot");
    
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
    
    % save('v_0.8.mat',"y1");
    
    % y1 = vfunc(x);
    
    %% Draw Graph
    
    figure()
    % hold on
    % plot(xslot,y0);
    plot(xslot,y1);
    % hold off
    % 
    % legend({'j(t)', 'V(t)'}, 'FontSize', 12);
    pause
    
    %% Analyze min, max, and average
    
    % inversevfunc = @(x) -1*(1.2+1/2*(A^2*(T^4)/4 + 2/omega*b_1*x_a*sin(omega*(x-t))-A*b_1*T^3/omega*sin(omega*(x-t-T))...
    %     +2*A*b_1/(omega^4)*((T^2*omega^2-2)*cos(omega*(x-t-T))+2*omega*T*sin(omega*(x-t-T))...
    %     +2*cos(omega*(x-t)))-b_1^2/(8*omega^4)*((2*T^2*omega^2-1)*cos(2*omega*(x-t-T))...
    %     +2*T*omega*sin(2*omega*(x-t-T))+cos(2*omega*(x-t)))));
    % 
    % [xmin, minval] = fminbnd(vfunc,0,period);
    % [xmax, maxval] = fminbnd(inversevfunc,0,period);
    % maxval = -maxval;
    
    v1 = 0;
    v2 = 0;
    for i=1:length(xslot)
        v1 = v1 + y1(i);
        v2 = v2 + y1(i)^(3/2);
    end
    v1 = v1/length(xslot);
    v2 = v2/length(xslot);
    v2 = v2^(2/3);
    j1 = 4*sqrt(2)/9*v1^(3/2);
    j2 = 4*sqrt(2)/9*v2^(3/2);
    
    j_avg = mean(y0);
    
    k1 = norm(j_avg/j1) - 1
    k2 = norm(j_avg/j2) - 1
    k1lst(j) = k1;
    k2lst(j) = k2;
% end

figure()
hold on
plot(omegalst,k1lst);
plot(omegalst,k2lst);
xlabel('\omega');
ylabel('r value')
legend({'r_1', 'r_2'}, 'FontSize', 12);

hold off


