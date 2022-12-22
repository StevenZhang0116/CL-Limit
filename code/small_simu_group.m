A = 1;
% b_1 = 0.8;
step = 0.001;

omegalst = 3;
blst = 0.95;
% blst(end) = [];

[X,Y] = meshgrid(omegalst,blst);

k1lst = [];
k2lst = [];
maxk2lst = [];

for iii = 1:length(blst)
    b_1 = blst(iii);
    disp(b_1)
    for jjj = 1:length(omegalst)
        omega = omegalst(jjj);
        disp(omega)
        period = 2*pi/omega;
        
        jfunc = @(x) A + b_1*cos(omega*x);
        xfunc = @(x) 1/2*(A*T^3/3+b_1/omega^3*((2-omega^2*T^2)*sin(omega*(x-T))+2*omega*T*cos(omega*(x-T))-2*sin(omega*x)));
        
        x_a = 1;
        
        j0 = jfunc(0);
        va = j0 * 9/(4 * sqrt(2));
        va = va^(2/3);
        T0 = 3 / sqrt(2*va);
        tslot = 0:step:period;
        Tslot = [];
        Tslot(1) = T0;
        
        for i = 1:length(tslot)-1
            t = tslot(i);
            T = Tslot(i);
           
            integral = 0;
            for k=t-T:step:t
                integral = step*(t-k)*jfunc(k)+integral;
            end
        
            dTdt = 1-2/(T^2*jfunc(t-T))*integral;
            nextT = Tslot(i) + step*dTdt;
            Tslot(i+1) = nextT;
        end
        
        
        
        vfunc_solution = [];
        jfunc_solution = [];
        for i=1:length(tslot)
            t = tslot(i);
            T = Tslot(i);
            vfunc = @(x) 1/2*(A^2*T^4/4 + 2/omega*b_1*x_a*sin(omega*(x))-A*b_1*T^3/omega*sin(omega*(x-T))...
                +2*A*b_1/(omega^4)*((T^2*omega^2-2)*cos(omega*(x-T))+2*omega*T*sin(omega*(x-T))...
                +2*cos(omega*x))-b_1^2/(8*omega^4)*((2*T^2*omega^2-1)*cos(2*omega*(x-T))...
                +2*T*omega*sin(2*omega*(x-T))+cos(2*omega*x)));
            vfunc_solution(i) = vfunc(t);
            jfunc_solution(i) = jfunc(t);
        end
        
        v1 = 0;
        v2 = 0;
        for i=1:length(tslot)
            v1 = v1 + vfunc_solution(i);
            v2 = v2 + vfunc_solution(i)^(3/2);
        end
        v1 = v1/length(vfunc_solution);
        v2 = v2/length(vfunc_solution);
        v2 = v2^(2/3);
        j1 = 4*sqrt(2)/9*(v1^(3/2));
        j2 = 4*sqrt(2)/9*(v2^(3/2));
        
        j_avg = mean(jfunc_solution);
        
        k1 = norm(j_avg/j1) - 1;
        k2 = norm(j_avg/j2) - 1;
        k1lst(jjj) = k1/100;
        k2lst(jjj) = k2/100;
        disp(k2)

        maxk2lst(iii,jjj) = k2/100;
    end
    
end

% figure()
% hold on
% plot(omegalst,k1lst,'LineWidth',2,'LineStyle','-')
% plot(omegalst,k2lst,'LineWidth',2,'LineStyle','-')
% legend('r1','r2','FontSize',12)
% xlabel('\omega','FontSize',20)
% ylabel('Value','FontSize',20)
% hold off

% figure(1)
% contourf(omegalst,blst,maxk2lst,'ShowText','on')
% xlabel('\omega','FontSize',20)
% ylabel('b','FontSize',20)
% saveas(gcf,'contour.jpg')

% figure()
% plot(blst,maxk2lst,'LineWidth',2,'LineStyle','-')
% xlabel('b','FontSize',20)
% ylabel('max(k_2)','FontSize',20)

