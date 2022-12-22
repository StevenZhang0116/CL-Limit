close all
clear

t_max = 2;
step = 2 * 1e-3;
tint = 0:step:t_max;

Qlst = 0.50:0.001:1.50;

k1lst = [];
k2lst = [];

for nnn = 1:length(Qlst)
    Q = Qlst(nnn);
    disp(Q)
    newj = @(x) (x<Q);
    % A = 1;
    % b_1 = 0.8;
    % omega = 1.5;
    % period = 2*pi/omega;
    % tint = 0:step:period;
    % newj = @(x) A + b_1*cos(omega*x);
    
    j0 = newj(0);
    va = 9*j0/(4*sqrt(2));
    va = va^(2/3);
    newT0 = 3/sqrt(2*va);
    
    Tint = [];
    Tint(1) = newT0;
    Vint = [];
    jaint = [];
    
    % calculate dTdt
    for i=1:length(tint)-1
        t = tint(i);
        T = Tint(i);
        integral = 0;
        for k=t-T:step:t
            ind = 1;
            if k == t-T || k == t
                ind = 1/2;
            end
    
            integral = ind * step*(t-k)*newj(k)+integral;
        end
        dTdt = 1-2/(T^2*newj(t-T))*integral;
        nextT = Tint(i) + step*dTdt;
        Tint(i+1) = nextT;
    end
    
    newTint = [];
    lenn = double(Q/step)+1;
    for jjj = 1:length(tint)
        if jjj <= lenn
            newTint(end+1) = 2.03 + (jjj-1)*1.3/lenn;
        else
            newTint(end+1) = 3.33;
        end
    end
    
    % Tint = newTint;
    
    
    % calculate dVadt
    % for ii=1:length(tint)
    %     t = tint(ii);
    %     disp(t)
    %     T = Tint(ii);
    %     integral = 0;
    % 
    %     for k=t-T:step:t
    %         part = (t-k)^2*newj(k);
    %         small = 0;
    %         for kk=k:step:t
    %             small = small + step*newj(kk);
    %         end
    %         part = part * small;
    %         integral = integral + step*part;
    %     end
    %     Vint(ii) = 1/2 * integral;
    % end
    
    newVint = [];
    V0 = 0.43;
    Vmax = 0.84 * (Q/0.7)^1.1;
    lenn1 = double((Q/2)/step)+1;
    lenn2 = double((Q)/step)+1;
    for jjj=1:length(tint)
        if jjj <= lenn1
            newVint(end+1) = V0 + (jjj-1)*(Vmax-V0)/lenn1;
        elseif jjj > lenn1 && jjj <= lenn2
            newVint(end+1) = Vmax - (jjj-lenn1)*(Vmax-V0)/(lenn2-lenn1);
        else
            newVint(end+1) = V0;
        end
    end
    
    Vint = newVint;
    
    v1 = 0;
    v2 = 0;
    for i=1:length(tint)
        v1 = v1 + Vint(i);
        v2 = v2 + Vint(i)^(3/2);
    end
    v1 = v1/length(tint);
    v2 = v2/length(tint);
    v2 = v2^(2/3);
    j1 = 4*sqrt(2)/9*v1^(3/2);
    j2 = 4*sqrt(2)/9*v2^(3/2);
    
    j_avg = Q/t_max;
    
    k1 = norm(j_avg/j1) - 1;
    k2 = norm(j_avg/j2) - 1;

    k1lst(end+1) = k1
    k2lst(end+1) = k2

end

figure()
hold on
plot(Qlst,k2lst,'LineWidth',2,'LineStyle','-')
legend('r2','FontSize',12);
xlabel('Q','FontSize',20)
ylabel('value','FontSize',20)




% figure();
% hold on
% fplot(newj,[0 2],0.001,'LineStyle','-.');
% plot(tint,Vint,'LineWidth',2,'LineStyle','-');
% legend("j","V",'FontSize',12)
% xlabel('t','FontSize',20)
% ylabel('V','FontSize',20)
% % ylim([0 1.1])
% hold off




