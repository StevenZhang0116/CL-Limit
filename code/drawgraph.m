xval = load('tslot.mat').tslot;
jt = load('jt.mat').jfunc_solution;
vt = load('vt.mat').vfunc_solution;

figure()
hold on
plot(xval,jt)
plot(xval,vt)
legend('j(t)','V(t)');
xlabel('t');
ylabel('V');
hold off

