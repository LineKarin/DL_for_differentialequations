%tspan = [0 10];
a=0;
b=1;
N = 1003;
tspan = linspace(a,b,N);
u0 = 2;
U0 = [u0;0];
h = (b-a)/(N-1);
s = 4;
AT = [0,0,0,0; 1/2,0,0,0; 0,1/2,0,0;0,0,1,0]';
b = [1/6; 1/3; 1/3; 1/6];
c = [0, 1/2, 1/2, 1];
C = 1;

[t,u,FuncEval]=ExplicitRungeKutta(@funky,tspan,U0,h,s,AT,b,c,C);

% Coarsen
tc = [0 t(t==1/3) t(t==2/3) t(end)];
Uc = [u(1,:) u(t==1/3,:) u(t==2/3,:) u(end,:)]';

nn = length(Uc);
uc = Uc(1:2:nn);
vc = Uc(2:2:nn);

figure(1) 
clf
hold on
plot(t,u,'-','linewidth',2)
plot(t,u0*cos(C*t),'b--','linewidth',2)
plot(t,-C*u0*sin(C*t),'r--','linewidth',2)
plot(tc,uc,'kx',tc,vc,'kx','linewidth',2)
lgd1=legend({'u(t)','v(t)','u-ex','v-ex','uc','vc'},'location','bestoutside');
lgd1.FontSize=16;
hold off

%% Create data sets
I = 10;
u0i = rand(I,1);
a=0;
b=1;
N = 1003;
tspan = linspace(a,b,N);
h = (b-a)/(N-1);
s = 4;
AT = [0,0,0,0; 1/2,0,0,0; 0,1/2,0,0;0,0,1,0]';
b = [1/6; 1/3; 1/3; 1/6];
c = [0, 1/2, 1/2, 1];
C = 1;
steps = 4;
Uc = zeros(2*steps,I);

for i = 1:I
    u0 = u0i(i);
    U0 = [u0;0];
    
    [t,u,FuncEval]=ExplicitRungeKutta(@funky,tspan,U0,h,s,AT,b,c,C);

    % Coarsen
    tc = [0 t(t==1/3) t(t==2/3) t(end)];
    Uc(:,i) = [u(1,:) u(t==1/3,:) u(t==2/3,:) u(end,:)]';
end

TData = Uc(1:2:2*steps,:);

% Create data on csv-file
csvwrite('TrainData.csv',TData)


