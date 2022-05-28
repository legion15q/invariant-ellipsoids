function [] = worst()
clc;
clear ALL;
% 
% addpath(genpath('C:/course/matlab/YALMIP-master/YALMIP-master'))

addpath(genpath('C:\course\matlab\Yalmip\YALMIP-master'))
D=[0;
    0;
    0;
    1]


C=[0 0 1 0;
    0 0 0 1]



 A= [0, 0, 1, 0;
     0, 0, 0, 1;
     -2, 1, -0.2, 0;
     1, -1, 0, -0.2]

%Интегрирование системы

%задаем начальное возмущение
omega_0 = 0.8;
omega = [];
  t = 0:0.1:500;
  options = odeset('RelTol',1e-9,'AbsTol',1e-9);
  [T,mas]=ode45(@Kosh, t, [1, 1, 4, 4], options);

function [dx] =Kosh(t, x)
dx = zeros(4,1);
dx = A*[x(1); x(2); x(3); x(4)] + D*omega_0;
end   
P_trace = sdpvar(4,4);
[~,D_1] = eig(A);
alpha = -2*max(real(diag(D_1))) / 2
opt = sdpsettings('solver', 'sedumi');
F = [ A*P_trace + P_trace*A' + alpha*P_trace + (1/alpha)*D*D' == 0, P_trace>=0];
optimize(F, trace(C*P_trace*C'));
P_trace = value(P_trace)


 

% считаем наихудшее возмущение
for i=1:size(mas,1)
omega = [omega (D'*inv(P_trace) *[mas(i,1); mas(i,2);mas(i,3);mas(i,4)])/ norm(D'*inv(P_trace) *[mas(i,1); mas(i,2);mas(i,3);mas(i,4)])];
end


 options = odeset('RelTol',1e-9,'AbsTol',1e-9);
  [T,mas]=ode45(@Kosh_2, t, [1, 1, 4, 4], options);


function [dx] =Kosh_2(t, x)
dx = zeros(4,1);
if (mod(t, floor(t)) <= 0.5 )
   dx = A*[x(1); x(2); x(3); x(4)] + D*omega(1,floor(t)+1); 
else
    dx = A*[x(1); x(2); x(3); x(4)] + D*omega(1,ceil(t)+1);  
end

end   


 figure;
 plot(T,mas(:,1),'linewidth',1,'color','blue');
 grid on;
 xlabel('t','FontSize',18)
 ylabel('x_1','FontSize',18)
 figure;
 plot(T,mas(:,2),'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('x_2','FontSize',18)
 grid on;
 figure;
  plot(T,mas(:,3),'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('x_3','FontSize',18)
 grid on;
 figure;
  plot(T,mas(:,4),'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('x_4','FontSize',18)
 grid on;
 figure;
y_1 = sdpvar(1,1)
y_2 = sdpvar(1,1)
X = [y_1;y_2]
plot(X'*inv(value(C*P_trace*C'))*X <= 1,[],'r'); hold on;
grid on;
xlabel('y_1','FontSize',18)
ylabel('y_2','FontSize',18)
plot(mas(:,3),mas(:,4),'linewidth',1.2,'color','b' )




end