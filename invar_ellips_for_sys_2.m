function [] =  invar_ellips_for_sys_2() 
clc;
clear ALL;
clear classes;
% addpath(genpath('C:/course/matlab/YALMIP-master/YALMIP-master'))

addpath(genpath('C:\course\matlab\Yalmip\YALMIP-master'))

A = [0, 0, 1, 0;
    0, 0, 0, 1;
    -1, 1, 0, 0;
    1, -1, 0, 0];

B_1 = [0;
       0; 
       1; 
       0];
D = [0;
     0; 
     0; 
     1];
B_2 = [0; 
       1;];
C = [0, 1, 0, 0;
    0, 0, 0, 0];


omega_0 = 0.9;
alpha =  0.1
Y = sdpvar(1,4)
Z = sdpvar(1,1)
P = sdpvar(4,4)
F = [A*P + P*A' + B_1*Y + Y'*B_1' + alpha*P + (1/alpha)*D*D' <=0, P>=0];
optimize(F)
X = sdpvar(2,1)
value(P)
K = value(Y)*inv(value(P))
G = (C+B_2*K)*value(P)*(C + B_2*K)'
% figure;
% plot(X'*inv(G)*X<= 1,[],'b'); hold on;
% xlabel('x_2','FontSize',18)
% ylabel('u','FontSize',18)
% grid on;
u_plot = [];
omega = [];
 u_step_response = [];
   t_for_u = [];
  t = 0:0.1:30;
     options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [T,mas]=ode45(@Kosh, t, [0.1, 0.1, 0.1, 0.1], options);
   for i=1:size(mas,1)
       u_plot = [u_plot K*[mas(i,1); mas(i,2); mas(i,3); mas(i,4)]];
   end
% %    figure;
%     plot(mas(:,2),u_plot,'linewidth',1,'color','r');
%  figure;
%  plot(T,mas(:,1),'linewidth',1,'color','blue');
%  grid on;
%  xlabel('t','FontSize',18)
%  ylabel('x_1','FontSize',18)
 figure;
 plot(T,mas(:,2),'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('x_2','FontSize',18)
 grid on;
%  figure;
%   plot(T,mas(:,3),'linewidth',1,'color','blue');
%  xlabel('t','FontSize',18)
%  ylabel('x_3','FontSize',18)
%  grid on;
%  figure;
%   plot(T,mas(:,4),'linewidth',1,'color','blue');
%  xlabel('t','FontSize',18)
%  ylabel('x_4','FontSize',18)
%  grid on;
 figure;
 plot(t_for_u, u_step_response, 'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('u','FontSize',18)
 grid on;
%   figure;
%  plot(T, u_plot, 'linewidth',1,'color','blue');
%  xlabel('t','FontSize',18)
%  ylabel('u','FontSize',18)
%  grid on;
function [dx] =Kosh(t, x)
dx = zeros(4,1);
u = K*[x(1); x(2); x(3); x(4)];
dx = (A+B_1*K)*[x(1); x(2); x(3); x(4)] +  D*omega_0;
u_step_response =[u_step_response u*heaviside(t)];
t_for_u =[0 t_for_u + 0.01];
end   
disp('переходные характеристики для выхода u_2')
test = abs(u_plot);
Max = max(test);
Overshoot = (Max  - test(1,end))/ test(1,end);
disp(['Перерегулирование: ', num2str(Overshoot*100), ' %']);
eps = 1e-3;
k = 0;
for i=0:size(test,2)
    k = k + 1;
if abs(test(1,k)-test(1,end)) < eps && abs(test(1,k+5) - test(1,end)) < eps
  disp(['Время регулирования: ', num2str(T(k)), ' сек.'])
  break;
end
end
err = abs(omega_0 - test(1,end));
disp(['Установившаяся ошибка: ', num2str(err)]);
disp('----------------------------------------')
disp('переходные характеристики для выхода x_2')
test = abs(mas(:,2));
Max = max(test);
Overshoot = (Max  - test(end))/ test(end);
disp(['Перерегулирование: ', num2str(Overshoot*100), ' %']);
eps = 1e-3;
k = 0;
for i=0:size(test,1)
    k = k + 1;
if abs(test(k)-test(end)) < eps && abs(test(k+5) - test(end)) < eps
  disp(['Время регулирования: ', num2str(T(k)), ' сек.'])
  break;
end
end





end