function [] =  min_ellips_for_sys_3() 
clc;
clear ALL;
clear classes;
% addpath(genpath('C:/course/matlab/YALMIPe-master/YALMIP-master'))

addpath(genpath('C:\course\matlab\Yalmip\YALMIP-master'))

A = [0, 0, 1, 0;
    0, 0, 0, 1;
    -1, 1, 0, 0;
    1, -1, 0, 0];

B_1 = [0; 
       0; 
       1; 
       0];
D = [0, 0; 
     0, 0;
     1, 0 ;
     0, 1;];
B_2 = [0; 
       1;];
C = [0, 1, 0, 0;
    0, 0, 0, 0];

alpha = 0.45
omega_0_1 = 0.8;
omega_0_2 = 0.8;
P = sdpvar(4,4)
Y = sdpvar(1,4)
Z = sdpvar(1,1)
L = [Z Y;
    Y' P];
opt = sdpsettings('solver', 'sedumi');
F = [A*P + P*A' + B_1*Y + Y'*B_1' + alpha*P + (1/alpha)*D*D' <=0, L >=0, P>=0];
optimize(F, trace(C*P*C' + C*Y'*B_2' + B_2*Y*C' + B_2*Z*B_2'), opt)
K = value(Y)*inv(value(P))
%K_lqr = [-0.6285 , -0.0786 , -1.1212 , -0.8069 ];
%K = K_lqr;
x_2 = sdpvar(1,1);
u_ = sdpvar(1,1);
X = [x_2;u_];
G = C*value(P)*C' + C*value(Y)'*B_2' + B_2*value(Y)*C' + B_2*value(Z)*B_2'
figure;
plot(X'*inv(G)*X<= 1,[],'b'); hold on;
xlabel('x_2','FontSize',18)
ylabel('u','FontSize',18)
grid on;
u_plot = [];
omega = [];
 u_step_response = [];
   t_for_u = [];
  t = 0:0.1:40;
     options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [T,mas]=ode45(@Kosh, t, [0.5, 0.5, 0.5, 0.5], options);
   for i=1:size(mas,1)
       u_plot = [u_plot K*[mas(i,1); mas(i,2); mas(i,3); mas(i,4)]];
   end
    plot(mas(:,2),u_plot,'linewidth',1,'color','r');
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
 plot(t_for_u, u_step_response, 'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('u','FontSize',18)
 grid on;
  figure;
 plot(T, u_plot, 'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('u','FontSize',18)
 grid on;
function [dx] =Kosh(t, x)
dx = zeros(4,1);
u = K*[x(1); x(2); x(3); x(4)];
dx = (A+B_1*K)*[x(1); x(2); x(3); x(4)] +  D*[omega_0_1; omega_0_2];
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
if abs(test(1,k)-test(1,end)) < eps && abs(test(1,k+2) - test(1,end)) < eps
 disp(['Время регулирования: ', num2str(T(k)), ' сек.'])
  break;
end
end
err = abs(omega_0_1+omega_0_2 - test(1,end));
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
if abs(test(k)-test(end)) < eps && abs(test(k+2) - test(end)) < eps
  disp(['Время регулирования: ', num2str(T(k)), ' сек.'])
  break;
end
end

y = mas(:,2);
 y = y.^2;

 I = trapz(T,abs(y(:) - y(end)))
 figure;
 %plot(T,y);
 area(T,y, y(end))
  xlabel('t','FontSize',18)
 ylabel('x_2^2','FontSize',18)

 grid on;
 y = u_plot(:);
 y = y.^2;
 I = trapz(T,abs(y(:) - y(end)))
 figure;
  %plot(T,y);
 area(T,y, y(end))
  xlabel('t','FontSize',18)
 ylabel('u^2','FontSize',18)
 
 grid on;
 

end