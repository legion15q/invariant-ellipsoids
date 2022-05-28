function [] =  FindMin_2() 
clc;
clear ALL;
clear classes;
% addpath(genpath('C:/course/matlab/YALMIP-master/YALMIP-master'))

addpath(genpath('C:\course\matlab\Yalmip\YALMIP-master'))

A = [0, 0, 1, 0;
    0, 0, 0, 1;
    -1, 1, 0, 0;
    1, -1, 0, 0];

B_1 = [0; 0; 1; 0];
D = [0; 0; 0; 1;];
B_2 = [0; 1;];
C = [0, 1, 0, 0;
    0, 0, 0, 0];

alpha = 0.4
omega_0 = 0.8;
P = sdpvar(4,4)
Y = sdpvar(1,4)
Z = sdpvar(1,1)
L = [Z Y;
    Y' P];
opt = sdpsettings('solver', 'sedumi');
F = [A*P + P*A' + B_1*Y + Y'*B_1' + alpha*P + (1/alpha)*D*D' <=0, L >=0, P>=0];
optimize(F, trace(C*P*C' + C*Y'*B_2' + B_2*Y*C' + B_2*Z*B_2'), opt)
K = value(Y)*inv(value(P))
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
  t = 0:0.1:30;
     options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [T,mas]=ode45(@Kosh, t, [0,-2, -2, 1.5], options);
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
dx = (A+B_1*K)*[x(1); x(2); x(3); x(4)] +  D*omega_0;
u_step_response =[u_step_response u*heaviside(t)];
t_for_u =[0 t_for_u + 0.01];
end   
end