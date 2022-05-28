function [] = LQG()
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


nx = 4;    %Number of states
ny = 2;    %Number of outputs
sys = ss(A,B_1,C,B_2);
QXU = [1 0 0 0 0;
      0 1 0 0 0;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 1];
QWV = [1 0 0 0 0 0;
      0 1 0 0 0 0;
      0 0 1 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1];
QI = eye(ny);

         
KLQG1 = lqg (sys, QXU, QWV, QI, '1dof')

K = [-3.15     -0.6135      -2.881      -3.202]



 
  u_plot = [];
  t = 0:0.1:30;
     options = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [T,mas]=ode45(@Kosh, t, [0.5, 0.5, 0.5, 0.5], options);
   for i=1:size(mas,1)
       u_plot = [u_plot K*[mas(i,1); mas(i,2); mas(i,3); mas(i,4)]];
   end
 figure;
 plot(T,mas(:,1),'linewidth',1,'color','blue');
 grid on;
 xlabel('t, [ñ]','FontSize',18)
 ylabel('x_1, [ì/c]   ','FontSize',18)
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
 plot(T, u_plot, 'linewidth',1,'color','blue');
 xlabel('t','FontSize',18)
 ylabel('u','FontSize',18)
 grid on; 
%   
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
 
 
    
 function [dx] =Kosh(t, x)
dx = zeros(4,1);
u = K*[x(1); x(2); x(3); x(4)];
dx = (A+B_1*K)*[x(1); x(2); x(3); x(4)] +  D*[omega_0_1; omega_0_2];
end   
 
end