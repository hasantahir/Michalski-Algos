% This program implements Algorithm 7 from [1]
% The integrals should be of the type
% I_{\nv}(a, z, \rho) as defined in eq. (78)
clear all; close all
%% Global Parameters
% z = 1; % Second argument of the integral function
rho = 1; % 1 seems to be the optimal value
q = .4187; % Discrete increment
nu = 0;
tol = 1e-15;

% % lower limit of the integral
% a = 1.787;

% Call PE routine
% val = PE_Levin(a, tol, q);

%% Call PE_routine for sweep case (e.g. Figure 7a)
% global z i p
% z = linspace(1e-3, 1e1, 500);
% for i = 1 : length(z)
%     val(i) = PE_Levin(a, tol, q);
% end
%
% loglog(z,val)
% axis([ 1e-3 1e1 1e-3 1.5])

%% Call PE_routine for sweep case (e.g. Figure 7b)
% global z i p l
% z = .01;
% tau = 1;
% p = linspace(0, 3, 60);
% 
% for i = 1 : length(p)
%     a(i) = FindFirstZero(); % Find zeros of Bessel functions
%     if p(i) == tau
%         q = pi;
%     else
%         q = pi/abs(p(i) + tau); % This increment is very important, most of the times q = pi works
%     end
%     l = 0; % deactivate lucas transformation
%     val_1(i) = TanhSinhQuad(0, a(i), tol); % Integrate upto a through DE
%     l = 1; % set lucas transformation flag to J_plus
%     val_2(i) = PE_Levin(a(i), tol, q); % Tail through PE Levin with Lucas
%     l = 2; % set lucas transformation flag to J_minus
%     if p(i) == tau
%         q = pi;
%     else
%         q = pi/abs(p(i) - tau); % This increment is very important, most of the times q = pi works
%     end
%     val_3(i) = PE_Levin(a(i), tol, q); % Tail through PE Levin with Lucas
%     val(i) = val_1(i) + val_2(i) + val_3(i);
% end
% 
% figure (1)
% h1 = plot(p, val_1, 'linewidth',1.4,'color','black');
% hold on
% plot(p, val_1, 'ko', 'markersize',1.3);
% h2 = plot(p, val_2, 'linewidth',1.4,'color','red');
% plot(p, val_2, 'ro', 'markersize',1.3);
% h3 = plot(p, val_3, 'linewidth',1.4,'color','blue');
% plot(p, val_3, 'bo', 'markersize',1.3);
% xlabel('\rho')
% ylabel('I(z, \rho, \tau)')
% legend([h1 h2 h3],{'DE Rule', 'Lucas First', 'Lucas Second'});
% axis([ 0 3 -4 2])
% hold off
% 
% figure(2)
% plot(p, val, p, val, 'ko', 'markersize',1)
% axis([ 0 3 -4 2])
% xlabel('\rho')
% ylabel('I(z, \rho, \tau)')
% xlabel('\rho')
% ylabel('I(z, \rho, \tau)')

%% For eq. 78
a = 3.247;
% Call PE routine
tic
val = PE_Levin(a, tol, q)
toc