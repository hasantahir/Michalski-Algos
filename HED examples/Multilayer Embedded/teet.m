close all;
% % D = @(z) z.^50 +z.^12 - 5*sin(20*z).*cos(12*z) - 1;
% % EM constants
% mu0 = 4*pi*1e-7;
% ep0 = 8.854e-12;
% 
% lambda = 852e-9;
% c = 3e8;
% ep2 = -33.12 + 1.700i;
% ep1 = 1;
% omega = 2*pi*c./lambda;
% 
% k1 = omega * sqrt(mu0* ep0* ep1);
% k2 = omega * sqrt(mu0* ep0* ep2);
% k2 = conj(k2);
% 
% kz1 = @(kp) sqrt( k1^2 - kp.^2);
% kz2 = @(kp) sqrt( k2^2 - kp.^2);
% 
% D = @(kp) kz2(kp)./ep2 + kz1(kp)./ep1;
% D = @(z) log(z) -z + 1.195281 + .536353*1i;
% % D = @(z) (sqrt(z.^2 - 1));
% % p = linspace(k1/1.1,k1*1.1,30000);
% root = [];
% % for i = 1 : length(p)
% %     r = newtzero(D,p(i)+.001i);
% %     root = vertcat(root,r);
% % end
% % % % Sort the array
% % % root = sort(root);
% % % % % Clean up roots by weeding out too close values
% % % % if ~isempty(root)
% % % %     cnt = 1;  % Counter for while loop.
% % % %     
% % % %     while ~isempty(root)
% % % %         vct = abs(root - root(1)) < 1e-9; % Minimum spacing between roots.
% % % %         C = root(vct);  %C has roots grouped close together.
% % % %         [idx,idx] = min(abs(D(C)));  % Pick the best root per group.
% % % %         rt(cnt) = C(idx); %  Most root vectors are small.
% % % %         root(vct) = []; % Deplete the pool of roots.
% % % %         cnt = cnt + 1;  % Increment the counter.
% % % %     end
% % % %     root = sort(rt).';  % return a nice, sorted column vector
% % % % end
% % root = newtzero(D,k1 + 1000i,1000);
% root = newtzero(D,1+1i,30);
% % hold on
% % plot(real(k2),imag(k2),'o')
% % % plot(real(root),imag(root),'s')
%  hold on
% % plot(real(k1),imag(k1),'d')
% % % hold off
% % plot(real(root),imag(root),'o')
% % hold on
% % semilogx(p,real(D(p)))
% % semilogx(p,imag(D(p)))
% % plot(real(k1),imag(k1),'o')


% % Example frmo Delnitiz paper
% A = -.19435;
% B = 1000.41;
% C = 522463;
% T = .005;
% 
% f = @(z) A*z.^2 + B * exp(-T*z) + C;
% p = linspace(-15000,5000,3000);
% q = linspace(-15000,15000,3000);
% root = [];
% for i = length(p)
%     for j = length(q)
%     Z = p(i) + q(j);
%     r = newtzero(f,Z);
%     root = vertcat(root,r);
%     end
% end
% % root = newtzero(f,-15000,100);
% plot(real(root),imag(root),'s')
 options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-15, 'MaxIterations',1e6, 'Algorithm','trust-region-dogleg',...
            'OptimalityTolerance', 1e-15,'StepTolerance', 1e-4, 'MaxFunctionEvaluations',1e6);
x0 = linspace(0, k3,1e1);
root = [];
for i = 1 : length(x0)
    [r,fval,exitflag,output] = fsolve(fun,x0(i),options);  
    root = vertcat(root,r);
end

figure(3)
N = 2; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% Obtain the colors of the plot to apply it on marker faces
Colord = get(gca, 'ColorOrder');
%
plot(real(D(root)), 's', 'markersize',4,...
    'MarkerFaceColor',Colord(1,:));
hold on
plot(imag(D(root)), 's', 'markersize',4,...
    'MarkerFaceColor',Colord(2,:));

figure(1)
N = 5; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
Colord = get(gca, 'ColorOrder');

plot(real(root)/k1, (imag(root)/k1), 's', 'markersize',4,...
    'MarkerFaceColor',Colord(1,:));
hold on
plot(real(k2)/k1 , imag(k2)/k1, 'd', 'markersize',4,...
    'MarkerFaceColor',Colord(2,:));
plot(real(k1)/k1 , imag(k1)/k1, 'd', 'markersize',4,...
    'MarkerFaceColor',Colord(4,:));
plot(real(k3)/k1 , imag(k3)/k1, 'd', 'markersize',4,...
    'MarkerFaceColor',Colord(5,:));