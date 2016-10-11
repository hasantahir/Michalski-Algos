clear; close all
% Create Frequency Sweep variable
f =linspace(1e9,50e9,50);
% Initialize 
roots = [];
list = [];

% Time the code
tic

for i = 1 : length(f)
    % Call sweep function that finds the zeros of the programmed function
    r  = sweep(f(i));
    % Store each sweep's results aa a cell
    roots{i} = r';
    for j = 1 : length(roots{i})
        % Pair roots with its corresponding frequency
        temp = [f(i),roots{i}(j)];
        % Create a list for plotting
        list = vertcat(list, temp);
    end
end
toc

%% Plot Results
figure(1)
N = 1; % Number of colors to be used
% Use Brewer-map color scheme
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
Colord = get(gca, 'ColorOrder');
plot(list(:,1)/1e9, real(list(:,2)),'o', 'markersize',4,...
    'MarkerFaceColor',Colord(1,:));
xlabel('$\textrm{Frequency (GHz)}$','interpreter','latex')
ylabel('$\frac{k_{\rho}}{k_0}$','interpreter','latex')
% legend('Poles','Branch Point',...
%     'Location','southwest','Orientation','horizontal');
    title('TE Pole Locations');
% Decorations

box on
set(gcf,'color','white');
set(groot,'defaulttextinterpreter','latex');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,...
    'box','on',...
    'FontName','times new roman',...
    'FontSize',12);
hold off
ylim([1 3.5])
grid on;grid minor
matlab2tikz('filename',sprintf('figures/TM_pole_frequency_sweep_HEMT.tex'),'showInfo', false);