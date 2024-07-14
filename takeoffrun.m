clear
close all
g = 9.81;

% airport characteristics
delta_t = 0;
alt = [ 0 1000 2000]';
gamma = 0;
mu = 0.05;

[~,~,rho,~] = GetISA(alt,delta_t,0);

% aircraft characteristics
T = 0.9*335*10^3*4;
m = linspace(426*10^3,575*10^3,100);
W = m.*g;
S = 845;
CD0 = 0.0437;
K = 0.0192;
CL = (CD0/(3*K))^0.5;
CD = CD0 + K.*CL.^2;
CLmax = 1.5;

Vlo = 1.1*sqrt((2.*W)./(rho*S*CLmax));
A = T./W-mu*cos(gamma)-sin(gamma);
B = (rho*S*CL)./(2.*W).*(CD./CL-mu);

sg1 = -1./(2.*B*g).*(log(A-B.*Vlo.^2)-log(A));

%disp(['the required take off run is ',num2str(round(sg)),' m'])

T = 48.9*10^3;
m2 = linspace(7.17*10^3,8.07*10^3,100);
W = m2.*g;
S = 17.2;
CD0 = 0.0375;
K = 0.067;
CL = (CD0/(3*K))^0.5;
CD = CD0 + K.*CL.^2;
CLmax = 1.5;

Vlo = 1.1*sqrt((2.*W)./(rho*S*CLmax));
A = T./W-mu*cos(gamma)-sin(gamma);
B = (rho*S*CL)./(2.*W).*(CD./CL-mu);

sg2 = -1./(2.*B*g).*(log(A-B.*Vlo.^2)-log(A));

MDoverMDL = linspace(0.5,1,100);

grid on
hold on

n1 = {'-0 m','-1000 m','-2000 m'}; 
n2 = {'0 m','-1000 m','/2000 m'}; 

plot(MDoverMDL,sg1,'lineWidth',1.7)
plot(MDoverMDL,sg2,'linewidth',1.7)

% specify a legend
legend_data             = legend('A380 at 0 m','A380 at 1000 m','A380 at 2000 m','X-29A at 0m','X-29A at 1000m','X-29A at 2000m');
legend_data.Interpreter = 'latex';
legend_data.Location    = 'northwest';
box on

% change the plot features
ax                      = gca;     % gca = get current axes, and store this information in ax
ax.FontSize             = 14;      % set the property 'FontSize' to 14
ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
ax.XAxis.Exponent       = 0;       % force the x-axis exponent
ax.YAxis.Exponent       = 0;       % force the y-axis exponent
ax.TickLabelInterpreter = 'latex'; %

% now give the data labels and set the interpreter to LaTeX
xlabel('$\frac {DL}{MDL}$','Interpreter','latex');
ylabel('Distance to take-off [m]','Interpreter','latex');


% set the number of decimal places
ytickformat('%.0f'); % '%.2f' says 2 decimal places floating point

% current figure information 
fig = gcf; % gcf - get current figure, and store in 'fig'
fig.PaperPositionMode = 'auto';
fig.PaperSize         = [fig.PaperPosition(3) fig.PaperPosition(4)];
% set the figure name
figure_name = 'APD_TO_plus15';
 print(  fig,figure_name,'-dpdf' ,'-r600');