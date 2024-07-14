clear 
close all

% % State aircraft constants

K = 0.0276;         % Lift-dependent Drag Factor
CD0 = 0.0234;       % Zero-lift Drag Coefficient
S = 845;            % Wing Area, m^2
alt = 36000;        % Cruise Altitude

% Aircraft variables to be held fixed

c_tp = 1.5863*10^-5;     % Thrust-specific fuel consumption, kg/s/N


gamma = 1.4;        % Ratio of specific heats
g = 9.8065;         % acceleration due to gravity m/s/s

% Aircraft constants

MMO = 0.85;             % Initial constant operating Mach number
OEW = 277145;           % Operating Empty Weight kg
MTOW = 575000;          % Max Take-off Weight kg
MDL = MTOW - OEW;       % Max Disposable load (kg)
MSP = 84000;            % Max Structural Payload (kg)
MFC = 253983*0.9;       % Max Usable Fuel Capacity (kg)
res = 0.1; 
MRF = res*253983;       % Max Reserve Fuel (conservative) (kg)
 

% %X-29A Variables
% % State aircraft constants
% 
% K = 0.067;         % Lift-dependent Drag Factor
% CD0 = 0.069;       % Zero-lift Drag Coefficient
% S = 17.2;            % Wing Area, m^2
% alt = 33000;
% 
% % Aircraft variables to be held fixed
% 
% c_tp = 2.3*10^-5;     % Thrust-specific fuel consumption, kg/s/N
% gamma = 1.4;        % Ratio of specific heats
% g = 9.8065;        % acceleration due to gravity m/s/s
% 
% % Aircraft constants
% 
% MMO = 1.8;         % Initial constant operating Mach number
% OEW = 6260;        % Operating Empty Weight kg
% MTOW = 8074;       % Max Take-off Weight kg
% MDL = MTOW - OEW;   % Max Disposable load (kg)
% MSP = 100;        % Max Structural Payload (kg)
% MFC = 1804;        % Max Fuel Capacity (kg)
% res = 0.1; 
% MRF = res*MFC;      % Max Reserve Fuel (conservative) (kg)


[~,~,rho,a] = GetISA(alt*0.3048,0,0);

v = MMO*a;

% arbitrary point

m1 = OEW+MSP;
w = m1*g;

Cl = w/(0.5*rho*v^2*S);
Cd = CD0+K*Cl^2;

zeta_a = 0;

R_a = v/(g*c_tp)*Cl/Cd*log(1/(1-zeta_a));

% max payload range

w = MTOW*g;

Cl = w/(0.5*rho*v^2*S);
Cd = CD0+K*Cl^2;

zeta_p = ((MTOW-MSP-OEW-MRF))/MTOW;

R_p = v/(g*c_tp)*Cl/Cd*log(1/(1-zeta_p));

% max economic range

zeta_e = (MFC-MRF)/MTOW;

R_e = v/(g*c_tp)*Cl/Cd*log(1/(1-zeta_e));

% max ferry range

m = OEW+MFC;

w = m*g;

Cl = w/(0.5*rho*v^2*S);
Cd = CD0+K*Cl^2;

zeta_f = (MFC-MRF)/(OEW+MFC);

R_f = v/(g*c_tp)*Cl/Cd*log(1/(1-zeta_f));

%--------~~~~~----
PAYres = [ OEW+MSP OEW+MSP+(MFC-MRF)*0.1 MTOW-MFC OEW+MRF ];
PAY = [ OEW+MSP OEW+MSP MTOW-MFC-MRF OEW];
x = [R_a R_p R_e R_f];
y = [m1 MTOW MTOW OEW+MFC];

hold on
plot(x/1000,y/1000,'r','linewidth',2.5)
plot(x/1000,PAY/1000,'g','linewidth',2.5)
plot(x/1000,PAYres/1000,'r:','linewidth',2.5)
yline(MTOW/1000,'b:','linewidth',2)
yline(OEW/1000,'b','linewidth',2)
xline(R_p/1000,'--','color',[0 0.4470 0.7410],'linewidth',2)
xline(R_e/1000,'--','color',[0.8500 0.3250 0.0980],'linewidth',2)
xline(R_f/1000,'--','color',[0.9290 0.7940 0.1250],'linewidth',2)

    % specify a legend
legend_data             = legend('TOW','OEW+PAY','OEW+PAY+RES','MTOW','OEW','Max. payload range','Max. economic range','Max. ferry range');
legend_data.Interpreter = 'latex';
legend_data.Location    = 'southwest';
box on

% change the plot features
ax                      = gca;     % gca = get current axes, and store this information in ax
ax.FontSize             = 14;      % set the property 'FontSize' to 14
ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
ax.XAxis.Exponent       = 0;       % force the x-axis exponent
ax.YAxis.Exponent       = 0;       % force the y-axis exponent
ax.TickLabelInterpreter = 'latex'; %

% now give the data labels and set the interpreter to LaTeX
xlabel('$Range\:[km]$','Interpreter','latex');
ylabel('$TOW\:[t]$','Interpreter','latex');
title('$Payload\:range\:chart\:-\:A380$','Interpreter','latex');

% set the limits of the plot
xlim([000 15000]); % [a b] - 'a' is lower limit and 'b' is upper limit
ylim([10^2 MTOW*1.05/1000]);

% set the axis intervals
yticks(10^2:50:MTOW/1000*1.05)
%xticks(100000:1000:1.5*10^4);
%xticklabels({'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'})

% set the number of decimal places
ytickformat('%.0f'); % '%.2f' says 2 decimal places floating point

% current figure information 
fig = gcf; % gcf - get current figure, and store in 'fig'
% papersize required to print figure to pdf format
fig.PaperPositionMode = 'auto';
fig.PaperSize         = [fig.PaperPosition(3) fig.PaperPosition(4)];
% set the figure name
figure_name = 'APD_A380_rangechart';
% print figures...

 print(  fig,figure_name,'-dpdf' ,'-r600'); 

       
        







