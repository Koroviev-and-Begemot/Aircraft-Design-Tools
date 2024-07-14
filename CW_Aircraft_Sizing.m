clear
close all
format short

%% Constants

% Assumptions

PAX_weight = 80;                       % weight per passenger [kg]              
Baggage_weight = 20;                   % weight per bag [kg]    
We_A = 0.97;                           % Empty weight coef. []
We_exp = -0.06;                        % Empty weight exp. []
clb_credit = 200;                      % climb credit [NM]
unusablefuel = 0.05;                   % unusable fuel percentage []
ground_roll_pen = 0.7;                 % ground roll penalty []
clb_grad = 3;                          % climb gradiant 
mu = 0.3;                              % ground resistance []
LoverD_penalty_clb = 0.7;              % L/D penalty for climb []

WR_tax = 0.995;                        % Weight ratio taxi for take-off []
WR_TO = 0.99;                          % Weight ratio take-off []
WR_clb = 0.98;                         % Weight ratio climb []
WR_dec = 0.995;                        % Weight ratio decent []
WR_lnd = 0.995;                        % Weight ratio landing []
WR_tax2 = 0.999;                       % Weight ratio taxi after landing []
WR_alt_clb = 0.99;                     % Weight ratio taxi for alternate climb []
WR_alt_lnd = 0.995;                    % Weight ratio alternate landing []
WR_res = 1/(1+unusablefuel);           % Weight ratio taxi reserve []

% Aircraft constants
    
        AR = 9.5; 
        e = 0.85;
        C_D0 = 0.018;
        k = 1/(pi*e*AR);
        CL_TO = 2.1;
        CL_L = 2.3;

        a = 0.26;
        M_cr = 0.8;
        TSFC_cr = 0.429;
        TSFC_cr = TSFC_cr/3600;

        R_alt_cr = 1250;
        clb_credit = clb_credit*1852;
        R_alt_cr = R_alt_cr*1852-clb_credit;

        PAX_typ = 151;
        PAY_typ = PAX_typ*(PAX_weight+Baggage_weight);
        PAY_max = 20000;
        R_typ = 4200;
        R_typ = R_typ*1852;
        alt_cr = 35000;
        alt_cr = alt_cr/3.281;
        alt_max = 40000;                                    
        alt_max = alt_max/3.281;
        E_typ = 0.5*3600;
        g = 9.81;
        unusablefuel = 0.05;
        eng_nb = 2;
        OEI_penalty = eng_nb/(eng_nb-1);
        
        TODA = 1950;
        LDA = 1950;
          
        
%% TOW estimation

%cruise segments

CL = sqrt(C_D0/(3*k));                                             % CL for best SAR
    
CD = C_D0+k*CL^2;                                                  % CD for that CL
    
LoverD_cr = CL/CD;                                                 % L/D for cruise

[~,~,rho_cr,a_cr] = GetISA(alt_cr,0,0);                               % air density & speed of sound at cruise altitude [kg/m^3] & [m/s]

v_cr = M_cr*a_cr;                                                  % cruise velocity [m/s]                                              

WR_cr = exp(-((R_typ-clb_credit)*TSFC_cr)/(v_cr*LoverD_cr));       % weight ratio for the typical cruise segement

WR_alt_cr = exp(-(R_alt_cr*TSFC_cr)/(v_cr*LoverD_cr));             % weight ratio for the alt cruise segment 

%hold segment

CL = sqrt(C_D0/(k));                                            % CL for best Endurance

CD = C_D0+k*CL^2;                                               % CD for that CL

LoverD_end = CL/CD;                                             % L/D for hold segment

WR_hodl = exp(-(E_typ*TSFC_cr)/(LoverD_end));                   % weight ratio for hold segment

%totals

WR_tot = WR_tax*WR_TO*WR_clb*WR_dec*WR_lnd*WR_tax2*WR_alt_clb*WR_alt_lnd*WR_res*WR_cr*WR_hodl*WR_alt_cr;    % total flight weight ratio

zeta_tot = 1-WR_tot;                                                                                        % total flight block fuel ratio 


%iteration

MTOW_guess = 300000;                                                                                        % first guess for MTOW                                                                                 

dW = 1;                                                                                                     % to make sure the first iteration of the while loop runs
iter = 0;                                                                                                   % initializing iteration counter


while dW > 0.001                                                                                            % stopping criteria for the while loop - the relative difference between 2 consecutive iterations must be smaller than 0.1% 
    
 MTOW_guess_prev = MTOW_guess;                                                                              % logging previous guess  
    
 WR_e = We_A*MTOW_guess^We_exp;                                                                             % empty weight ratio
 
 MTOW_guess = PAY_typ/(1-zeta_tot-WR_e);                                                                    % new guess using previously calculated empty weight ratio
 
 dW = abs(MTOW_guess_prev-MTOW_guess)/MTOW_guess;                                                           % calculating relative change between 2 consecutive guesses
 
 iter = iter+1;                                                                                             % iteration counter
    
end

MTOW = MTOW_guess;
MTOW1 = MTOW_guess/1000 ;
disp(['MTOW: ' num2str(round(MTOW1,1)) ' T'])                                                               % displaying MTOW [T]

OEW = WR_e*MTOW;
OEW1 = OEW/1000;
disp(['OEW: ' num2str(round(OEW1,1)) ' T'])                                                                 % displaying OEW [T]

%% Payload range

% Max payload range 

WR_cr_pay = (1-(MTOW-OEW-PAY_max)/MTOW)/(WR_tax*WR_TO*WR_clb*WR_dec*WR_lnd*WR_tax2*WR_alt_clb*WR_alt_lnd*WR_res*WR_hodl*WR_alt_cr);     % calculating cruise weight ratio for max playload mission 

R_pay = (v_cr/TSFC_cr*LoverD_cr*log(1/WR_cr_pay)+clb_credit)/1852;                                                                      % calculating range for max payload mission [NM]

% Max economic range

R_eco = R_typ/1852;                                                                                                                     % range for max economic range / typical mission range [NM]

% Max ferry range

TOW = OEW+MTOW*(1-WR_tot);                                                                                                              % TOW for max ferry range [kg]

WR_cr_fer = (1-(TOW-OEW)/TOW)/(WR_tax*WR_TO*WR_clb*WR_dec*WR_lnd*WR_tax2*WR_alt_clb*WR_alt_lnd*WR_res*WR_hodl*WR_alt_cr);               % cruise weight ratio for max ferry range 

R_fer = (v_cr/TSFC_cr*LoverD_cr*log(1/WR_cr_fer)+clb_credit)/1852;                                                                      % range for max ferry range mission [NM]


%% constraints analysis

WL = linspace(0,20000,100);                                                                % vector for wing loading 

% cruise constraint

q_cr = 0.5*rho_cr*v_cr^2;                                                                  % impact pressure at cruise [Pa]

beta_cr = WR_tax*WR_TO*WR_clb;                                                             % weight ratio for cruise segment

ToverW_cr = (beta_cr/a).*((C_D0*q_cr)./(beta_cr.*WL)+((beta_cr*k)/q_cr).*WL);              % minimum T/W for cruise as a function of wing loading

% ceiling constraint

[~,~,rho_ceil,a_ceil] = GetISA(alt_max,0,0);                                                  % air density & speed of sound at max altitude

v_ceil = M_cr*a_ceil;                                                                      % velcoity at max altitude

q_ceil = 0.5*rho_ceil*v_ceil^2;                                                            % impact pressure at cruise [Pa]

beta_ceil = WR_tax*WR_TO*WR_clb*WR_cr;                                                     % weight ratio for cruise segment

ToverW_ceil = (beta_ceil/a).*((C_D0*q_ceil)./(beta_ceil.*WL)+((beta_ceil*k)/q_cr).*WL);    % minimum T/W for cruise as a function of wing loading

% TO constraint

[~,~,rho_TO,~] = GetISA (0,0,0);                                                               % air desity at take-off  [kg/m^3] 

S_TO = TODA*ground_roll_pen;                                                                % take-off distance [m]

ToverW_TO = 1.1^2/(CL_TO*rho_TO*g*S_TO).*WL;                                                % minimum T/W for take off as a fucntion of wing loading

% Climb angle constraint

gamma = atan(clb_grad/100);                                                                 % required climb angle [rad]

LoverD_clb = LoverD_end*LoverD_penalty_clb;                                                 % L/D during climb

ToverW_clb = (gamma+1/(LoverD_clb)).*OEI_penalty;                                           % minimum T/W for take-off at OEI

% Landing constraint

beta_lnd = WR_tax*WR_TO*WR_clb*WR_cr_pay*WR_dec;                                            % weight ratio for landing with max payload

S_L = LDA;                                                                                  % landing distance [m]

WL_lnd = (1/beta_lnd*(0.6*(S_L-300)*mu*g*rho_TO*CL_L))/(1.2^2);                             % max wing loading for landing 

% load factor constrait

n_max = 2.5;

ToverW_n = q_cr/a.*1./WL.*(k*n_max*beta_cr./q_cr.*WL).^2+C_D0;


%% plots

%% Constraints
figure
hold on
grid on

design_point = [ 6.6 0.235 ];

scatter(design_point(1),design_point(2),150,'h','filled');

plot(WL/1000,ToverW_cr,'color','r','linewidth',1.5)
plot(WL/1000,ToverW_ceil,'color','g','linewidth',1.5)
plot(WL/1000,ToverW_TO,'color','b','linewidth',1.5)
plot(WL/1000,ToverW_n,'color','c','linewidth',1.5)
yline(ToverW_clb,'m','linewidth',1.5)
xline(WL_lnd/1000,'y','linewidth',1.5)

area(WL/1000,ToverW_cr,'facealpha',0.1,'facecolor','r');
area(WL/1000,ToverW_ceil,'facealpha',0.1,'facecolor','g');
area(WL/1000,ToverW_TO,'facealpha',0.1,'facecolor','b');
area(WL/1000,ToverW_n,'facealpha',0.1,'facecolor','c')
fill([0 20 20 0],[0 0 ToverW_clb ToverW_clb],'m','facealpha',0.1);
fill([WL_lnd/1000 20 20 WL_lnd/1000],[0 0 10 10],'y','facealpha',0.1);
ylim([0 0.8])
xlim([2 9])
xlabel('W/S')
ylabel('T/W')
legend('Design Point','Cruise','Ceiling','Take Off','Load factor','Climb','Landing')
%title('Constraints Analysis Diagram')

hold off

% specify a legend
legend_data             = legend('Design Point','Cruise','Ceiling','Take Off','Load factor','Climb','Landing');
legend_data.Interpreter = 'latex';
legend_data.Location    = 'northeast';
box on

% change the plot features
ax                      = gca;     % gca = get current axes, and store this information in ax
ax.FontSize             = 14;      % set the property 'FontSize' to 14
ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
ax.XAxis.Exponent       = 0;       % force the x-axis exponent
ax.YAxis.Exponent       = 0;       % force the y-axis exponent
ax.TickLabelInterpreter = 'latex'; %

% now give the data labels and set the interpreter to LaTeX
xlabel('$\frac{W}{S}\:\left[\frac{kN}{m^2}\right]$ ','Interpreter','latex');
ylabel('$\frac{T}{W}$','Interpreter','latex');

% current figure information 
fig = gcf; % gcf - get current figure, and store in 'fig'
% papersize required to print figure to pdf format
fig.PaperPositionMode = 'auto';
fig.PaperSize         = [fig.PaperPosition(3) fig.PaperPosition(4)];
% set the figure name
figure_name = 'APD_Constraints';
 print(  fig,figure_name,'-dpdf' ,'-r600'); % print to pdf - ignore any warnings
% print(  fig,figure_name,'-dtiff','-r600'); % print to tiff
 %% Payload Range Chart

figure
hold on
grid on

plot([0 R_pay R_eco R_fer],[PAY_max PAY_max PAY_typ 0]/1000,'r','linewidth',1.5)
plot([0 R_eco],[PAY_typ PAY_typ]/1000,'g--','linewidth',1.5)
ylim([0 PAY_max/1000+5])
xline(R_pay,'c--','linewidth',1.5)
xline(R_eco,'c--','linewidth',1.5)
t = {' . Max Payload Range'; '  Max Economic Range'};
text([R_pay R_eco],[PAY_max PAY_typ]/1000,t)
xlabel('Range [Nautical Miles]')
ylabel('Payload [T]')
legend('Fuel limit','Typical Payload')
%title('Payload Range Diagram')

% specify a legend
legend_data             = legend('Fuel limit','Typical Payload');
legend_data.Interpreter = 'latex';
legend_data.Location    = 'northeast';
box on

% change the plot features
ax                      = gca;     % gca = get current axes, and store this information in ax
ax.FontSize             = 14;      % set the property 'FontSize' to 14
ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
ax.XAxis.Exponent       = 0;       % force the x-axis exponent
ax.YAxis.Exponent       = 0;       % force the y-axis exponent
ax.TickLabelInterpreter = 'latex'; %

% now give the data labels and set the interpreter to LaTeX
xlabel('Range [Nautical Miles]','Interpreter','latex');
ylabel('Payload [T]','Interpreter','latex');

% current figure information 
fig = gcf; % gcf - get current figure, and store in 'fig'
% papersize required to print figure to pdf format
fig.PaperPositionMode = 'auto';
fig.PaperSize         = [fig.PaperPosition(3) fig.PaperPosition(4)];
% set the figure name
figure_name = 'Payload_range';
% print figures...
print(  fig,figure_name,'-dpdf' ,'-r600'); % print to pdf - ignore any warnings


%% V-n diagram

WS_n = 6600;

rho_n = rho_TO;

v_n = linspace(0,350,350);

TW_max = 0.23;

n_max_thrust = sqrt((0.5*rho_n.*v_n.^2)./(k*WS_n).*(TW_max-0.5*rho_n.*v_n.^2*C_D0/WS_n));

n_max = 0.5*rho_n.*v_n.^2.*CL_L/WS_n;


figure(3)
plot(v_n,n_max_thrust,'linewidth',1.5,'color','r')
hold on
grid on
plot(v_n,n_max,'linewidth',1.5,'color','b')
xline(66,'-.m','linewidth',1.8)
xline(272,'-.','linewidth',1.8)
yline(2.5,'-.c','linewidth',1.8)


xa1 = linspace(66,104,39);
xb1 = flip(xa1);
x1 = [xa1 xb1];
ya1 = 0.5*ones(1,39);
yb1 = flip(n_max(66:104));


xa2 = linspace(104,113,10);
xb2 = flip(xa2);
x2 = [xa2 xb2];
ya2 = 0.5*ones(1,10);
yb2 = flip(n_max_thrust(104:113));


x = [xa1 xa2 272 272 xb2 xb1];
y = [ya1 ya2 0.5 2.5 yb2 yb1];
fill(x,y,'g','facealpha',0.1);

yline(0.5,'-.c','linewidth',1.8)

ylim([0 5])
xlim([0 350])
xlabel('Velocity [m/s]')
ylabel('Load Factor n')
legend('Thrust limit','CL max limit')

% specify a legend
legend_data             = legend('Thrust limit','CL max limit','Stall Speed','V Never Exceed','Parctical Limits','Practical V-n envelope');
legend_data.Interpreter = 'latex';
legend_data.Location    = 'northwest';
box on

% change the plot features
ax                      = gca;     % gca = get current axes, and store this information in ax
ax.FontSize             = 13;      % set the property 'FontSize' to 14
ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
ax.XAxis.Exponent       = 0;       % force the x-axis exponent
ax.YAxis.Exponent       = 0;       % force the y-axis exponent
ax.TickLabelInterpreter = 'latex'; %

% now give the data labels and set the interpreter to LaTeX
xlabel('Velocity [m/s]','Interpreter','latex');
ylabel('Load Factor n','Interpreter','latex');

% current figure information 
fig = gcf; % gcf - get current figure, and store in 'fig'
% papersize required to print figure to pdf format
fig.PaperPositionMode = 'auto';
fig.PaperSize         = [fig.PaperPosition(3) fig.PaperPosition(4)];
% set the figure name
figure_name = 'Vn';
 print(  fig,figure_name,'-dpdf' ,'-r600'); % print to pdf - ignore any warnings