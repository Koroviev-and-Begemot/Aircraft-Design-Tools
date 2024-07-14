clear

g = 9.8065;

gamma = 1.4;
roh_k = 8.075;
roh_GTL = 7.784;

% asumming max fuel + max payload at MTOW

% OEW = 277145;
% ctp = 1.5863*10^-5;
% pf = 0.077;
% min_zeta = 0.22;
% max_zeta = 0.44;
% S = 845;
% CD_0 = 0.0239;
% K = 0.0247;
% M_MO = 0.85;
% dT = 0;
% dL = 0;

OEW = 6260;
ctp = 2.95*10^-5;
pf = 0.012;
min_zeta = 0.11;
max_zeta = 0.22;
S = 17.2;
CD_0 = 0.0375;
K = 0.067;
M_MO = 1.8;
dT = 0;
dL = 0;


zeta = linspace(min_zeta,max_zeta,300);
C_LmaxSAR = sqrt(CD_0/(3*K));
roh_fuel=0.84;



 for z=1:length(zeta)  
     Cd = CD_0+K*C_LmaxSAR^2;
    
     TOW = OEW./(1 - pf - zeta);
     
     p1= (2.*TOW*g)./(gamma*M_MO^2*S*C_LmaxSAR);
     
     [p,T,~,~] = GetISA(20000,dT,dL);
     
     [H1,a1,T1] = GetISAinverse(p1,dT,dL);
     
     v1(z)= a1(z)*M_MO;

     R1(z) = ((v1(z)/(g*ctp))*(C_LmaxSAR./Cd)*log(1/(1-zeta(z))))/10^3;
     
     FB1(z) = (zeta(z)*TOW(z))/R1(z);
     
     FE(z) = 100*(100*(FB1(z)/roh_fuel)/(pf*TOW(z)));
 
          
 end
 figure(1);
 grid on;
 plot(R1,FB1,'b','Linewidth', 2);
% set(gca,'xlim',[0 15000],'ylim',[0 4]);
 xlabel('Aircraft Range,(km)');
 ylabel('Fuel Burn, (kg/km)');
 legend('$C_L=C_{L,maxSAR},t_p, M=M_{MO}, H=H_{opt}$','Location','southeast','interpreter','latex');
 
 figure(2);
 grid on;
 plot(zeta,R1,'b', 'Linewidth', 2);
 xlabel('Block Fuel Ratio, \zeta, (-)');
 ylabel('Aircraft Range,(km)');
 legend('$C_L=C_{L,maxSAR},t_p, M=M_{MO}, H=H_{opt}$','Location','southeast','interpreter','latex');
 
 figure(3);
 grid on;
 plot(R1,FE,'b', 'Linewidth', 2);
 xlabel('Aircraft Range,(km)');
 ylabel('Fuel Economy, (L Fuel per 100 km per 100 kg of payload)');
 legend('$C_L=C_{L,maxSAR},t_p, M=M_{MO}, H=H_{opt}$','Location','southeast','interpreter','latex');
 
 
 


% for k=1:length(H)
% 
%     Cl=w./((1/2)*roh_(k)*(a_(k)*M(k))^2*S);
%     Cd=CD_0+K.*Cl.^2;
%     TR(k,:)=0.5*roh_(k).*(a_(k)*M(k))^2*S.*Cd;
%     SAR(k,:)=(a_(k)*M(k))./(TR(k,:).*ctp);
%     
%     
%    
% end

% figure;
% grid on;
% hold on;
% plot(v,SAR(1,:),'b','LineWidth',2);
% plot(v,SAR(2,:),'g','LineWidth',2);
% plot(v,SAR(3,:),'r','LineWidth',2);
% xlabel('Trim Velocity, V_t_r_i_m, (m/s)');
% ylabel('Specific Air Range, SAR_t_p, (m/kg)');
% legend('10000 ft','36000 ft','Location','northwest');
% 
% 
% figure;
% grid on;
% hold on;
% plot(v,TR(1,:)./1000,'b','LineWidth',2);
% plot(v,TR(2,:)./1000,'g','LineWidth',2);
% plot(v,TR(3,:)./1000,'r','LineWidth',2);
% set(gca,'xlim',[0 300],'ylim',[0 150])
% xlabel('Trim Velocity, V_t_r_i_m, (m/s)');
% ylabel('Thrust Required, TR, or Drag, D, (kN)');
% legend('10000 ft','36000 ft','Location','northwest');



