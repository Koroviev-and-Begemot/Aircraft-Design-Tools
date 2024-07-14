clear
close all

altitude = [ 0 10000 36000];
s = 845;
m = 575000;
g = 9.81;
w = m*g;
v = 1:300;
n = altitude.*0.3048;
TSFC = 1.5863*10^-5;
CD0 = 0.0239;
k = 0.0247;

TR  = zeros(length(altitude),length(v));
SAR = zeros(length(altitude),length(v));

for x = 1:length(n) 
    
[~,~,rho,~] = GetISA(n(x),0,0);

CL = w./(0.5*rho.*v.^2*s);
CD = CD0 + k.*CL.^2;

TR(x,:) = 0.5*rho.*v.^2*s.*CD;

SAR(x,:) = v./(TR(x,:).*TSFC);

end

        
        figure
        hold on
        plot(v,TR(1,:)./1000,'b-','LineWidth',2)
        plot(v,TR(2,:)./1000,'g-','LineWidth',2)
        plot(v,TR(3,:)./1000,'r-','LineWidth',2)
        grid on
        legend('0 ft','10000 ft','36000 ft')
        title('Thrust Required vs Trim Velocity ')
        xlabel('Trim Velocity m/s')
        ylabel('Thrust Required or Drag kN')
        set(gca,'xlim',[0 300],'ylim',[250 400])
        hold off
       
        
        figure
         hold on
        plot(v,SAR(1,:),'b-','LineWidth',2)
        plot(v,SAR(2,:),'g-','LineWidth',2)
        plot(v,SAR(3,:),'r-','LineWidth',2)
        grid on
        legend('0 ft','10000 ft','36000 ft')
        title('Specific Air Range vs Trim Velocity ')
        xlabel('Trim Velocity m/s')
        ylabel('Specific Air Range m/kg')
        set(gca,'xlim',[0 300],'ylim',[0 100])
        



