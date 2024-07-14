
function [p,t,rho,a] = GetISA(Altitude,DeltaT,DeltaL)


%This function takes the inputs Altitide, Delta Temp. and Delta Lapse rate and 
% returns 4 variables: pressure, temprature, density and the speed of sound
% in that order for the given inputs.

    h = Altitude; 
    t0 = 288.15+DeltaT;
    L = -0.0065+DeltaL;
    p0 = 101325;
    p11 = p0*(1+(L/t0)*11000)^(-9.81/(287*L));

   if h<11000
        
        p = p0.*(1+(L/t0).*h).^(-9.81/(287*L));
        t = t0+L.*h;
        
    else
        
        p = p11.*exp((-1.576*10^(-4).*(h-11000)));
        t = t0+L*11000;
    end
    
    rho = p./(287.*t);
    a = sqrt(1.4*287.*t);
    
end

