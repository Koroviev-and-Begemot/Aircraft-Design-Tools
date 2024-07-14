function [h,a,t] = GetISAinverse(p,dT,dL)

p0 = 101325;
t0 = 288.15+dT;
L = -0.0065+dL;
g = 9.81;
R = 287;

h = t0/L.*((p/p0).^((-L.*R)./g)-1);

t = t0+L.*h;
a = sqrt(1.4*287.*t);

%[~,t,~,a] = ISA(h,dT,dL);

end

