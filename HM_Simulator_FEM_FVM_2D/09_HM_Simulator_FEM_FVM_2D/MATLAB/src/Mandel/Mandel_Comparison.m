function [Mandel] = Mandel_Comparison(a,b,cf,E,F,k,Ks,mu,phi,t,vd)

%Compressible Mandel from Castelletto 2015

%% - Inputs
%beta points
n = 100;                        %[-]

%% - Basic Calculations
%x-location
x = linspace(0,a,100);
y = linspace(0,b,100);
%Mobility
kmob = k/mu;

%% - Poroelastic calculations - with drained parameters
%Shear Modulus (same drained and undrained)
G = E/(2*(1+vd));
%Unixial Drained bulk modulus
Kv = 2*G*(1-vd)/(1-2*vd);
%Bulk Modulus
K = 2*G*(vd + 1)/(3*(1-2*vd));
%Biot coefficient
alpha = 1 - K/Ks;
%Fluid Bulk modulus
Kf = 1/cf;
%Biot's Modulus
M = ((alpha - phi)/Ks + phi/Kf)^-1;
%Skempton Coefficient
B = alpha*Kf/(alpha - phi*(1-alpha)*Kf + phi*K);
B = alpha*M/(K + alpha^2*M);
%Undrained Poisson's Ratio
vu = (B*alpha*(1-2*vd) + 3*vd)/(3 - B*alpha*(1 - 2*vd));
%Fluid diffusivity coefficient
cv = kmob*M*Kv/(Kv + alpha^2*M);
%Auxillary elastic constant
nu = (1-vd)/(vu-vd);


%% - Beta
beta = zeros(n,1) + pi()/2;
for i = 1:n
    Error = 1;
    while Error > 1e-8
        h = i-1;
        betaold = beta(i);
        beta(i) = atan(beta(i)*2*nu) + (h)*pi();
        Error = abs((betaold - beta(i))/beta(i));
    end
end

%% - Solve for pore pressure
Pd = zeros(length(x),length(t));
for i = 1:length(t)
    for j = 1:length(x)
        sum = 0;
        for kk = 1:n
            sum = sum + sin(beta(kk))*...
                (cos(beta(kk)*x(j)/a) - cos(beta(kk)))*...
                exp(-beta(kk)^2*cv*t(i)/a^2)/...
                (beta(kk) - sin(beta(kk))*cos(beta(kk)));
        end
        Pd(j,i) = 2*sum;
    end
end

%% - Solve for horizontal displacement
ux = zeros(length(x),length(t));
for i = 1:length(t)
    for j = 1:length(x)
        sum = 0;
        sum1 = 0;
        for kk = 1:n
            sum = sum + sin(beta(kk))*cos(beta(kk))/(beta(kk) - sin(beta(kk))*cos(beta(kk)))*exp(-beta(kk)^2*cv*t(i)/a^2);
            sum1 = sum1 + sin(beta(kk)*x(j)/a)*cos(beta(kk))/(beta(kk) - sin(beta(kk))*cos(beta(kk)))*exp(-beta(kk)^2*cv*t(i)/a^2);
        end
        ux(j,i) = (F*vd/(2*G*a) - F*vu/(G*a)*sum)*x(j) + F/G*sum1;
    end
end


%% - Solve for vertical displacement
uy = zeros(length(x),length(t));
for i = 1:length(t)
    for j = 1:length(x)
        sum = 0;
        for kk = 1:n
            sum = sum + sin(beta(kk))*cos(beta(kk))*exp(-beta(kk)^2*cv*t(i)/(a^2))/(beta(kk) - sin(beta(kk))*cos(beta(kk)));
        end
        uy(j,i) = (-F*(1-vd)/(2*G*a) + F*(1-vu)*sum/(G*a))*y(j);
    end
end

%% - Outputs
Mandel.x = x;
Mandel.y = y;
Mandel.Pd = Pd;
Mandel.t = t;
Mandel.td = cv*t/(a^2);
Mandel.p0 = -B*(1+vu)*F/(3*a);
Mandel.tnorm = cv/(a^2);

Mandel.ux0 = F*vu/(2*G);
Mandel.ux = ux;
Mandel.uy0 = -F*(1-vu)/(2*G);
Mandel.uy = uy;


end

