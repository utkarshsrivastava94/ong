%Utkarsh Srivastava
imax = 20;
dx = 20; %ft
dy = 30; %ft
dz = 30; %ft
dt = 0.005; %day
mu_o = 0.8; %cP
mu_w = 0.5; %cP
km = 10; %mD
phim = 0.2;
phif = 0.005;
ut = 14; %ft/day
gammaw = 0.45; %psi/ft
gammao = 0.35; %psi/ft
swrf = 0.05;
sorf = 0.05;
swrm = 0.35;
sorm = 0.30;
alpham1 = -0.55;
swmx = 0.55;
somx = 1 - swmx;
pcmsf = 2; %psia
pcfsf = 0.1; %psia
epsilon = 0.0001;
krwfs = 0.9;
krofs = 0.9;
krwms = 0.05;
kroms = 0.7;
nwf = 2;
nwm = 2;
nof = 3;
nom = 3;
t = 300/dt;
swf = zeros(imax+1, t);
sof = zeros(imax+1, t);
swm = zeros(imax+1, t);
som = zeros(imax+1, t);
fwf = zeros(imax+1, t);
fwm = zeros(imax+1, t);
k_rwf = zeros(imax+1, t);
k_rof = zeros(imax+1, t);
k_rwm = zeros(imax+1, t);
k_rom = zeros(imax+1, t);
lambdawf = zeros(imax+1, t);
lambdaof = zeros(imax+1, t);
lambdaom = zeros(imax+1, t);
hwf = zeros(imax+1, t);
hwm = zeros(imax+1, t);
pcwom = zeros(imax+1, t);
pcwof = zeros(imax+1, t);
tau = zeros(imax+1, t);
SWF = zeros(imax+1,1);
SWM =zeros(imax+1 , 1);
KROM = zeros(imax+1, 1);
KROF = zeros(imax+1, 1);
KRWM = zeros(imax+1, 1);
KRWF = zeros(imax+1, 1);
Lambdawf = zeros(imax+1, 1);
Lambdaof = zeros(imax+1, 1);
FWF = zeros(imax+1, 1);
Pcwom = zeros(imax+1,1);
Pcwof = zeros(imax+1, 1);
lx = 3; ly = 3; lz = 3; %ft
sigma = 4*((1/lx^2)+(1/ly^2)+(1/lz^2)); %psi^-2
sigmaz = 4/lz^2; %psi^-2
alpham2 = -((swmx - swrm)/(1-swmx-sorm))*alpham1;
%part a
beta = ((1-swrf)-swrf)/imax;
delta = ((1-sorm)-swrm)/imax;
for i = 1:imax+1
    if i == 1
        SWF(i,1) = swrf;
        SWM(i,1) = swrm;
    else
        SWF(i,1) = SWF(i-1,1) + beta;
        SWM(i,1) = SWM(i-1,1) + delta;
    end
end
for i = 1:imax+1
    KRWF(i,1) = krwfs*((SWF(i,1)-swrf)/(1-swrf-sorf))^nwf;
    KROF(i,1) = krofs*((1-SWF(i,1)-sorf)/(1-swrf-sorf))^nof;
    KRWM(i,1) = krwms*((SWM(i,1)-swrm)/(1-swrm-sorm))^nwm;
    KROM(i,1) = kroms*((1-SWM(i,1)-sorm)/(1-swrm-sorm))^nom;
    Lambdawf(i,1) = KRWF(i,1)/mu_w;
    Lambdaof(i,1) = KROF(i,1)/mu_o;
    FWF(i,1) = Lambdawf(i,1)/(Lambdawf(i,1)+Lambdaof(i,1));
    if (SWM(i,1) < swmx) && (SWM(i,1) > swrm)
        Pcwom(i,1) = pcmsf*alpham2*log((1-somx-swrm)/(SWM(i,1)-swrm+epsilon));
    elseif (SWM(i,1) < (1-sorm)) && (SWM(i,1) > swmx)
        Pcwom(i,1) = alpham1*log((1-swmx-sorm)/(1-SWM(i,1)-sorm+epsilon));
    elseif SWM(i,1) <= swrm
        Pcwom(i,1) = pcmsf*alpham2*log((1-somx-swrm)/(swrm-swrm+epsilon));
    end
    if SWF(i,1) <= 1-sorf
        Pcwof(i,1) = pcfsf*reallog((1-sorf-swrf)/(SWF(i,1)-swrf + epsilon));
    end
end
%part b
for j = 1:t+1
    swf(1,j) = 1; %initial node of injection water saturation
    swm(1,j) = 1;
    sof(1,j) = 1 - swf(1,j);
    som(1,j) = 1 - swm(1,j);
    fwf(1,j) = 1;
end
for i = 2:imax+1
    swf(i,1) = swrf;
    sof(i,1) = 1 - swf(i,1);
    swm(i,1) = swrm;
    som (i,1) = 1 - swm(i,1);
end
for j =1:t
    for i = 1:imax+1
        %Fracture rel perm
        if swf(i,j) < swrf
            k_rwf(i,j) = 0;
            k_rof(i,j) = krofs;
        else
            k_rwf(i,j) = krwfs*((swf(i,j)-swrf)/(1-swrf-sorf))^nwf;
            k_rof(i,j) = krofs*((sof(i,j)-sorf)/(1-swrf-sorf))^nof;
        end
        if swf(i,j) >= (1-sorf)
            k_rwf(i,j) = krwfs;
            k_rof(i,j) = 0;
        else
            k_rwf(i,j) = krwfs*((swf(i,j)-swrf)/(1-swrf-sorf))^nwf;
            k_rof(i,j) = krofs*((sof(i,j)-sorf)/(1-swrf-sorf))^nof;
        end
        %Matrix rel perm
        if swm(i,j) < swrm
            k_rwm(i,j) = 0;
            k_rom(i,j) = kroms;
        else
            k_rwm(i,j) = krwms*((swm(i,j)-swrm)/(1-swrm-sorm))^nwm;
            k_rom(i,j) = kroms*((som(i,j)-sorm)/(1-swrm-sorm))^nom;
        end
        if swm(i,j) >= (1-sorm)
            k_rwm(i,j) = krwms;
            k_rom(i,j) = 0;
        else
            k_rwm(i,j) = krwms*((swm(i,j)-swrm)/(1-swrm-sorm))^nwm;
            k_rom(i,j) = kroms*((som(i,j)-sorm)/(1-swrm-sorm))^nom;
        end
        lambdawf(i,j) = k_rwf(i,j)/mu_w;
        lambdaof(i,j) = k_rof(i,j)/mu_o;
        lambdaom(i,j) = k_rom(i,j)/mu_o;
        if i == 1
            fwf(i,j) = 1;
        else
            fwf(i,j) = lambdawf(i,j)/(lambdawf(i,j)+lambdaof(i,j));
        end
        %capillary pressure calculation
        if (swm(i,j) < swmx) && (swm(i,j) > swrm)
            pcwom(i,j) = pcmsf*alpham2*log((1-somx-swrm)/(swm(i,j)-swrm+epsilon));
        elseif (swm(i,j) < (1-sorm)) && (swm(i,j) > swmx)
            pcwom(i,j) = alpham1*log((1-swmx-sorm)/(1-swm(i,j)-sorm+epsilon));
        elseif swm(i,j) <= swrm
            pcwom(i,j) = pcmsf*alpham2*log((1-somx-swrm)/(swrm-swrm+epsilon));
        end
        if swf(i,j) < 1-sorf
            pcwof(i,j) = pcfsf*reallog((1-sorf-swrf)/(swf(i,j)-swrf + epsilon));
        end
        hwf(i,j) = lz*((swf(i,j)-swrf)/(1-sorf-swrf));
        hwm(i,j) = lz*((swm(i,j)-swrm)/(1-sorm-swrm));
        %transfer function
        tau(i,j) = 0.006328*sigma*km*(lambdawf(i,j)/(lambdawf(i,j)+lambdaom(i,j)))*(lambdaom(i,j)*((pcwom(i,j)-pcwof(i,j))+(sigmaz/sigma)*(gammaw-gammao)*(hwf(i,j)-hwm(i,j))));
    end
    %saturation calculation
    for i = 2:imax+1
        swf(i,j+1) = swf(i,j) + dt*(-ut*((fwf(i,j)-fwf(i-1,j))/dx)-(tau(i,j)/phif));
        swm(i,j+1) = swm(i,j) + dt*(tau(i,j)/phim);
        sof(i,j+1) = 1 - swf(i,j+1);
        som(i,j+1) = 1 - swm(i,j+1);
    end
end
%part c
for i = 1:imax
    x(i,1) = i*dx;
end
for i = 2:imax+1
    Swf_300(i-1,1) = swf(i,t);
    Swf_200(i-1,1) = swf(i,t-2000);
    Swf_100(i-1,1) = swf(i,t-4000);
    Swm_300(i-1,1) = swm(i,t);
    Swm_200(i-1,1) = swm(i,t-2000);
    Swm_100(i-1,1) = swm(i,t-4000);
end
%part d
qt = ut*dz*dy*phif*0.17811; %cu ft/day to bbl/day
qw = zeros(t,1); 
qo = zeros(t,1); 
cumm = zeros(t,1); %bbl/day
WOR = zeros (t,1);
time = zeros (t,1); %days

for i = 1:t
    time(i,1) = i*0.05; %days
    qw(i,1) = fwf(21,i)*qt; %bbl/day
    qo(i,1) = (1-fwf(21,i))*qt; %bbl/day
    WOR(i,1) = qw(i,1)/qo(i,1);
    if  i == 1
        cumm(i,1) = qo(i,1);
    else
        cumm(i,1) = qo(i,1) + cumm(i-1,1);
    end
end
% 
% figure(1)
% plot(x, Swf_300);
% axis ([20 400 0 1]);
% xlabel('Distance (ft)');
% ylabel('Fracture Water saturation (S_{wf})');
% hold on;
% plot (x, Swf_100);
% hold on;
% x1=0:20:400;
% y1(1:21) = swrf;
% y2(1:21)= (1-sorf);
% plot (x1,y1,'--');
% hold on;
% plot (x1,y2,'--');
% legend ('t = 300 days','t = 100 days');
% hold off;
% 
% figure(2)
% plot(x, Swm_300);
% axis ([20 400 0 1]);
% xlabel('Distance (ft)');
% ylabel('Matrix Water saturation (S_{wm})');
% hold on;
% plot (x, Swm_100);
% hold on;
% x1=0:20:400;
% y1(1:21) = swrm;
% y2(1:21)= (1-sorm);
% plot (x1,y1,'--');
% hold on;
% plot (x1,y2,'--');
% legend ('t = 300 days','t = 100 days');
% hold off;
% 
% figure(3)
% plot(time, qo);
% axis([0 300 0 15]);
% xlabel('Time (days)');
% ylabel('Flow rate (bbls/day)');
% hold on;
% plot(time, qw);
% hold on;
% legend('Oil', 'Water');
% hold off;
% 
% figure(4)
% plot(time, WOR);
% axis([0 300 0 43]);
% xlabel('Time (days)');
% ylabel('WOR (bbl water/bbl oil)');
% 
% figure(5)
% plot(time, cumm);
% axis([0 300 0 60000]);
% xlabel('Time (days)');
% ylabel('Oil production (bbls)');
% 
% figure(6)
% plot(SWF, KRWF);
% axis([0 1 0 1]);
% xlabel('Water saturation');
% ylabel('Relative permeability');
% hold on;
% plot(SWF, KROF);
% legend('k_{rw,f}', 'k_{ro,f}');
% hold off;
% 
% figure(7)
% plot(SWM, KRWM);
% axis([0 1 0 1]);
% xlabel('Water saturation');
% ylabel('Relative permeability');
% hold on;
% plot(SWM, KROM);
% legend('k_{rw,m}', 'k_{ro,m}');
% hold off;
% 
% figure(8)
% plot(SWF, FWF);
% axis([0 1 0 1]);
% xlabel('Water saturation');
% ylabel('Fractional flow');
% 
% figure(9)
% plot(SWM, Pcwom);
% axis([0 1 -5 12]);
% xlabel('Water saturation');
% ylabel('Capillary pressure');
% hold on;
% plot(SWF, Pcwof);
% legend('Matrix capillary pressure', 'Fracture capillary pressure');
% hold off;









