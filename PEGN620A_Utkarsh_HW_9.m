% clear all
% close all

%%%%%%Utkarsh Srivastava
imax = 20;
dx = 20; %ft
dy = 30; %ft
dz = 30; %ft
dt = 0.005; %day
mu_o = 0.8; %cP
mu_w = 0.5; %cP
km = 10; %mD
kfeff = 100; %mD
phim = 0.2;
phif = 0.005;
gammaw = 0.45; %psi/ft
gammao = 0.35; %psi/ft
qw = 350; %ft3/day
pout = 1600; %psi
swrf = 0.05;
sorf = 0.05;
swrm = 0.35;
sorm = 0.30;
alpham1 = -0.55;
swmx = 0.55;
somx = 0.45;
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
chi = 0.006328; %conversion
time = 300;
t = time/dt;
swf = zeros(imax, t);
sof = zeros(imax, t);
swm = zeros(imax, t);
som = zeros(imax, t);
fwf = zeros(imax, t);
fwm = zeros(imax, t);
k_rwf = zeros(imax, t);
k_rof = zeros(imax, t);
k_rwm = zeros(imax, t);
k_rom = zeros(imax, t);
lambdawf = zeros(imax, t);
lambdaof = zeros(imax, t);
lambdaom = zeros(imax, t);
hwf = zeros(imax, t);
hwm = zeros(imax, t);
pcwom = zeros(imax, t);
pcwof = zeros(imax, t);
p = zeros(imax,t);
e_11 = zeros(imax, t);
e_12 = zeros(imax, t);
e_13 = zeros(imax, t);
e_21 = zeros(imax, t);
e_22 = zeros(imax, t);
e_23 = zeros(imax, t);
e_31 = zeros(imax, t);
e_33 = zeros(imax, t);
e_34 = zeros(imax, t);
e_41 = zeros(imax, t);
e_43 = zeros(imax, t);
e_44 = zeros(imax, t);
d_11 = zeros(imax, t);
d_21 = zeros(imax, t);
f_11 = zeros(imax, t);
f_21 = zeros(imax, t);
pom = zeros(imax, t);
pof = zeros(imax, t);
pwm = zeros(imax, t);
pwf = zeros(imax, t);
Rpof = zeros(imax, t);
Rpom = zeros(imax, t);
Rswf = zeros(imax, t);
Rswm = zeros(imax, t);
A = zeros(4*imax, 4*imax);
R = zeros(4*imax,1);
p = zeros(4*imax,t);
pe = 1600; %psi
lx = 3; ly = 3; lz = 3; %ft
sigma = 4*((1/lx^2)+(1/ly^2)+(1/lz^2)); %psi^-2
sigmaz = 4/lz^2; %psi^-2
alpham2 = -((swmx - swrm)/(1-swmx-sorm))*alpham1;
%%%%part a
beta = ((1-swrf)-swrf)/imax;
delta = ((1-sorm)-swrm)/imax;
WIw = (chi*kfeff*dy*dz)/(dx/2);
WIo = (chi*kfeff*dy*dz)/(dx/2);
%%%%%initialization

for i = 1:imax
    swf(i,1) = swrf;
    
    swm(i,1) = swrm;
    
end
for i = 1:imax
    pof(i,1) = pe;
    pom(i,1) = pe;
    
end
for j =1:t
    for i = 1:imax
        %%%%%        Fracture rel perm
        sof(i,j) = 1 - swf(i,j);
        som (i,j) = 1 - swm(i,j);
        k_rwf(i,j) = krwfs*((swf(i,j)-swrf)/(1-swrf-sorf)).^nwf;
        k_rof(i,j) = krofs*((sof(i,j)-sorf)/(1-swrf-sorf)).^nof;
        
        
        %%%%%        Matrix rel perm
        
        k_rwm(i,j) = krwms*((swm(i,j)-swrm)/(1-swrm-sorm)).^nwm;
        k_rom(i,j) = kroms*((som(i,j)-sorm)/(1-swrm-sorm)).^nom;
    end
    for i=1:imax
        %%%%%%%%%%    mobility calculation
        lambdawf(i,j) = k_rwf(i,j)/mu_w;
        lambdaof(i,j) = k_rof(i,j)/mu_o;
        lambdaom(i,j) = k_rom(i,j)/mu_o;
        lambdawm(i,j) = k_rwm(i,j)/mu_w;
    end
    for i=1:imax
        %%%%%%%%%%     capillary pressure calculation
        if (swm(i,j) < swmx) && (swm(i,j) >= swrm)
            pcwom(i,j) = pcmsf*alpham2*log((1-somx-swrm)/(swm(i,j)-swrm+epsilon));
        elseif (swm(i,j) <= (1-sorm)) && (swm(i,j) >= swmx)
            pcwom(i,j) = alpham1*log((1-swmx-sorm)/(1-swm(i,j)-sorm+epsilon));
        end
        
        pcwof(i,j) = pcfsf*log((1-sorf-swrf)/(swf(i,j)-swrf + epsilon));
    end
    for i=1:imax
        pwf(i,j) = pof(i,j) - pcwof(i,j);
        pwm(i,j) = pom(i,j) - pcwom(i,j);
    end
    for i=1:imax
        hwf(i,j) = lz*((swf(i,j)-swrf)/(1-sorf-swrf));
        hwm(i,j) = lz*((swm(i,j)-swrm)/(1-sorm-swrm));
        
    end
    
    %coefficient input
    
    for i=1:imax
        if i > 1
            if  pof(i,j) > pof(i-1,j)
                d_11(i,j) = (dy * dz * chi * kfeff * lambdaof(i,j))/dx;
            else
                d_11(i,j) = (dy * dz * chi * kfeff * lambdaof(i-1,j))/dx;
            end
        end
        if i > 1
            if pwf(i,j) > pwf(i-1,j)
                d_21(i,j) = (dy * dz * chi * kfeff * lambdawf(i,j))/dx;
            else
                d_21(i,j) = (dy * dz * chi * kfeff * lambdawf(i-1,j))/dx;
            end
        end
    end
    for i=1:imax
        if i < imax
            if pof(i,j) >= pof(i+1,j)
                f_11(i,j) = (dy * dz * chi * kfeff * lambdaof(i,j))/dx;
            else
                f_11(i,j) = (dy * dz * chi * kfeff * lambdaof(i+1,j))/dx;
            end
        end
        if i < imax
            if pwf(i,j) >= pwf(i+1,j)
                f_21(i,j) =  (dy * dz * chi * kfeff * lambdawf(i,j))/dx;
            else
                f_21(i,j) =  (dy * dz * chi * kfeff * lambdawf(i+1,j))/dx;
            end
        end
    end
    
    
    
    %
    for i=1:imax
        hyd_w(i,1)=(pof(i,j) - pom(i,j))+(sigmaz/sigma)*gammaw*(hwf(i,j)-hwm(i,j))-(pcwof(i,j)-pcwom(i,j));
        hyd_o(i,1)=(pof(i,j) - pom(i,j))+(sigmaz/sigma)*gammao*(hwf(i,j)-hwm(i,j));
    end
    for i=1:imax
        if hyd_w(i,1) >= 0
            e_23(i,j) = (dx * dy * dz * sigma * chi * km * lambdawf(i,j)); %MFR
            e_41(i,j) = chi * sigma * km * lambdawf(i,j);
        else
            e_23(i,j) = (dx * dy * dz * sigma * chi * km * lambdawm(i,j)); %MFR
            e_41(i,j) = chi * sigma * km * lambdawm(i,j);
        end
        if hyd_o(i,1) >= 0
            e_31(i,j) = chi * sigma * km * lambdaof(i,j);
            e_13(i,j) = dx * dy * dz *chi * sigma * km * lambdaof(i,j); %MFL
        else
            e_31(i,j) = chi * sigma * km * lambdaom(i,j); %MFL
            e_13(i,j) = dx * dy * dz *chi * sigma * km * lambdaom(i,j);
        end
        e_33(i,j) = -e_31(i,j);
        e_43(i,j) = -e_41(i,j);
        e_44(i,j) = - phim/dt;
        e_34(i,j) = phim/dt;
    end
    for i = 1:imax
        %LHS
        if i == imax
            %
            e_11(i,j)= -d_11(i,j)-e_13(i,j)-WIo * lambdaof(i,j);
            %
        elseif  i == 1
            %
            e_11(i,j)= -f_11(i,j)-e_13(i,j);
        else
            %
            e_11(i,j)= -d_11(i,j)-f_11(i,j)-e_13(i,j);
        end
        e_12(i,j) = (dx * dy * dz * phif)/dt;
        
    end
    for i=1:imax
        %
        if i == imax
            e_21(i,j) = -d_21(i,j)-e_23(i,j)- (WIw * lambdawf(i,j));
        elseif i==1
            e_21(i,j) = -f_21(i,j)-e_23(i,j);
        else
            e_21(i,j) = -d_21(i,j)-f_21(i,j)-e_23(i,j);
        end
        e_22(i,j) = - e_12(i,j);
    end
    
    %RHS
    for i=1:imax
        if i == imax
            Rpof(i,j) =  (e_13(i,j) * (sigmaz/sigma) * gammao * (hwf(i,j) - hwm(i,j)))  + (dx * dy * dz * swf(i,j) * (phif/dt)) - (WIo * lambdaof(i,j)* pout);
        else
            Rpof(i,j) =  (e_13(i,j) * (sigmaz/sigma) * gammao * (hwf(i,j) - hwm(i,j)))  + (dx * dy * dz * swf(i,j) * (phif/dt));
        end
        Rpom(i,j) = (phim * swm(i,j)/dt) - (e_31(i,j) * (sigmaz/sigma) * gammao * (hwf(i,j) - hwm(i,j)));
        Rswm(i,j) = -(phim * swm(i,j)/dt) - (sigma * chi * km * lambdawf(i,j) * ( (sigmaz/sigma) * gammaw * (hwf(i,j) - hwm(i,j)) - (pcwof(i,j) - pcwom(i,j))));
        if i == 1
            Rswf(i,j) = dy * dz * chi *((kfeff * lambdawf(i,j))*((pcwof(i+1,j) - pcwof(i,j))/dx) - (kfeff * 0)*(pcwof(i,j) - 0)/dx)+(dx * dy * dz * sigma * km * chi * lambdawf(i,j) * ((sigmaz/sigma)*gammaw*(hwf(i,j) - hwm(i,j)) - (pcwof(i,j) - pcwom(i,j)))) - qw - (dx * dy * dz * phif * swf(i,j)/dt);
        elseif i == imax
            Rswf(i,j) = dy * dz * chi *((kfeff * 0)*((0 - pcwof(i,j))/dx) - (kfeff * lambdawf(i-1,j))*(pcwof(i,j) - pcwof(i-1,j))/dx)+(dx * dy * dz * sigma * km * chi * lambdawf(i,j) * ((sigmaz/sigma)*gammaw*(hwf(i,j) - hwm(i,j)) - (pcwof(i,j) - pcwom(i,j)))) - ( WIw*lambdawf(i,j)*pcwof(i,j)+WIw*lambdawf(i,j)*pout ) -((dx*dy*dz)*(phif*swf(i,j)/dt));
        else
            Rswf(i,j) = dy * dz * chi *((kfeff * lambdawf(i,j))*((pcwof(i+1,j) - pcwof(i,j))/dx) - (kfeff * lambdawf(i-1,j))*(pcwof(i,j) - pcwof(i-1,j))/dx)+(dx * dy * dz * sigma * km * chi * lambdawf(i,j) * ((sigmaz/sigma)*gammaw*(hwf(i,j) - hwm(i,j)) - (pcwof(i,j) - pcwom(i,j)))) - (dx * dy * dz * phif * swf(i,j)/dt);
        end
    end
    %Constructing the A matrix
    for i = 1:imax
        b = (i-1)*4+1;
        A(b+1,b) = e_11(i,j);
        A(b+1,b+1) = e_12(i,j);
        A(b,b) = e_21(i,j);
        A(b,b+1) = e_22(i,j);
        A(b+1,b+2) = e_13(i,j);
        A(b,b+2) = e_23(i,j);
        A(b+3,b) = e_31(i,j);
        A(b+2,b) = e_41(i,j);
        A(b+3,b+2) = e_33(i,j);
        A(b+3,b+3) = e_34(i,j);
        A(b+2,b+2) = e_43(i,j);
        A(b+2,b+3) = e_44(i,j);
    end
    for i = 1:imax-1
        b = (i-1)*4 + 1;
        A(b+1,b+4) = f_11(i,j);
        A(b,b+4) = f_21(i,j);
    end
    for i = 2:imax
        b = (i-1)*4 +1;
        A(b+1,b-4) = d_11(i,j);
        A(b,b-4) = d_21(i,j);
    end
    for i = 1:imax
        b = (i-1)*4 + 1;
        R(b+1,j) = Rpof(i,j);
        R(b,j) = Rswf(i,j);
        R(b+3,j) = Rpom(i,j);
        R(b+2,j) = Rswm(i,j);
    end
    
    %calcualting the x matrix
    p(:,1) = A\R(:,j);
    %Seperating the values from the p matrix
    for i =1:imax
        b = (i-1)*4+1;
        pof(i,j+1) = p(b,1);
        swf(i,j+1) = p(b+1,1);
        pom(i,j+1) = p(b+2,1);
        swm(i,j+1) = p(b+3,1);
        sof(i,j+1) = 1 - swf(i,1);
        som(i,j+1) = 1 - swm(i,1);
        %
    end
    
    %Saturation values limit
    for i= 1:imax
        if swf(i,j+1)<swrf
            swf(i,j+1) = swrf;
        elseif swf(i,j+1)> (1-sorf)
            swf(i,j+1) = 1-sorf;
        end
        if swm(i,j+1)<swrm
            swm(i,j+1) = swrm;
        elseif swm(i,j+1) > (1-sorm)
            swm(i,j+1) = 1-sorm;
        end
    end
end
for j = 1:t
    Swm1(j,1) = swm(1,j);
    Som1(j,1) = som(1,j);
    Swm20(j,1) = swm(20,j);
    Som20(j,1) = som(20,j);
    Pof1(j,1) = pof(1,j);
    Pof20(j,1) = pof(20,j);
    time(j,1) = j * 0.005;
end

figure(1)
plot(time, Swm1);
axis([0 300 0 1]);
xlabel('Number of days');
ylabel('Saturation of water over time in matrix');
hold on;
plot(time, Som1);
hold on;
y1(1:t) = swrm;
plot(time, y1, '--');
legend('Swm','Som');
hold off;

figure(2)
plot(time, Swm20);
axis([0 300 0 1]);
xlabel('Number of days');
ylabel('Saturation of water over time in matrix');
hold on;
plot(time, Som20);
hold on;
y1(1:t) = swrm;
plot(time, y1, '--');
legend('Swm','Som');
hold off;

figure(3)
plot(time, Pof1)
xlabel('Number of days');
ylabel('Pressue of oil over time in matrix (psi)');

figure(4)
plot(time, Pof20)
xlabel('Number of days');
ylabel('Pressue of oil over time in matrix (psi)');




