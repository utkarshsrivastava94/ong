% l_log function credits Ilkay Iker
clear all
close all
xe = 250; %ft
ye = 250; %ft
imax = 12;
jmax = 10;
deltax = xe/imax; %ft
deltax_hf = 0.0164; %ft half hydraulic fracture width
deltay = ye/jmax; %ft
delz = 40; %ft
delt = 0.01; %day
h = 40; %ft
k_m = 0.02; %mD permeabilty matrix
k_hf = 10000; %mD peremabilty hydraulic fracture
phi_m = 0.05; %porosity matrix
phi_hf = 0.26; %porosity hydraulic fracture
mu = 0.3; %cP
ct = 0.0001; %psi^-1
chi = 0.006328; %conversion constant
t_step = 3000;
p_well = 3000; %psi
p_i = 7000; %psi
q_calc = zeros(imax*jmax,1); %rcf/day
p = zeros(imax*jmax,t_step);%psia
A = zeros(imax*jmax,imax*jmax);
delxa = zeros(imax,jmax);
C = zeros(imax*jmax,1);
G = zeros(imax*jmax,1);
D = zeros(imax*jmax,1);
F = zeros(imax*jmax,1);
E = zeros(imax*jmax,1);
R = zeros(imax*jmax,t_step);
k = zeros(imax*jmax,1);
k_hy = zeros(imax*jmax,1);
k_hx = zeros(imax*jmax,1);
phi = zeros(imax*jmax,1);
V = zeros(imax*jmax,1);
%Values for all node
deltx = l_log(deltax_hf,imax,xe);
for j = 1:jmax
    for i = 1:imax
        delxa(i,j) = deltx(i,1);
    end
end
delx = reshape(delxa,[],1);
for i = 1: imax*jmax
    k(i,1) = k_m;
%   delx(i,1) = deltax;
    phi(i,1) = phi_m;
    V(i,1) = delx(i,1)*deltay*delz;
end
%Values only for fractures
for i = 1:(jmax)
    k((((i-1)*imax)+1),1) = k_hf;
%   delx((((i-1)*imax)+1),1) = deltax_hf;
    phi((((i-1)*imax)+1),1) = phi_hf;
    V((i-1)*imax+1,1) = delx((i-1)*imax+1,1)*deltay*delz;
end
%filling of matrix elements
for i = (imax+1):imax*jmax
    k_hy(i,1) = H_avg(k(i),k(i-imax),deltay,deltay);
    C(i,1) = chi*(k_hy(i,1)/mu)*delx(i,1)*delz/deltay;
end
for i = 1:((imax*jmax)-imax)
    k_hy(i,1) = H_avg(k(i),k(i+imax),deltay,deltay);
    G(i,1) = chi*(k_hy(i,1)/mu)*delx(i,1)*delz/deltay;
end
for i = 2:imax*jmax
    k_hx(i,1) = H_avg(k(i),k(i-1),delx(i),delx(i-1));
    D(i,1) = chi*(k_hx(i,1)/mu)*deltay*delz/((delx(i-1,1)+delx(i,1))/2);
end
for i = 1:jmax
    D((i-1)*imax+1,1) = 0;
end
for i = 1:(imax*jmax-1)
    k_hx(i,1) = H_avg(k(i),k(i+1),delx(i),delx(i+1));
    F(i,1) = chi*(k_hx(i,1)/mu)*deltay*delz/((delx(i,1)+delx(i+1,1))/2);
end
for i = 1:jmax
    F(i*imax,1) = 0;
end   
for i = 1:imax*jmax
    if i == 1
        wi = chi*(k(i,1)/mu)*(delx(i,1)*delz/(deltay/2));
        E(i,1) = -(D(i,1)+F(i,1)+C(i,1)+G(i,1)+((V(i,1)*phi(i,1)*ct)/delt)+(wi));
    else
        E(i,1) = -(D(i,1)+F(i,1)+C(i,1)+G(i,1)+((V(i,1)*phi(i,1)*ct)/delt));
    end
end
%Assining values to the A matrix
for i = 1:imax*jmax
    A(i,i) = E(i,1);
    if i < imax*jmax
        A(i,i+1) = F(i,1);
    end
    if i>1
        A(i,i-1) = D(i,1);
    end
    if i < ((imax*jmax)-imax+1)
        A(i,i+imax) = G(i,1);
    end
    if i > imax
        A(i,i-imax) = C(i,1);
    end
end
%Initial pressure for the matrix
for i = 1:imax*jmax
    p(i,1) = p_i;
end
%Solving for pressure and flow rate
for j = 1:t_step
   for i = 1:imax*jmax
    if i == 1
        R(i,j) = -(wi*p_well+(((V(i,1)*phi(i,1)*ct)/delt)*p(i,j)));
    else
        R(i,j) = -(((V(i,1)*phi(i,1)*ct)/delt)*p(i,j));
    end
   end
   p(:,j+1) = A\R(:,j);
   q_calc(j,1) = -4*(wi*(p(1,j+1)-p_well))/5.615;%flow from all adjacent matrix into fracture
end
sqrt_t = zeros(t_step,1);
t = zeros(t_step,1);
p4t = zeros(t_step,1);
for i = 1:t_step
    sqrt_t(i,1) = sqrt(delt*24*i);
    t(i,1) = delt*i;
    p4t(i,1) = t(i,1)^0.25;
end
q_cumm = zeros(t_step,1);
q_cumm(1,1) = -q_calc(1,1);
for i = 2:t_step
    q_cumm(i,1) = -(-q_cumm(i-1,1)+q_calc(i,1));
end
del_pq = zeros(t_step,1);
for i = 1:t_step
    del_pq(i,1) = (p(1,1)-p_well)/(-q_calc(i,1));
end
y = -q_calc;
figure(1)
plot (t,y);
xlabel('Time (days)');
ylabel('Flow rate (bbl/day)');
title({'Flow rate from well using logarithmic griding';' '});

















