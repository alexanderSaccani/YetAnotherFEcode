close all
clc
clearvars

%enter parameters
a = 0.2; %horizontal length of beam
theta = 1e-3; %initial inclination of beam in degs

E = 70e9;%Young modulus
alpha = 23e-6;%23e-6;%linear thermal expansion coeff
r = 1e-3; %radius of circular cross section

%dependent parameters
b = a*tand(theta); %arc height
A = 2*pi*r^2; %cross section area
l = sqrt(a^2+b^2);
k = A*E/l;
k_add = k/1e3;

%internal force
f_int = @(T,u) k/(4*l^2)*(2*u.^3+6*b*u.^2+4*b^2*u) + k_add*u - alpha*k*T/2*(2*u+2*b);

dl = @(u) 0.5*(b^2+a^2)^(-0.5)*(u.^2+2*u*b);
V_el = @(T,u) k*dl(u).^2/2 - alpha*k*l*T*dl(u) +1/2*k_add*u.^2;

%T samples
T = [0,20,50,100];

u = linspace(-4*b,2*b,1e3);
u = linspace(-a/10,a/10,1e3);
%internal force
figure
for ii = 1:length(T)
    
  plot(u,f_int(T(ii),u),'linewidth',2);
  hold on;
  grid on
  
end

plot(u,u*0,'color','k','linestyle','--','linewidth',1.5);
legii = cell(length(T),1);
for ii = 1:length(T)
  legii{ii} = num2str(T(ii));
end
legend(legii)

title('internal force')
  
%potential energy
figure
for ii = 1:length(T)
    
  plot(u,V_el(T(ii),u),'linewidth',2);
  hold on;
  grid on
  
end

%plot(u,u*0,'color','k','linestyle','--','linewidth',1.5);
legii = cell(length(T),1);
for ii = 1:length(T)
  legii{ii} = num2str(T(ii));
end
legend(legii)

title('elastic potential energy')


