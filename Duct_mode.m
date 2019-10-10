clear
clc


%Circular duct with mean flow
c0 = 343;         %speed of sound
f = 500;          %frequency of the fan shaft: frequency of the m-lobe pattern
omega = f*2*pi;   %angular frequency
k = omega/c0;     %total wavenumber


Ma = 0.2;         %Mach number of the mean flow in the duct
Uin = Ma*c0;      %sign depends on the rotation direction
a = 1;             %radius of the duct

%Pre-calculation of the zeros of the derivative of bessel function
% tol = 1e-12;
% M = [0:100];     %Bessel mode number
% N = [1:6];     %radial mode number
% [jprimemk,JofJ] = BessDerivZerosBisect2(M,N,tol);
% save jprimemk

%total mode order in the decomposition
m = 36;   %bessel mode, should be smaller than 50
n = 6;   %radsial mode, should be smaller than 6

load jprimemk

%zero matrix for designated modes
alpha_mn = jprimemk(1:m+1,1:n);   %n starts from 0
miu_mn = alpha_mn/a;   %normalize to the problem scale

kz_mn_plus1 = sqrt(k^2-miu_mn.^2);
kz_mn_minus1 = -sqrt(k^2-miu_mn.^2);
kz_mn_plus = 1/(1-Ma^2)*(-k*Ma + sqrt(k^2 - (1-Ma^2)*miu_mn.^2));   %follow the equation in website
kz_mn_minus = 1/(1-Ma^2)*(-k*Ma - sqrt(k^2 - (1-Ma^2)*miu_mn.^2));  %follow the equation in website

%For decomposition
Nr = 20;     %number of points in radial direction
Nt = 80;     %number of points in azimuthal direction

r = a/Nr:a/Nr:a;
theta = 0:2*pi/Nt:2*pi-2*pi/Nt;
%coefficient matrix
A = zeros(Nr*Nt,n*(m+1)+1);
for ir=1:Nr
    for it=1:Nt
        for jm=1:m+1
            for jn=1:n
                A((ir-1)*Nt+it,(jm-1)*n+jn) = besselj(jm-1,miu_mn(jm,jn)*r(ir))*cos(jm*theta(it));
            end
        end
    end
end
%Add plane wave mode   A00 is added to the last column
for ii=1:Nr*Nt
    A(ii,end) = 1;
end

%solve for A_mn
p_in = zeros(Nr*Nt,1);   %this should be known at the beginning; provided by simulation

%over-determined equations Ax=b
x = A\p_in;

%in x stores A_mn with n=0,m=1,2,3,...  at last ,it's A00
%check! fini!



