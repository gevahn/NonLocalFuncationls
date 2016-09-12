block = 1000;
A = importdata('194.dat');
xyz = A(1:block,1:3);
tic
R = (A(1:block,5)+A(1:block,6))*(A(1:block,5)+A(1:block,6))';
R_A = A(1:block,5)*A(1:block,5)';
R_B = A(1:block,6)*A(1:block,6)';
W = A(1:block,4)*A(1:block,4)';
D = sqrt(squareform(pdist(xyz)));
D(D==0) = 1;
D1 = 1./D;
D1(D1 == 1) = 0;

kF_A = (3*pi)^(1/3)*R_A.^(1/3);
kF_B = (3*pi)^(1/3)*R_B.^(1/3);

kF_A_sym = (0.5*(A(1:block,5).^5)*ones(1,block)+0.5*ones(block,1)*(A(1:block,5)'.^5)).^(1/5);
kF_B_sym = (0.5*(A(1:block,6).^5)*ones(1,block)+0.5*ones(block,1)*(A(1:block,6)'.^5)).^(1/5);

kFAdis = kF_A_sym.* D;
kFBdis = kF_B_sym.* D;

xcholeA = -9*(sin(kFAdis)./(kFAdis.^3)-cos(kFAdis)./(kFAdis.^2));
xcholeB = -9*(sin(kFBdis)./(kFBdis.^3)-cos(kFBdis)./(kFBdis.^2));

integral = 0.5*(sum(sum(R_A.*W.*D1.*xcholeA))+ sum(sum(R_B.*W.*D1.*xcholeB)))
colmb = sum(sum(R.*W.*D1))

toc
% 0.000161