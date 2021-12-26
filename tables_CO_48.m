%% Script to plot safety factors for Solar panel arrays on hillside %%
%% Copyright 2019 Human Powered Future PLLC All Rights Reserved %%
%% INPUTS %%
%% a - angle from horizontal,    degrees, 0, 5, ... 90
%% q - angle from horizontal,    radians, 0, ... pi/2
%% h - height of front of array, inches,  12, 24 ... 120
%% l - panel EW dimension, inches, 40 inch default
%% d - panel NS dimension, inches, 80 inch default
%% D - array NS dimension, inches, 160 inch (2*d) default
%% L - array EW dimension, inches, L = N * l / 2 default
%% N - number of modules
%% A - array area, square inches (A = L * D)
%% H - array center height H = h due to slope
%% U - wind speed, mph, (0, 5, 180)
%% u - wind speed, m/s  (0, ... ~80)
%% n - number of posts, (1, 2, ... n = (2 * L / l)
%% B - beam matrix c I S s Y
%% c - beam size, inches, 6, 8, 10
%% I - area moment, in^4, 6, 8, 10
%% S - section modulus, in^3
%% s - yield strength, psi
%% M - bending moment, lb-in
%% Y - yield moment, lb-in
%% f - force per unit length, N/m per Fage 1927 converted to lb/in
%% F - total force on array (N,U,a)
%% SF - saftey factor matix (h,n,b,a,U,N)
%% SP - span matix (h,n,b,a,U,N)

% a = linspace(0,90,19);  %% populate angle vector
% a = linspace(0,90,7);  %% populate angle vector
a = 40;  %% populate angle vector
q = deg2rad(a);  %% convert array angle to radians
d = 66; %% panel NS dimension, inches
l = 40; %% panel EW dimension, inches
D = 2*d/39.1; %% array height in meters
U = 105; %% populate wind vector
% U = linspace(60,180,7); %% populate wind vector
% U = linspace(0,180,37); %% populate wind vector
u = U/2.23694; %% wind speed vector in meters per second

f = zeros(length(U),length(a));  %% create force matrix
p = 1.225; %% air density, kg/m^3
b = 1;     %% infinite plate height, m

h = 18; %% initialize front (south) height vector (in)
% h = linspace(12,120,10); %% initialize front (south) height vector (in)
H = zeros(length(h), length(q)); %% initialize centroid height array

for r = 1:length(h)
    for s = 1:length(q)
%         H(r,s)= d * sin(q(s)) + h(r); %% populate centroid height vector
        H(r,s)= h(r) + 12 ; %% 12 inch rise for Muldoon
    end
end

H

N = [48];  %% 48 panels 
% N = linspace(4,100,25);  %% initialize number of panels vector
F = zeros(length(N),length(U),length(a)); %% total array force, lbs
M = zeros(length(h),length(N),length(U),length(a)); %% total array moment in-lb
m = zeros(length(h),length(U),length(a));  %% array moment, in-lb/in

L = zeros(1, length(N));  %% initialize array length vector
for kk = 1:length(N)
    L(kk) = N(kk) .* l ./ 2;    %% find EW length of solar arrays, in    
end

n = linspace(4,12,9);   %% initialize number of posts vector
B = ones(3,5);  %% initialize beam array, three cross sections, five params
X = size(B); x=X(1); %% number of beam types
SF = zeros(length(h), length(n), x, length(q), length(U), length(N)); % 3Big, 3Small, X,Y,Z,x,y,z
SP = zeros(length(h), length(n), x, length(q), length(U), length(N)); % 3Big, 3Small, X,Y,Z,x,y,z

%% define B matrix

B = [6, 16.4, 5.6, 57400, 1;8, 30, 8, 62000, 1;10, 54, 11, 58800, 1];
for ii=1:3
    B(ii,5) = B(ii,3)*B(ii,4); %% Populate yeild moment field
end



%% calculate bending force on array

% f = p * U * U * b * pi * sin(q)/(4 + pi * sin(q))  %% debug force eqn
% f = p * U(12) * U(12) * b * pi * sin(q(5))/(4 + pi * sin(q(5))) %% debug N/m, 1-m plate
% f(12,5) = p * U(12) * U(12) * b * pi * sin(q(5))/(4 + pi * sin(q(5))) %% debug N/m, 1-m plate
for i = 1:length(u)
    for j = 1:length(a)
        f(i,j) = p * u(i) * u(i) * b * pi * sin(q(j))/(4 + pi * sin(q(j))); %% N/m, 1-m plate
    end
end

f = f * D;      %% force per meter of actual array, N/m
f = f * 0.2248; %% force per meter of actual array, lbf/m 
f = f / 39.1;   %% force per inch of actual array,  lbf/in

for k = 1:length(L)
    F(k,:,:) = f .* L(k);      %% total force on array, lbs
end
F

%% Calculate all bending moments for h, N, U, a  %%
for mm = 1:length(h)
    for nn = 1:length(N)
        for oo = 1:length(U)
            for pp = 1:length(a)
                M(mm,nn,oo,pp) = F(nn,oo,pp) .* H(mm);    %% total moment on array in-lb
            end
        end
    end
end
M

%% Calculate all safety factors for h, n, b, a, U, N
for mm = 1:length(h)
    for nn = 1:length(n)
        for oo = 1:x
            for pp = 1:length(a)
                for qq = 1:length(U)
                    for rr = 1:length(M)
                        SF(mm,nn,oo,pp,qq,rr) = n(nn)*B(oo,5)/M(rr);
                    end
                end
            end
        end
    end
end

SF

%% Calculate all spans for h, n, b, a, U, N
for mm = 1:length(h)
    for nn = 1:length(n)
        for oo = 1:x
            for pp = 1:length(a)
                for qq = 1:length(U)
                    for rr = 1:length(N)
                        SP(mm,nn,oo,pp,qq,rr) = L(rr) ./ n(nn);
                    end
                end
            end
        end
    end
end

SP