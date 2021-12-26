# Photovoltaic
MATLAB Code to calculate loads on photovoltaic arrays
This MATLAB-based code, which I surprisingly could not find a template for on GitHub when I initially uploaded on 12/26/21 was written to optimize the design of photovoltaic arrays.
%% Script to plot safety factors for a ground-mount PV system %%
%% INPUTS %%
%% a - wall angle from horizontal,    degrees, 90
%% q - wall angle from horizontal,    radians, pi/2
%% g - roof angle from horizontal,    degrees, 10
%% G - roof angle from horizontal,    radians, pi
%% h - height of bottom of house, inches,  24 inch
%% l - house EW dimension, inches, 57 feet (684 inches)
%% d - house vert dimension, inches, 8 feet (96 inches)
%% w - house NS dimension, inches, 24 feet (192 inches)
%% D - equal to d
%% L - equal to l
%% W - equal to w
%% A - house wall area, square inches (A = L * D)
%% H - house wall center height (H = h + D / 2)
%% U - wind speed, mph, (0, 10, 180)
%% u - wind speed, m/s  (0, ... ~80)
%% n - number of posts, (32)
%% B - beam matrix c I S s Y
%% c - beam height, inches, 6, 8, 10
%% I - area moment, in^4, 6, 8, 10
%% S - section modulus, in^3
%% s - yield strength, psi
%% M - bending moment, lb-in
%% Y - yield moment, lb-in
%% f - force per unit length, N/m per Fage 1927 converted to lb/in
%% F - total force on house (N,U,a)
%% SF - saftey factor matix (h,n,b,a,U,N)

a = 90;  %% populate angle vector (degrees)
q = deg2rad(a);  %% convert array angle to radians
g = 10;  %% roof angle (degrees)
G = deg2rad(g);  %% roof angle (radians)

d = 8 * 12;  %%  house height, inches
l = 57 * 12; %%  house EW dimension, inches
w = 24 * 12; %%  house NS dimension, inches
D = d/39.06; %%  find structure height in meters
W = w/39.06; %%  find structure NS width in meters

U = linspace(0,180,37); %% populate wind vector
% U = 110; %% populate wind vector, mph
u = U/2.23694; %% wind speed vector in meters per second

f = zeros(length(U),length(a));  %% create force matrix
p = 1.225; %% air density, kg/m^3
b = 1;     %% infinite plate height, m

h = 24; %% initialize front (south) height vector (in)
H = zeros(length(h), length(q)); %% initialize centroid height array

for r = 1:length(h)
    for s = 1:length(q)
        H(r,s)= d * sin(q(s)) / 2 + h(r); %% populate centroid height vector (inches)
    end
end

H

N = [32];  %% eight piers per each of four I-beams
F = zeros(length(N),length(U),length(a)); %% total array drag force, lbs
lift = zeros(length(U),length(a)); %% specific lift force, N/m 
Lift = zeros(length(N),length(U),length(a)); %% total array lift force, lbs
M = zeros(length(h),length(N),length(U),length(a)); %% total array moment in-lb
m = zeros(length(h),length(U),length(a));  %% array moment, in-lb/in

L = zeros(1, length(N));  %% initialize array length vector
for kk = 1:length(N)
    L(kk) = l;    %% find EW length of structure, in    
end

%% calculate drag force and lift force

for i = 1:length(u)
    for j = 1:length(a)
        f(i,j) = p * u(i) * u(i) * b * pi * sin(q(j))/(4 + pi * sin(q(j))); %% N/m, 1-m plate
        lift(i,j) = p * u(i) * u(i) * b * pi * sin(G); %% N/m, 1-m plate
    end
end

f = f * D;      %% drag force per meter of structure, N/m
f = f * 0.2248; %% drag force per meter of structure, lbf/m 
f = f / 39.1;   %% drag force per inch of structure,  lbf/in

lift = lift * D;      %% lift force per meter of structure, N/m
lift = lift * 0.2248; %% lift force per meter of structure, lbf/m
lift = lift / 39.1;   %% lift force per inch of structure, lbf/in

for k = 1:length(L)
    F(k,:,:) = f .* L(k);           %% total drag force on structure, lbs
    Lift(k,:,:) = lift .* L(k);     %% total lift force on structure, lbs
end

F             % horizontal force, lbs
Lift          % vertical force, lbs


%% Calculate all bending moments for h, N, U, a  %%
% for mm = 1:length(h)
%     for nn = 1:length(N)
%         for oo = 1:length(U)
%             for pp = 1:length(a)
%                 M(mm,nn,oo,pp) = F(nn,oo,pp) .* H(mm);    %% total moment on array in-lb
%             end
%         end
%     end
% end
% M;

%% Calculate CSA of piers

NumPosts = 32 ;
PostOuterDiameter = 24 ; %% inches
PostRadius = PostOuterDiameter / 2;
SinglePostArea = 3.14159 * PostRadius * PostRadius  %% square inches
SystemPostArea = NumPosts * SinglePostArea %% square inches

%% Calculate Shear-Failure Safety Factor
PostShearStrength = 40000;  %% psi
PostShearForce = PostShearStrength * SystemPostArea  %% lbf 
ShearSafetyFactor = PostShearForce / F  % shear safety factor
