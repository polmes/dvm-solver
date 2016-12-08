clear variables;
% clc;

airfoil = input('4-digit NACA: ','s');
if (length(airfoil) ~= 4)
    error('Invalid NACA airfoil');
end

% Airfoil variables
f = str2double(airfoil(1)) / 100; % max camber 
p = str2double(airfoil(2)) / 10; % max camber position

alpha = input('Angle of attack (º): ');
if (alpha < -10 || alpha > 10)
    warning('Thin Airfoil Theory may not provide accurate results for the chosen angle of attack');
end

M = input('Number of panels: ');
if (M <= 0)
    error('Invalid number of panels');
elseif (0 < M && M < 10)
    warning('Discrete Vortex Method may not provide accurate results for such a low number of panels');
end

% Number of data points
N = M + 1;

dist = input('(a) Uniform\n(b) Full cosine\n(c) Smart\nDistribution for geometry discretization (choose from the above a, b, c): ','s');
if (dist == 'a')
    x = linspace(0,1,N); % uniform
elseif (dist == 'b')
    k = 1:N;
    x = 1/2 * (1 - cos((k-1) / M * pi)); % full cosine
elseif (dist == 'c')
    disp('Coming real soon')
else
    error('Invalid type of distribution'); % make sure x [0:1]
end

flap = input('Want a flap? [y/n] ','s');
if (flap == 'y')
    flap = true;
    x_h = input('Flap hinge location (% chord): ') / 100;
    eta = input('Flap deflection angle (º): ');
else
    flap = false;
end

% Data points (foreach x)
z = zeros(1,N);
if (p > 0)
    reg1 = (x <= p);
    z(reg1) = f/(p^2) * (2*p*x(reg1) - x(reg1).^2);

    reg2 = (p < x);
    z(reg2) = f/((1-p)^2) * (1 - 2*p + 2*p*x(reg2) - x(reg2).^2);
else
    z(:) = 0; % flat plate
end

% Flap
if (flap)
    z_flap = zeros(1,N);
    [~,x_h] = min(abs(x_h-x));
    x_h = x(x_h);
    
    reg3 = x >= x_h;
    z_flap(reg3) = -tan(deg2rad(eta)) * (x(reg3) - x_h);
    plot(x,z_flap);
    z = z + z_flap;
end

% i and i+1
x1 = x(1:M);
x2 = x(2:N);
z1 = z(1:M);
z2 = z(2:N);

% Deltas
delta_x = x2 - x1;
delta_z = z2 - z1;

% Panel chords (length)
c = sqrt((x2 - x1).^2 + (z2 - z1).^2);

% Tangent vectors
t = zeros(2,M);
t(1,:) = (x2 - x1) ./ c; % x's
t(2,:) = (z2 - z1) ./ c; % z's

% Normal vectors
n = zeros(2,M);
n(1,:) = -t(2,:); % x's
n(2,:) = t(1,:); % z's

% Vortex points (aerodynamic centers)
VP = zeros(2,M);
VP(1,:) =  x1 + (x2 - x1) * 1/4; 
VP(2,:) =  z1 + (z2 - z1) * 1/4; 

% Control points
CP = zeros(2,M);
CP(1,:) =  x1 + (x2 - x1) * 3/4; 
CP(2,:) =  z1 + (z2 - z1) * 3/4; 

% System of equations
RHS = zeros(1,M);
A = zeros(M);
for i = 1:M % CP's
    for j = 1:M % vortices
        % Induced velocity at CP (i) due to vortex (j)
        r_sq = (CP(1,i) - VP(1,j))^2 + (CP(2,i) - VP(2,j))^2;
        u = 1/(2*pi) * (CP(2,i) - VP(2,j)) / r_sq;
        w = -1/(2*pi) * (CP(1,i) - VP(1,j)) / r_sq;
        A(i,j) = [u w] * n(:,i);
    end
    RHS(i) = -(cos(deg2rad(alpha)) * n(1,i) + sin(deg2rad(alpha)) * n(2,i));
end

% Induced velocity at CP due to vortex
% v = zeros(2,M);
% for i = 1:M % CP's
%     r_sq = (CP(1,i) - VP(1,:)).^2 + (CP(2,i) - VP(2,:)).^2;
%     v(1,:) = (CP(2,i) - VP(2,:)) ./ r_sq;
%     v(2,:) = -(CP(1,i) - VP(1,:)) ./ r_sq;
%     v = 1/(2*pi) * v;
% end
% AA = v' * n;
% RHSS = -[cos(deg2rad(alpha)) sin(deg2rad(alpha))] * n;

gamma = linsolve(A,RHS');

% Results
Cl = 2 * sum(gamma);
Cm_LE = -2 * sum(gamma' .* VP(1,:) * cos(deg2rad(alpha)));
Cm_AC = -2 * sum(gamma' .* (VP(1,:) - 1/4) * cos(deg2rad(alpha)));

% Plots
% plot(x, z); % mean camber line discretization
% plot(x1,gamma);
