clear variables;
clc;

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
alpha = deg2rad(alpha);

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
if (flap == 'y' || flap == 'Y')
    flap = true;
    
    x_h = input('Flap hinge location (in tenths of chord): ') / 10;
    if (x_h > 10)
        error('Invalid hinge position');
    end
    
    eta = deg2rad(input('Flap deflection angle (º): '));
else
    flap = false;
end

% Data points (foreach x)
z = zeros(1,N);
if (p > 0)
    reg1 = (x <= p);
    z(reg1) = f/(p^2) * (2*p*x(reg1) - x(reg1).^2);

    reg2 = (x > p);
    z(reg2) = f/((1-p)^2) * (1 - 2*p + 2*p*x(reg2) - x(reg2).^2);
else
    z(:) = 0; % flat plate
end

% Flap
if (flap)
    [~,x_h] = min(abs(x_h-x));
    x_h = x(x_h);
    
    z_flap = zeros(1,N);
    reg3 = (x >= x_h);
    z_flap(reg3) = -tan(eta) * (x(reg3) - x_h);
    
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
c = sqrt(delta_x.^2 + delta_z.^2);

% Normal vectors
n = zeros(2,M);
n(1,:) = -delta_z ./ c; % x's
n(2,:) = delta_x ./ c; % z's

% Vortex points (aerodynamic centers)
VP = zeros(2,M);
VP(1,:) =  x1 + delta_x * 1/4; 
VP(2,:) =  z1 + delta_z * 1/4; 

% Control points
CP = zeros(2,M);
CP(1,:) =  x1 + delta_x * 3/4; 
CP(2,:) =  z1 + delta_z * 3/4; 

% Induced velocity at control point (i) due to vortex (j)
in = combvec(1:M,1:M); % i's and j's

% Airfoil (and flap) without angle of attack
r_sq = (CP(1,in(2,:)) - VP(1,in(1,:))).^2 + (CP(2,in(2,:)) - VP(2,in(1,:))).^2;
v(1,:) = 1/(2*pi) * (CP(2,in(2,:)) - VP(2,in(1,:))) ./ r_sq;
v(2,:) = -1/(2*pi) * (CP(1,in(2,:)) - VP(1,in(1,:))) ./ r_sq;

% Flat plate with angle of attack
r_sq_alphae = (CP(1,in(2,:)) - VP(1,in(1,:))).^2;
v_alpha = -1/(2*pi) * (CP(1,in(2,:)) - VP(1,in(1,:))) ./ r_sq_alphae;

% Systems of equations
A = reshape(sum(v .* n(:,in(2,:))),[M M])';
RHS = -n(1,:);
A_alpha = reshape(v_alpha,[M M])';
RHS_alpha(1:M) = -sin(alpha);

% Circulation
gamma = linsolve(A,RHS') + linsolve(A_alpha,RHS_alpha'); % superposition

% Results
Cl = 2 * sum(gamma);
Cm_LE = -2 * sum(gamma' .* VP(1,:) * cos(alpha));
Cm_AC = -2 * sum(gamma' .* (VP(1,:) - 1/4) * cos(alpha));

% Plots
% plot(x, z); % mean camber line discretization
% plot(x1,gamma);
