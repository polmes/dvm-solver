clear variables;
% clc;

airfoil = input('4-digit NACA: ','s');
if (length(airfoil) ~= 4)
    error('Invalid NACA airfoil');
end

% Airfoil variables
f = str2double(airfoil(1)) / 100; % max camber 
p = str2double(airfoil(2)) / 10; % max camber position

alpha = input('Angle of attack (�): ');
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
if (flap == 'y' || flap == 'Y')
    flap = true;
    
    x_h = input('Flap hinge location (in tenths of chord): ') / 10;
    if (x_h > 10)
        error('Invalid hinge position');
    end
    
    eta = input('Flap deflection angle (�): ');
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
    z_flap(reg3) = -tan(deg2rad(eta)) * (x(reg3) - x_h);
    
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
tic
RHS = zeros(1,M);
RHS_alpha = zeros(1,M);
A = zeros(M);
A_alpha = A;
for i = 1:M % CP's
    for j = 1:M % vortices
        % Induced velocity at CP (i) due to vortex (j)
        r_sq = (CP(1,i) - VP(1,j))^2 + (CP(2,i) - VP(2,j))^2;
        u = 1/(2*pi) * (CP(2,i) - VP(2,j)) / r_sq;
        w = -1/(2*pi) * (CP(1,i) - VP(1,j)) / r_sq;
        A(i,j) = [u w] * n(:,i);
        
        r_sq_alpha = (CP(1,i) - VP(1,j))^2;
        A_alpha(i,j) = -1/(2*pi) * (CP(1,i) - VP(1,j)) / r_sq_alpha; % * 1
    end
    RHS(i) = -n(1,i);
    RHS_alpha(i) = -sin(deg2rad(alpha));
end
toc

% Induced velocity at control point (i) due to vortex (j)
tic
in = combvec(1:M,1:M); % i's and j's

% Airfoil (and flap) without angle of attack
r_sq = (CP(1,in(2,:)) - VP(1,in(1,:))).^2 + (CP(2,in(2,:)) - VP(2,in(1,:))).^2;
v(1,:) = 1/(2*pi) * (CP(2,in(2,:)) - VP(2,in(1,:))) ./ r_sq;
v(2,:) = -1/(2*pi) * (CP(1,in(2,:)) - VP(1,in(1,:))) ./ r_sq;

% Flat plate with angle of attack
r_sq_alphae = (CP(1,in(2,:)) - VP(1,in(1,:))).^2;
v_alpha = -1/(2*pi) * (CP(1,in(2,:)) - VP(1,in(1,:))) ./ r_sq_alphae;

% Systems of equations
Ae = reshape(sum(v .* n(:,in(2,:))),[M M])';
RHSe = -n(1,:);
A_alphae = reshape(v_alpha,[M M])';
RHS_alphae(1:M) = -sin(deg2rad(alpha));
toc

% Circulation
gamma = linsolve(A,RHS');
gammae = linsolve(Ae,RHSe');
gamma_alpha = linsolve(A_alpha,RHS_alpha');
gamma_tot = gamma + gamma_alpha; % sum linsolves

% Results
Cl = 2 * sum(gamma_tot);
Cm_LE = -2 * sum(gamma_tot' .* VP(1,:) * cos(deg2rad(alpha)));
Cm_AC = -2 * sum(gamma_tot' .* (VP(1,:) - 1/4) * cos(deg2rad(alpha)));

% Plots
% plot(x, z); % mean camber line discretization
% plot(x1,gamma);
