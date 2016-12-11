function [Cl,Cm_LE,Cm_AC] = dvm(airfoil,alpha,M,dist,flap,x_h,eta)
    % clear variables;
    % clc;
    
    if (nargin == 0)
        airfoil = input('4-digit NACA: ','s');
        alpha = input('Angle of attack (º): ');
        M = input('Number of panels: ');
        dist = input('(a) Uniform\n(b) Full cosine\n(c) Optimal\nDistribution for airfoil geometry discretization (choose from the above a, b, c): ','s');
        flap = input('Include a flap? [y/n] ','s');
    elseif (nargin < 5)
        error('Not enough input arguments');
    elseif (nargin > 7)
        error('Too many input arguments');
    end
        
    if (length(airfoil) ~= 4)
        error('Invalid NACA airfoil');
    end

    % Airfoil variables
    f = str2double(airfoil(1)) / 100; % max camber 
    p = str2double(airfoil(2)) / 10; % max camber position

    % Angle of attack
    if (alpha < -10 || alpha > 10)
        warning('Thin Airfoil Theory may not provide accurate results for the chosen angle of attack');
    end
    alpha = deg2rad(alpha);

    % Number of panels
    M = floor(M);
    if (M <= 0)
        error('Invalid number of panels');    
    elseif (M < 10)
        warning('Discrete Vortex Method may not provide accurate results for such a low number of panels');
    end

    % Number of data points
    N = M + 1;

    % Geometry discretization distributio
    if (dist == 'a')
        x = linspace(0,1,N); % uniform
    elseif (dist == 'b')
        k = 1:N;
        x = 1/2 * (1 - cos((k-1) / M * pi)); % full cosine
    elseif (dist == 'c')
        % Integral of the second derivative
        dz = linspace(2 * f/p,-2 * f/(1 - p),N);
        x = zeros(1,N);

        reg1_opt = (dz >= 0);
        x(reg1_opt) = p * (1 - p * dz(reg1_opt) / (2*f));

        reg2_opt = (dz < 0);
        x(reg2_opt) = p - (p - 1)^2 / (2*f) * dz(reg2_opt);
    else
        error('Invalid type of distribution');
    end
    
    % Flap parameters
    if (flap == 'y' || flap == 'Y')
        flap = true;

        if (M < 2)
            error('Flap analysys requires at least two panels');
        end
        
        if (nargin == 0)
            x_h = input('Flap hinge location (in tenths of chord): ');
            eta = input('Flap deflection angle (º): ');
        elseif (nargin < 7)
            error('Not enough input arguments for flap analysis');
        end
        
        % Hinge location
        x_h = x_h / 10;
        if (x_h > 10)
            error('Invalid hinge position');
        end
        
        % Deflection angle
        eta = deg2rad(eta);
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
    
    % Assemble system of equations
    function [v,n,VP] = assemble(x,z)
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

        % Airfoil with angle of attack
        r_sq = (CP(1,in(2,:)) - VP(1,in(1,:))).^2 + (CP(2,in(2,:)) - VP(2,in(1,:))).^2;
        v(1,:) = 1/(2*pi) * (CP(2,in(2,:)) - VP(2,in(1,:))) ./ r_sq;
        v(2,:) = -1/(2*pi) * (CP(1,in(2,:)) - VP(1,in(1,:))) ./ r_sq;
    end

    % Matrices
    [v,n,VP] = assemble(x,z);
    A = reshape(sum(v .* n(:,in(2,:))),[M M])';
    RHS = -[cos(alpha) sin(alpha)] * n;

    % Circulation (solve the system of equations)
    gamma = linsolve(A,RHS');

    % Coefficients
    Cl = 2 * sum(gamma);
    Cm_LE = -2 * sum(gamma' .* VP(1,:) * cos(alpha));
    Cm_AC = -2 * sum(gamma' .* (VP(1,:) - 1/4) * cos(alpha));

    % Flap analysis
    if (flap)
        % Instead of using airfoil discretization, we will take two uniform distributions
        N_flap = floor(x_h / (x_h + (1 - x_h) / cos(eta)) * N); % whole plates %
        x_flap = union(linspace(0,x_h,N_flap),linspace(x_h,1,N-N_flap+1));

        % Flat plates
        z_flap = zeros(1,N);
        reg_flap = (x_flap >= x_h);
        z_flap(reg_flap) = -tan(eta) * (x_flap(reg_flap) - x_h);

        % Matrices
        [v_flap,n_flap,VP_flap] = assemble(x_flap,z_flap);
        A_flap = reshape(sum(v_flap .* n_flap(:,in(2,:))),[M M])';
        RHS_flap = -n_flap(1,:);

        % Circulation (solve the system of equations)
        gamma_flap = linsolve(A_flap,RHS_flap');

        % Coefficients by superposition
        Cl = Cl + 2 * sum(gamma_flap);
        Cm_LE = Cm_LE - 2 * sum(gamma_flap' .* VP_flap(1,:));
        Cm_AC = Cm_AC - 2 * sum(gamma_flap' .* (VP_flap(1,:) - 1/4));
    end

    % Plots
    % plot(x, z); % mean camber line discretization
    % plot(x1,gamma);
end