function [Cl,Cm] = dvm(airfoil,alpha,M,dist,flap,x_h,eta)
    if (nargin == 0)
        clc;
        airfoil = input('4-digit NACA: ','s');
        alpha = input('Angle of attack (º): ');
        M = input('Number of panels: ');
        dist = input('(a) Uniform\n(b) Full cosine\n(c) Experimental\nDistribution for geometry discretization (choose from the above a, b, c): ','s');
        flap = input('Include a flap? (y/n) ','s');
    elseif (nargin < 5)
        error('Not enough input arguments');
    elseif (nargin > 7)
        error('Too many input arguments');
    end
    
    if (length(airfoil) ~= 4)
        error('Invalid NACA airfoil');
    end

    % Airfoil camber
    f = str2double(airfoil(1)) / 100; % max camber
    p = str2double(airfoil(2)) / 10; % max camber position
    
    % Airfoil thickness
    t = str2double(airfoil(3:4)) / 100;
    if (t > 0.15)
        warning('Thin Airfoil Theory may not apply for the chosen airfoil thickness');
    end
    
    % Angle of attack
    if (alpha < -10 || alpha > 10)
        warning('Thin Airfoil Theory may not apply for the chosen angle of attack');
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

    % Geometry discretization distribution
    if (dist == 'a')
        x = linspace(0,1,N); % uniform
    elseif (dist == 'b')
        k = 1:N;
        x = 1/2 * (1 - cos((k-1) / M * pi)); % full cosine
    elseif (dist == 'c')
        if (f > 0)
            % Integral of the second derivative
            dz = linspace(2 * f/p,-2 * f/(1 - p),N);
            x = zeros(1,N);

            reg1_opt = (dz >= 0);
            x(reg1_opt) = p * (1 - p * dz(reg1_opt) / (2*f));

            reg2_opt = (dz < 0);
            x(reg2_opt) = p - (p - 1)^2 / (2*f) * dz(reg2_opt);
        else
            x = linspace(0,1,N); % uniform flat plate
        end
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
    
    % Flap distribution
    if (flap)
        % Using this discretization, the first flap data point might not be aligned with the rest
        reg_flap = (x >= x_h);
        z(reg_flap) = z(reg_flap) - tan(eta) * (x(reg_flap) - x_h);
    end
    
    % % i's and j's
    ij = ndgrid(1:M,1:M); % combvec equivalent
    in = zeros(2,M*M); % indices
    in(1,:) = reshape(ij,[1 M*M]);
    in(2,:) = reshape(ij',[1 M*M]);
    
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
    r_sq = (CP(1,in(2,:)) - VP(1,in(1,:))).^2 + (CP(2,in(2,:)) - VP(2,in(1,:))).^2;
    v(1,:) = 1/(2*pi) * (CP(2,in(2,:)) - VP(2,in(1,:))) ./ r_sq;
    v(2,:) = -1/(2*pi) * (CP(1,in(2,:)) - VP(1,in(1,:))) ./ r_sq;

    % Matrices
    A = reshape(sum(v .* n(:,in(2,:))),[M M])';
    RHS = -[cos(alpha) sin(alpha)] * n;

    % Circulation (solve the system of equations)
    gamma = linsolve(A,RHS');

    % Results
    Cl = 2 * sum(gamma); % lift coefficient
    Cp = 2 * gamma' ./ c; % pressure coeficient
    Cm = -2 * sum(gamma' .* VP(1,:) * cos(alpha)); % pitching moment coefficient about the leading edge

    % Plots
    figure;
    
    subplot(3,1,1);
    plot(x,z);
    xlim([0 1]);
    title('Mean camber line discretization');
    xlabel('$$\frac{x}{c}$$','Interpreter','latex');
    ylabel('$$z$$','Interpreter','latex');
    
    subplot(3,1,2);
    bar(x1,gamma,'histc');
    xlim([0 1]);
    title('Circulation distribution');
    xlabel('$$\frac{x}{c}$$','Interpreter','latex');
    ylabel('$$\Gamma$$','Interpreter','latex');
    
    subplot(3,1,3);
    plot(x1,Cp);
    xlim([0 1]);
    title('Pressure distribution');
    xlabel('$$\frac{x}{c}$$','Interpreter','latex');
    ylabel('$$\Delta C_p$$','Interpreter','latex');
end
