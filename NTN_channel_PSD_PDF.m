%% Marginal PDF Plotting for alpha and beta for various beta_ele values
% This script computes the normalized marginal pdfs of alpha and beta based
% on the joint pdf p_{beta_ele}(alpha,beta) for different beta_ele values (10,20,...,90 deg)
clear; clc; close all;

%%
beta_ele_vec = 0:10:90;
H = 65; % Building height in meters
c_0 = 3e8; % Speed of light in m/s
max_relative_delay = [1.0373, 0.8680, 0.6382, 0.4575, 0.3474, 0.2750, ...
 0.2237, 0.1852, 0.1548, 0.1300] * 1e-6; % in seconds

% Define delay spread values 
delay_spread = [250, 183.7667, 125.1762, 85.4138, 63.7133, 50.0438, ...
                40.9588, 34.9798, 31.5052, 30] * 1e-9; % in seconds

% Minimum c=b, which happens at beta=90
c_min = compute_c_at_beta90(H, delay_spread(end));
fprintf('For H = %.1f m and σ_τ = %.1f ns at β = 90°, c_min = %.2f m\n', H, delay_spread(end)*1e9, c_min);
% Initialize arrays to store results
a_values = zeros(size(beta_ele_vec));
b_values = zeros(size(beta_ele_vec));
c_values = zeros(size(beta_ele_vec));

for i = 1:length(beta_ele_vec)
    beta_ele = beta_ele_vec(i); % elevation angle in degrees
    delta_tau_max = max_relative_delay(i); % maximum delay in seconds
    sigma_tau=delay_spread(i);
    % Compute a and c for current elevation angle
    [a, c] = compute_a_c(delta_tau_max, H, beta_ele, c_min);
    b = compute_b(sigma_tau, a, c, beta_ele);

    % Store results
    a_values(i) = a;
    b_values(i) = b;
    c_values(i) = c;
    
    % Display results for each elevation angle
    fprintf('β = %d°: Computed a = %.2f m, b= %.2f m,, c = %.2f m\n', beta_ele, a, b, c);
end

% Create a plot of a and c vs elevation angle
beta_ele_vec_fine=0:0.1:90;
a_values_fine = interp1(0:10:90, a_values, beta_ele_vec_fine, 'pchip');
b_values_fine = interp1(0:10:90, b_values, beta_ele_vec_fine, 'pchip');
c_values_fine = interp1(0:10:90, c_values, beta_ele_vec_fine, 'pchip');
figure;
plot(beta_ele_vec_fine, a_values_fine, 'LineWidth', 2, 'DisplayName', 'Semi-major axis a');
hold on;
plot(beta_ele_vec_fine, b_values_fine, 'LineWidth', 2, 'DisplayName', 'Semi-minor axis c');
plot(beta_ele_vec_fine, c_values_fine, 'LineWidth', 2, 'DisplayName', 'Semi-minor axis b');
grid on;
xlabel('Elevation Angle β (degrees)');
ylabel('Axis Length (meters)');
title('Semi-Ellipsoid Parameters vs Elevation Angle');
legend('Location', 'best');
saveas(gcf, 'a_b_c_DS.png');



% sqrt(15/16)
% Define beta_ele values in degrees:
beta_ele_deg_vec = [0 30 60 90];
nBetaEle = length(beta_ele_deg_vec);

% Integration limits (in radians)
alpha_min = deg2rad(-180);  alpha_max = deg2rad(180);      % alpha: -180 to 180 deg
beta_min  = 0;    beta_max  = pi/2;      % beta: 0 to 90 deg

% Define evaluation grids for alpha and beta marginals
alpha_vals = linspace(alpha_min, alpha_max, 200);  % For marginal p(alpha)
beta_vals  = linspace(beta_min, beta_max, 200);      % For marginal p(beta)
% Define a grid for alpha and beta
alpha_grid = linspace(alpha_min, alpha_max, 1000);
beta_grid = linspace(beta_min, beta_max, 1000);
[Alpha, Beta] = meshgrid(alpha_grid, beta_grid);
    
% Preallocate matrices for marginal pdfs 
p_alpha_all = zeros(nBetaEle, length(alpha_vals));
p_beta_all  = zeros(nBetaEle, length(beta_vals));

%% Autocorrelation Para
% Define the function G(alpha,beta). For example, unity gain:
G = @(alpha,beta) 1;

% Define tau range (in seconds, for example) for R_E will be computed:
hz = 1000;  % Sampling rate: 1 kHz (1000 Hz)
period = 1 / hz;  % Sampling period: 1 ms (0.001 seconds)
tau_range = 0:period:1-period;  % Time range from 0 to 2 seconds, with steps of 1 ms

% Preallocate R_E(tau)
R_E = zeros(nBetaEle, length(tau_range));
fd = 100;  % Doppler frequency: 100 Hz
omega_d = fd*2*pi;  % Convert Doppler frequency to radians per second (100 rad/s)



%% Loop over beta_ele values
for k = 1:nBetaEle
    % Set beta_ele for current iteration (in radians)
    beta_ele = deg2rad(beta_ele_deg_vec(k) );
    
    
    % %ec
    % %ec
    % a_low = lowa(H, eb, beta_ele); 
    % a_high = higha(H, eb, beta_ele); 
    % if beta_ele_deg_vec(k)==0
    % a(1)=300;
    % else
    % a=(a_low+a_high)/2;
    % end

    % c_mean = compute_value(H, a, beta_ele);
    a=H;
    c = 39;
    b=39;
    % Define the function based on the provided equation
    fprime = @(alpha, beta) ( ...
        ((cos(alpha).*cos(beta).*cos(beta_ele) + sin(beta).*sin(beta_ele)).^2) /a^2 + ...
        (sin(alpha).^2 .* cos(beta).^2) / b^2 + ...
        ((-cos(alpha).*cos(beta).*sin(beta_ele) + sin(beta).*cos(beta_ele)).^2) / c^2 ...
        ).^(-3/2);

    % Define the joint pdf p_joint(alpha,beta)
    p_joint = @(alpha, beta) ( cos(beta) / (2*pi*a*b*c)) .* fprime(alpha, beta);
    
    %% Compute marginal pdf for alpha: p_alpha(alpha) = ∫_{beta=beta_min}^{beta_max} p_joint(alpha,beta) d beta
    p_alpha = zeros(size(alpha_vals));
    parfor i = 1:length(alpha_vals)
        a_val = alpha_vals(i);
        % Integrate over beta for fixed alpha 'a'
        p_alpha(i) = integral(@(beta) p_joint(a_val, beta), beta_min, beta_max, 'RelTol',1e-6);
    end
    
    %% Compute marginal pdf for beta: p_beta(beta) = ∫_{alpha=alpha_min}^{alpha_max} p_joint(alpha,beta) d alpha
    p_beta = zeros(size(beta_vals));
    parfor j = 1:length(beta_vals)
        b = beta_vals(j);
        % Integrate over alpha for fixed beta 'b'
        p_beta(j) = integral(@(alpha) p_joint(alpha, b), alpha_min, alpha_max, 'RelTol',1e-6);
    end
    
    %% Normalize the marginals so that their integrals equal 1
    norm_alpha = trapz(alpha_vals, p_alpha);
    p_alpha = p_alpha / norm_alpha;
    
    norm_beta = trapz(beta_vals, p_beta);
    p_beta = p_beta / norm_beta;
    
    %% Store the computed marginals for plotting later
    p_alpha_all(k,:) = p_alpha;
    p_beta_all(k,:)  = p_beta;
    
    %% Autocorrelation
    % Precompute parts of the integrand that don't depend on tau
    BaseIntegrand = G(Alpha, Beta) .* p_joint(Alpha, Beta);
    
    parfor m = 1:length(tau_range)
        tau = tau_range(m);
        % Compute the tau-dependent part
        Integrand = BaseIntegrand .* exp(1j * omega_d * tau .* cos(Alpha).*cos(Beta));
        % Double integration using trapz
        int_beta = trapz(beta_grid, Integrand, 2);
        R_E(k,m) = trapz(alpha_grid, int_beta);
    end

end

%% Plot the marginal pdf of alpha and beta for each beta_ele value
sz_label = 18; % Font size for labels
sz_legend = 14; % Font size for legend

% Create figure with a slightly larger size
figure;
set(gcf, 'Units', 'inches', 'Position', [0, 0, 10, 8]); % Figure size 10x8 inches

% First Subplot
subplot(2,1,1); hold all
for k = 1:nBetaEle
    plot(alpha_vals*180/pi, p_alpha_all(k,:), 'LineWidth', 4); % Increased line width
end
xlabel('\alpha_l (deg)', 'FontSize', sz_label, 'FontWeight', 'bold');
ylabel('p(\alpha_l)', 'FontSize', sz_label, 'FontWeight', 'bold');
xticks(-180:30:180);
title('Marginal PDF of \alpha_l', 'FontSize', sz_label, 'FontWeight', 'bold');
lgd = legend(arrayfun(@(x) sprintf('\\beta_{ele} = %d^{\\circ}', x), beta_ele_deg_vec, 'UniformOutput', false), ...
    'FontSize', sz_legend, 'FontWeight', 'bold');
lgd.NumColumns = 2;
xlim([-180 180]);
set(gca, 'FontSize', 14, 'FontWeight', 'bold'); % Axis tick labels
grid on;
grid minor;

% Second Subplot
subplot(2,1,2); hold all
for k = 1:nBetaEle
    plot(beta_vals*180/pi, p_beta_all(k,:), 'LineWidth', 4); % Increased line width
end
xlabel('\beta_l (deg)', 'FontSize', sz_label, 'FontWeight', 'bold');
ylabel('p(\beta_l)', 'FontSize', sz_label, 'FontWeight', 'bold');
xticks(0:10:90);
title('Marginal PDF of \beta_l', 'FontSize', sz_label, 'FontWeight', 'bold');
lgd = legend(arrayfun(@(x) sprintf('\\beta_{ele} = %d^{\\circ}', x), beta_ele_deg_vec, 'UniformOutput', false), ...
    'FontSize', sz_legend, 'FontWeight', 'bold', 'Location', 'best');
lgd.NumColumns = 2;
xlim([0 90]);
set(gca, 'FontSize', 14, 'FontWeight', 'bold'); % Axis tick labels
grid on;
grid minor;

% Remove whitespace by adjusting subplot positions
subplot(2,1,1);
set(gca, 'LooseInset', get(gca, 'TightInset') * 0.2);
pos1 = get(gca, 'OuterPosition');
pos1(1) = 0;
pos1(3) = 0.97;
pos1(2) = 0.55;
pos1(4) = 0.45;
set(gca, 'OuterPosition', pos1);

subplot(2,1,2);
set(gca, 'LooseInset', get(gca, 'TightInset') * 0.2);
pos2 = get(gca, 'OuterPosition');
pos2(1) = 0;
pos2(3) = 0.97;
pos2(2) = 0;
pos2(4) = 0.5;
set(gca, 'OuterPosition', pos2);

% Set PaperPosition to match figure size and save
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 10, 8]);
saveas(gcf, 'Marginal_PDFs.png');


color_order = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
color_order = repmat(color_order, 1, ceil(nBetaEle/length(color_order)));

%% PSD
L = 4096;
f = [-L/2:1:L/2-1]./L;
figure; hold all
for k = 1:nBetaEle
    PSD = fftshift(abs(fft(R_E(k,:), L)));
    plot(f*hz/fd, PSD./max(PSD), 'LineWidth', 4, ...
         'Color', color_order(k), ...
         'DisplayName', sprintf('\\beta_{ele} = %d^{\\circ}', beta_ele_deg_vec(k)))
end
% for k = 1:nBetaEle
%     PSD = fftshift(abs(fft(R_E(k,:), L)));
%     plot(f*hz/fd, PSD(end:-1:1)./max(PSD), '--', 'LineWidth', 3, ...
%          'Color', color_order(k), ...
%          'DisplayName', sprintf('\\beta_{ele} = %d^{\\circ}', 180-beta_ele_deg_vec(k)))
% end
xlabel('f_{normalized} = f_l/fd', 'FontSize', 16);
ylabel('normalized S(f_{normalized})', 'FontSize', 16);
xticks(-1:0.2:1)

grid on
xlim([-1 1]);
set(gca, 'FontSize', 16); % Increase axis tick label size


lgd=legend('FontSize', 12,'FontWeight', 'bold');
lgd.NumColumns = 2;
lgd.Position = [0.82, 0.9, 0, 0]; % Adjust these values as needed
set(gcf, 'Position', [1000,100, 800, 500]); % [left, bottom, width, height] in pixels
set(gca, 'FontSize', 18); % Increase axis tick label size
set(gca, 'FontWeight', 'bold'); % Make axis ticks bold
% Save the figure
set(gca, 'LooseInset', max(get(gca,'TightInset'), 0.02)); % Reduce whitespace

% Save the figure
saveas(gcf, 'Power_Spectral_Density.png');



function c = compute_c_at_beta90(H, sigma_tau)
    % compute_c_at_beta90 - Calculate semi-minor axis c at β_ele = 90°
    % 
    % Inputs:
    %   H - Maximum building height in meters
    %   sigma_tau - RMS delay spread in seconds
    %   c0 - Speed of light (default: 3e8 m/s)
    %
    % Output:
    %   c - Value of semi-minor axis at β_ele = 90° in meters
    
    % Set default for speed of light if not provided
    c0 = 3e8; % m/s
    
    
    % At β_ele = 90°, a = H and b = c
    a = H;
    
    % Use optimization to find c that gives the required RMS delay spread
    options = optimset('Display', 'off', 'TolFun', 1e-12);
    
    % Initial guess for c (scale based on delay spread)
    c_init = c0 * sigma_tau * 0.6; 
    
    % Find c value that gives the target RMS delay spread
    c = fminsearch(@(c_val) abs(compute_delay_spread(a, c_val, c_val, 90) - sigma_tau), c_init, options);
end

function sigma_tau = compute_delay_spread(a, b, c, beta_ele_deg)
    c_0=3e8;
    % Convert elevation angle to radians
    beta_ele = deg2rad(beta_ele_deg);
    
    % Set up numerical integration grid
    n_alpha = 100; % Number of points for azimuth integration
    n_beta = 100;  % Number of points for elevation integration
    
    alpha_vals = linspace(0, 2*pi, n_alpha);
    beta_vals = linspace(0, pi/2, n_beta);
    
    d_alpha = 2*pi/n_alpha;
    d_beta = (pi/2)/n_beta;
    
    % Precompute trigonometric values
    cos_beta_ele = cos(beta_ele);
    sin_beta_ele = sin(beta_ele);
    
    % Initialize arrays for numerical integration
    moment1_integrand = zeros(n_alpha, n_beta);
    moment2_integrand = zeros(n_alpha, n_beta);
    
    % Compute integrands at each point in the grid
    for i = 1:n_alpha
        cos_alpha = cos(alpha_vals(i));
        sin_alpha = sin(alpha_vals(i));
        
        for j = 1:n_beta
            cos_beta = cos(beta_vals(j));
            sin_beta = sin(beta_vals(j));
            
            % Calculate r_max for this direction using ellipsoid formula
            denominator = (cos_alpha*cos_beta*cos_beta_ele + sin_beta*sin_beta_ele)^2 / a^2 + ...
                          (sin_alpha*cos_beta)^2 / b^2 + ...
                          (sin_beta*cos_beta_ele - cos_alpha*cos_beta*sin_beta_ele)^2 / c^2;
            
            r_max = 1/sqrt(denominator);
            
            % Calculate term inside expectations
            term = 1 - (cos_alpha*cos_beta*cos_beta_ele + sin_beta*sin_beta_ele);
            
            % First moment integrand
            moment1_integrand(i,j) = term * cos_beta * r_max^4;
            
            % Second moment integrand
            moment2_integrand(i,j) = term^2 * cos_beta * r_max^5;
        end
    end
    
    % Perform numerical integration
    moment1 = (3/(8*pi*a*b*c)) * sum(moment1_integrand(:)) * d_alpha * d_beta;
    moment2 = (3/(10*pi*a*b*c)) * sum(moment2_integrand(:)) * d_alpha * d_beta;
    
    % Calculate RMS delay spread
    sigma_tau = (1/c_0) * sqrt(moment2 - moment1^2);
end

function b = compute_b(sigma_tau, a, c, beta_ele_deg)
% COMPUTE_B Computes semi-axis b using RMS delay spread
%
% Inputs:
%   sigma_tau     - RMS delay spread (s)
%   a             - Semi-axis a (m)
%   c             - Semi-axis c (m)
%   beta_ele      - Elevation angle (radians, 0 <= beta_ele < pi/2)
%   c_0           - Speed of light (m/s)
%
% Output:
%   b             - Semi-axis b (m)
c_0=3e8;
beta_ele=deg2rad(beta_ele_deg);
% Validate inputs
if sigma_tau <= 0
    error('sigma_tau must be positive.');
end
if a <= 0
    error('a must be positive.');
end
if c <= 0
    error('c must be positive.');
end
if beta_ele < 0 || beta_ele >pi/2
    error('beta_ele must be in the range [0, pi/2] radians.');
end

% Compute trigonometric terms
cos_beta = cos(beta_ele);
sin_beta = sin(beta_ele);

% Function to compute r_max(alpha, beta, b)
r_max = @(alpha, beta, b) ( ...
    (cos(alpha) .* cos(beta) * cos_beta + sin(beta) * sin_beta).^2 / a^2 + ...
    (sin(alpha) .* cos(beta)).^2 / b^2 + ...
    (sin(beta) * cos_beta - cos(alpha) .* cos(beta) * sin_beta).^2 / c^2 ...
).^(-0.5);

% Function to compute xi(alpha, beta)
xi = @(alpha, beta) cos(alpha) .* cos(beta) * cos_beta + sin(beta) * sin_beta;

% Function to compute sigma_tau given b
function sigma = compute_sigma(b)
    % Integration limits
    alpha_range = [0, 2*pi];
    beta_range = [0, pi/2];
    
    % Compute I1 = integral of (1 - xi) * cos(beta) * r_max^4
    integrand1 = @(alpha, beta) (1 - xi(alpha, beta)) .* cos(beta) .* r_max(alpha, beta, b).^4;
    I1 = integral2(integrand1, alpha_range(1), alpha_range(2), beta_range(1), beta_range(2));
    
    % Compute I2 = integral of (1 - xi)^2 * cos(beta) * r_max^5
    integrand2 = @(alpha, beta) (1 - xi(alpha, beta)).^2 .* cos(beta) .* r_max(alpha, beta, b).^5;
    I2 = integral2(integrand2, alpha_range(1), alpha_range(2), beta_range(1), beta_range(2));
    
    % Compute E[r_l - x']
    E_rx1 = (3 / (8 * pi * a * b * c)) * I1;
    
    % Compute E[(r_l - x')^2]
    E_rx2 = (3 / (10 * pi * a * b * c)) * I2;
    
    % Compute sigma_tau
    sigma = (1 / c_0) * sqrt(E_rx2 - E_rx1^2);
end

% Objective function to minimize
objective = @(b) abs(compute_sigma(b) - sigma_tau);

% Initial guess for b (use a as a starting point)
b_initial = 0.6*a;

% Optimization options
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

% Bounds for b
lb = 1e-6; % Small positive value to avoid numerical issues
ub = Inf;

% Optimize to find b
[b_opt, fval] = fmincon(objective, b_initial, [], [], [], [], lb, ub, [], options);

% Assign result
b = b_opt;

% Verify computation
sigma_tau_computed = compute_sigma(b);
if abs(sigma_tau_computed - sigma_tau) > 1e-10
    warning('Computed sigma_tau (%.2e s) does not match input (%.2e s).', ...
        sigma_tau_computed, sigma_tau);
end

end


% Define the function for delta_tau_max
function [a,c]=compute_a_c(delta_tau_max,H,beta_ele_deg,c_min)
beta_ele=deg2rad(beta_ele_deg);
c_0 = 3e8; 
if beta_ele_deg ~=0
    fun = @(c) (c * (1 + cos(beta_ele)) / (c_0 * sin(beta_ele) * H)) * ...
        sqrt(H^2 - c^2 * cos(beta_ele)^2) - delta_tau_max;
    
    % Constraint: H^2 - c^2 * cos(beta_ele)^2 >= 0
    c_max = H / abs(cos(beta_ele)); % Maximum possible c, ~65.9851
    c_guess = (H+c_min)/2; % Initial guess, closer to 61.408465 than 24.195104
    
    % Use fzero to find c
    try
        c = fzero(fun, c_guess);
        if c <= c_min
            error('Computed c = %.6f is less than or equal to %d.', c, c_min);
        end
        if c > c_max
            error('Computed c = %.6f exceeds maximum possible value %.6f.', c, c_max);
        end
    catch e
        error('Failed to find c: %s', e.message);
    end
    
    % Step 2: Compute a
    if abs(sin(beta_ele)) < 1e-10
        error('sin(beta_ele) is too close to zero, cannot compute a.');
    end
    
    % a = sqrt((H^2 - c^2 * cos(beta_ele)^2) / sin(beta_ele)^2)
    a = sqrt((H^2 - c^2 * cos(beta_ele)^2) / sin(beta_ele)^2);
else
    % Special case for beta = 0°
    % From equation (22) when beta = 0:
    % delta_tau_max = (1/c_0) * (1 + 1) * (1/a)^(-1/2) = 2*a/c_0
    a = c_0 * delta_tau_max / 2;
    
    % Since at beta = 0, a is the semi-major axis aligned with x-axis
    % and c is the semi-minor axis aligned with z-axis (vertical)
    c = H; % The height constraint directly gives us c
        
end
% Check if a is positive
if a <= 0
    error('Computed a is not positive.');
end

% Verify delta_tau_max
% delta_tau_check = (c * (1 + cos(beta_ele)) / (c_0 * sin(beta_ele) * H)) * ...
%     sqrt(H^2 - c^2 * cos(beta_ele)^2);
% fprintf('Computed delta_tau_max = %.6e s, Expected = %.6e s\n', ...
%     delta_tau_check, delta_tau_max);
end