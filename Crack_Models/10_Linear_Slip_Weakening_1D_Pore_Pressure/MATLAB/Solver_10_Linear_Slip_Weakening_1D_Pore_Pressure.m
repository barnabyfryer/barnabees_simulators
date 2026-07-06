clc
clear all
close all

%% - About
% Nucleation and propagation on pressurized fault with production phase
% Model used in [Fryer et al., JGR: Solid Earth, 2023]
% doi: 10.1029/2022JB025443

%% - Compute Static Equilibrium
% Implementation following:
% Garagash & Germanovich (2012), Journal of Geophysical Research
% doi:10.1029/2012JB009209

%% - Notation

% Material properties
% mu_star  : mu/(1-nu) for Mode II rupture, or mu for Mode III rupture
% mu        : Shear modulus of the host medium
% nu        : Poisson's ratio
% alpha     : Fluid diffusivity along the fault
% f_p       : Peak friction coefficient
% f_r       : Residual friction coefficient
% w         : Slip-weakening rate
% delta_w   : Weakening slip distance, delta_w = f_p / w
% delta_r   : Critical slip distance, delta_r = (f_p - f_r) / w

% Loading properties
% sigma_n    : Normal stress acting on the fault
% p_0        : Initial pore pressure
% Delta_p_inj  : Increase in pore pressure due to injection
% Delta_p_prod : Decrease in pore pressure due to production
% dt_prod      : Duration of the production phase
% sigma_0    : Initial effective (Terzaghi) stress, sigma_n - p_0
% tau_p      : Initial peak shear strength, f_p * sigma_0
% tau_r      : Residual shear strength, f_r * sigma_0

%% - Normalization

% S : Shear stress normalized by tau_p
% P : Pressure normalized by sigma_0
% D : Slip normalized by delta_w
% X : Distance normalized by a_w = (mu_star / tau_p) * delta_w
% T : Time normalized by a_w^2 / alpha

%% - Numerical Method

% Piecewise Constant Slip Method (Appendix A2) is used to account for the
% influence of residual friction.

%% - Inputs
%Decrease of pore pressure during production
dP_prod = -0.7;
%Duration of production phase
dt_prod = 30;
%Increase of pore pressure due to injection
dP_inj = 0.075;
%Residual friction
f_r = 0.6;
%Background stress
tau_b = 0.975;
%Number of points along each crack, discretization
N = 101;

%% - Computed crack lengths
aa = linspace(0,2,843);
bb = linspace(2,30,843);
%Crack lengths to evaluate
a = [aa bb];

%Remove duplicates
a = unique(sort(a));

%% - Basic calculations
%Critical slip
D_r = 1 - f_r;

%% - Define pressure function and derivative
Pi = @(x) erfc(x);
Pi_prime = @(x) -(2/sqrt(pi))*exp(-x.^2);

%% - Preparation
disp('----------- Preparation -----------')
disp('Preliminary computations for Piecewise Method')

%% - Discretization

%Cell indices
i = (0:N-1)';                   %[N,1]
j = (0:N)';                     %[N+1,1]

%Cell size  
dX = 1/N;                       %[1,1]

%Evaluation points
X_i = i*dX;                     %[N,1]
X_j = j*dX;                     %[N+1,1]


a_over_t = [0, logspace(-3,-1,10), logspace(-1,3,130), logspace(3,4,10)];
%Remove duplicates
a_over_t = unique(sort(a_over_t));

k_diff_j = zeros(length(a_over_t), N+1);
dk_diff_j = zeros(length(a_over_t), N+1);

for m = 1:length(a_over_t)

    for n = 1:N+1

        integrand_k  = @(X) Pi(a_over_t(m)*abs(X)) ./ sqrt(1 - X.^2);
        integrand_dk = @(X) Pi_prime(a_over_t(m)*abs(X)) .* abs(X) ./ sqrt(1 - X.^2);

        a1 = max(0, (j(n)-0.5)*dX);
        b1 = min(1, (j(n)+0.5)*dX);

        k_diff_j(m,n)  = (2/pi) * integral(integrand_k, a1, b1);
        dk_diff_j(m,n) = (2/pi) * integral(integrand_dk, a1, b1);

    end
end

%% - Define interpolations
interp_k_diff_j = @(aot) interp1(a_over_t, k_diff_j, aot, 'pchip', 'extrap').';
interp_dk_diff_j = @(aot) interp1(a_over_t, dk_diff_j, aot, 'pchip', 'extrap').';


%% - Elastic influence matrix
%Create row/column index matrices for assembling the influence matrix
[I,J] = ndgrid(i,j);

K_ij = ...
    -1./(2*pi*dX*((I-J).^2-0.25)) ...
    -1./(2*pi*dX*((I+J).^2-0.25));      %[N,N+1]

K_ij(:,1) = -1./(2*pi*dX*(i.^2-0.25));  %[N,N+1]



%% - Constant pressure influence coefficients

k_const_j = zeros(N+1,1);               %[N+1,1]

integrand = @(X) 1./sqrt(1-X.^2);

for n = 1:N+1

    xmin = max(-1,(j(n)-0.5)*dX);
    xmax = min( 1,(j(n)+0.5)*dX);

    k_const_j(n) = (1/pi)*integral(integrand,xmin,xmax);

    if n>1

        xmin = max(-1,(-j(n)-0.5)*dX);
        xmax = min( 1,(-j(n)+0.5)*dX);

        k_const_j(n) = k_const_j(n) + ...
            (1/pi)*integral(integrand,xmin,xmax);

    end

end


%% - Resolution

disp('----------- Resolution ------------')

% Allocate memory
slip = nan(length(a),N);                  %[Na,N]
time = nan(length(a),1);                  %[Na,1]
slip_center = nan(length(a),1);           %[Na,1]

% Solver options
options = optimoptions('fsolve',...
    'Display','off',...
    'SpecifyObjectiveGradient',true);

% Loop over crack lengths
for n = 1:length(a)

    fprintf('Computing Piecewise Method: %3.0f %%\r',100*n/length(a));

    %% - Initial guess

    if n == 1

        % Initial time guess
        t_guess = 1e-2*a(1);

        success = false;

        while ~success && t_guess <= 1e3

            % Initial guess: zero slip everywhere
            x_guess = [zeros(N,1); t_guess];

            % Solve equilibrium equations
            [x,~,exitflag] = fsolve(@(x) Piecewise_in_t(x,...
                a(n),...
                tau_b,...
                dP_inj,...
                dP_prod,...
                dt_prod,...
                K_ij,...
                k_const_j,...
                interp_k_diff_j,...
                interp_dk_diff_j,...
                X_j,...
                D_r,...
                f_r,...
                Pi,...
                Pi_prime), x_guess, options);

            success = exitflag > 0;

            % Increase time guess if solver failed
            t_guess = 2*t_guess;

        end

        if ~success
            fprintf('\nInitialization failed for a = %.3f\n',a(n));
            continue
        end

    else

        % Use previous successful solution as initial guess

        if isnan(time(n-1))
            n_guess = find(~isnan(time),1,'last');
        else
            n_guess = n-1;
        end

        x_guess = [slip(n_guess,:)'; time(n_guess)];

        [x,~,exitflag] = fsolve(...
            @(x) Piecewise_in_t(...
                x,...
                a(n),...
                tau_b,...
                dP_inj,...
                dP_prod,...
                dt_prod,...
                K_ij,...
                k_const_j,...
                interp_k_diff_j,...
                interp_dk_diff_j,...
                X_j,...
                D_r,...
                f_r,...
                Pi,...
                Pi_prime),...
                x_guess,...
                options);

        success = exitflag > 0;

    end

    %% - Check solution

    if success

        slip_solution = x(1:end-1);

        % Uniform slip sign
        valid_sign = all(sign(slip_solution)==sign(slip_solution(1)));

        % Slip decreases away from crack centre
        valid_growth = all(diff(slip_solution)<=0);

        if valid_sign && valid_growth

            slip(n,:) = slip_solution';

            time(n) = x(end);

            slip_center(n) = slip_solution(1);

        end

    end

end

fprintf('Computing Piecewise Method: SUCCESS\n');


%% - Function Piecewise_in_t

function [F,J] = Piecewise_in_t( ...
    x,...
    a,...
    tau_b,...
    dP_inj,...
    dP_prod,...
    dt_prod,...
    K_ij,...
    k_const_j,...
    interp_k_diff_j,...
    interp_dk_diff_j,...
    X_j,...
    D_r,...
    f_r,...
    Pi,...
    Pi_prime)

%% ------------------------------------------------------------------------
% Decompose unknown vector
%% ------------------------------------------------------------------------

D_i = x(1:end-1);                     %[N,1]
D_j = [D_i;0];                        %[N+1,1]
t   = x(end);                         %[1,1]

%% ------------------------------------------------------------------------
% Check time
%% ------------------------------------------------------------------------

check_t_prod = t > -dt_prod;
check_t_inj  = t > 0;

%% ------------------------------------------------------------------------
% Pressure kernel coefficients
%% ------------------------------------------------------------------------

if check_t_prod

    xi = a/sqrt(t^2 + dt_prod^2);

    k_j = k_const_j ...
        - dP_prod*interp_k_diff_j(xi);

    dk_j_dt = ...
        dP_prod*a*t/(t^2+dt_prod^2)^(3/2) ...
        * interp_dk_diff_j(xi);

else

    k_j = k_const_j;

    dk_j_dt = zeros(size(k_j));

end

if check_t_inj

    xi = a/t;

    k_j = k_j ...
        - (dP_inj-dP_prod)*interp_k_diff_j(xi);

    dk_j_dt = dk_j_dt ...
        + (dP_inj-dP_prod)*a/t^2 ...
        * interp_dk_diff_j(xi);

end

%% ------------------------------------------------------------------------
% Slip weakening filters
%% ------------------------------------------------------------------------

filter_residual  = D_i > D_r;
filter_weakening = ~filter_residual;

%% ------------------------------------------------------------------------
% Overpressure
%% ------------------------------------------------------------------------

if check_t_prod

    xi = a*abs(X_j)/sqrt(t^2+dt_prod^2);

    P_j = dP_prod*Pi(xi);

    dP_j_dt = ...
        -dP_prod*a*t/(t^2+dt_prod^2)^(3/2) ...
        .*abs(X_j).*Pi_prime(xi);

else

    P_j = zeros(size(D_j));

    dP_j_dt = zeros(size(D_j));

end

if check_t_inj

    xi = a*abs(X_j)/t;

    P_j = P_j ...
        + (dP_inj-dP_prod)*Pi(xi);

    dP_j_dt = dP_j_dt ...
        - (dP_inj-dP_prod)*a/t^2 ...
        .*abs(X_j).*Pi_prime(xi);

end

P_i    = P_j(1:end-1);
dP_i_dt = dP_j_dt(1:end-1);

%% ------------------------------------------------------------------------
% Elastic equilibrium
%% ------------------------------------------------------------------------

KD_i = K_ij*D_j;

shear_i = tau_b - KD_i/a;

normal_i = 1 - P_i;

%% ------------------------------------------------------------------------
% Friction law
%% ------------------------------------------------------------------------

friction_i = max(f_r,1-D_i);

%Make smooth version of friction, limited to min. fr and reduced by slip
eps = 1e-6;
H = 0.5*(1 - tanh((D_i - D_r)/eps));

friction_i = f_r + (1 - D_i - f_r).*H;

friction_j = [friction_i;1];

%% ------------------------------------------------------------------------
% Residual vector
%% ------------------------------------------------------------------------

F = zeros(length(D_i)+1,1);

F(1:end-1) = shear_i - friction_i.*normal_i;

F(end) = tau_b - k_j'*friction_j;

%% ------------------------------------------------------------------------
% Jacobian
%% ------------------------------------------------------------------------

N = length(D_i);

J = zeros(N+1,N+1);

J(1:end-1,1:end-1) = ...
    -(1/a)*K_ij(:,1:end-1) ...
    + diag(filter_weakening.*normal_i);

J(end,1:end-1) = ...
    (k_j(1:end-1).*filter_weakening)';

J(1:end-1,end) = ...
    friction_i.*dP_i_dt;

J(end,end) = ...
    -dk_j_dt'*friction_j;

end