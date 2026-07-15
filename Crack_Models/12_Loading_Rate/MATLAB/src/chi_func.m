function chi = chi_func(zeta, rs_type, kappa0)

%Check we are in physical domain (i.e., not vr > cs)
if any(~isreal(zeta)) || any(zeta <= 0)
    error('zeta left physical domain');
end

switch rs_type
    case 'slip'
        V = V_func_slip(zeta, kappa0);
    case 'aging'
        V = V_func_aging(zeta, kappa0);
    otherwise
        error('rs_type must be "slip" or "aging".');
end


% χ function of G21's Eq. (2.17)
chi = V .* exp(-zeta);
end


function V = V_func_slip(zeta, kappa0)
V = kappa0 * exp(zeta) ./ my_cube_root(zeta).^2;
end

function V = V_func_aging(zeta, kappa0)
L = dilog_func(1 - exp(zeta));
V = kappa0 * (-exp(zeta) .* my_cube_root(1 - exp(-zeta)) .* my_cube_root(L) ./ my_cube_root(zeta).^4);
end

%Cube root allowing negatives
function root = my_cube_root(x)
root = sign(x).*abs(x).^(1/3);
end

