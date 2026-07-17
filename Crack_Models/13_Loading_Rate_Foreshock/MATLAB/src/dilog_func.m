function y = dilog_func(x)

% DILOG  Compute real part of the dilogarithm function Li_2(x)
    y = real(polylog(2, x));
end