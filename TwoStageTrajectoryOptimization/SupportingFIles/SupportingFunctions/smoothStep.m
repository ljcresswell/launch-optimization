function s = smoothStep(x, x0, eps)
% smooth transition 0 → 1 around x0
s = 0.5 * (1 + tanh((x - x0) / eps));
end