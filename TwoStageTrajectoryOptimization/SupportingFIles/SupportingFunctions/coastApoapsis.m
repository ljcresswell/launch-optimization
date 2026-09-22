function X_apo = coastApoapsis(T_all, X_all)

rdot = X_all(:,3);

first_neg = find(rdot < 0, 1, 'first');
if ~isempty(first_neg)
    back_pos = find(rdot(first_neg:end) > 0, 1, 'first');
    if ~isempty(back_pos)
        cut   = first_neg + back_pos - 2;
        T_all = T_all(1:cut);
        X_all = X_all(1:cut,:);
        rdot  = X_all(:,3);
    end
end

idx = find(rdot(1:end-1) > 0 & rdot(2:end) <= 0, 1, 'first');

if isempty(idx)
    [~, idx] = max(X_all(:,1));
    X_apo = X_all(idx,:);
    return;
end

a1    = rdot(idx);
a2    = rdot(idx+1);
alpha = max(0, min(1, a1/(a1-a2)));
X_apo = X_all(idx,:) + alpha*(X_all(idx+1,:) - X_all(idx,:));

end