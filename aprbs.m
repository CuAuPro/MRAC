function y=aprbs(N, Tp, max_u)
Band = 1/Tp;
u = idinput(N,'prbs',[0 Band],[-1 1]);
d = diff(u);
idx = find(d) + 1;
idx = [1;idx];
for ii = 1:length(idx) - 1
     amp = randn;
     u(idx(ii):idx(ii+1)-1) = amp*u(idx(ii));
end
y = normalize(u, 'range').*max_u;
end
