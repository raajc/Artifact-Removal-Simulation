function [Phi, spatial] = DMD( X1 , X2 , n)
% A = X2/X1;
% [w, d] = eig(A);
% X2 = A * X1
[U,S,V] = svd( X1, "econ" );
UU = U(:,1:n);
VV = V(:,1:n);
SS = S(1:n,1:n);
iSS = inv(SS);
AA = UU' * X2 * VV * iSS;
[W,D] = eig( AA );
Phi = conj(X2 * VV * iSS * W)';
spatial = inv(conj(W')*iSS'*VV');
% lamdba = diag(D);       % eigen value
% omega = log(lamdba)/(t(2)-t(1)); % log of eigen value
% x1 = X1(:, 5);
% b = pinv(Phi)*x1;
% t_dyn = zeros(n, length(t));
% 
% for ii = 1:size(X1,2)
%    t_dyn(:, ii) = (b.*exp(omega*t(ii))); 
% end
% 
% f_dmd = Phi*t_dyn;
