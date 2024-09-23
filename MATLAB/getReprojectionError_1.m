function energy = getReprojectionError_1(gamma, gamma_bar, F21,x)
gamma_hat     = [x(1,1:2)';1];
gamma_hat_bar = [x(1,3:4)';1];
lambda        = x(1,5);
energy = sum((abs(gamma_hat - gamma)).^2) + ...
         sum((abs(gamma_hat_bar - gamma_bar)).^2) + ...
            lambda*transpose(gamma_hat_bar)*F21*gamma_hat;

end