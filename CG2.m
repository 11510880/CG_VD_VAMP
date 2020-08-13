function x_hat = CG2(z, fxnW, m, iter_numb, init)

if init == 0
    x_hat = zeros(m,1);
    r = z;
else
    x_hat = init;
    r = z - fxnW(x_hat);
end

p = r;

alpha = zeros(1, iter_numb);
beta = zeros(1, iter_numb);

for i = 1:iter_numb
    W_p = fxnW(p);
    
    alpha(i) = (r' * r)/(p' * W_p);
    x_hat = x_hat + alpha(i) * p;
    
    r_prev = r;
    r = r - alpha(i) * W_p;
    
    beta(i) = (r' * r)/(r_prev' * r_prev);
    p = r + beta(i) * p;
    
end