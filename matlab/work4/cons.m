function [c, ceq] = cons(x)
    c = zeros(40, 1);
    alpha = 0.2;
    beta = 20;
    lambda_t = 2*pi/3;  
    mx = 6;
    for k=1:40  
        lambda_k = x(1+(k-1)*mx);
        e_k = x(5+(k-1)*mx);
        c(k) = alpha*exp(-beta*(lambda_k-lambda_t)^2) - e_k;
    end
    ceq = [];
end
    