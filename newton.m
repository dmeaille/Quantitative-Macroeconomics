function k_new = newton(k,c,Z,alpha,theta,mu,delta,gamma,cash)

h = 1e-4;
f_x = 1;

while abs(f_x) > 0.0001
    
    f_x = Z*((1-alpha)/theta*Z*k^alpha*c^(-gamma))^((1-alpha)/(mu+alpha))*k^alpha+(1-delta)*k-cash;
    fprime_x = (Z*((1-alpha)/theta*Z*(k+h)^alpha*c^(-gamma))^((1-alpha)/(mu+alpha))*(k+h)^alpha+(1-delta)*(k+h)-cash-f_x)/h;
    k_new = k - f_x/fprime_x;
    k = k_new;
end
end