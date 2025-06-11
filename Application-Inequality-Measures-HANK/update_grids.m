
par.Gamma = Gamma; 
par.mu_g = (1+par.Gamma)*par.mu_b;
par.wage=w;
income = [par.mu_g * par.wage * (1-par.tau), par.mu_b * par.wage * (1-par.tau), par.mu_u * par.wage];

par.y = income;
[par.Amesh,par.Ymesh] = meshgrid(par.agrid,par.y);