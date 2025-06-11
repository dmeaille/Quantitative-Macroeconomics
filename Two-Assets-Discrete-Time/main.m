clear all

par = parameters();
par.r_a = 0.02;
par.r_b = 0.0053;


%% Two assets VFI (it takes a long time to run)

tic
[A_policy, B_policy, ~] = vfi_TA(par);
C_policy = (1+par.r_a)*par.AAA + (1+par.r_b) * par.BBB - par.g(par.AAA, A_policy) + par.SSS*par.w - A_policy - B_policy;
toc

Anext = A_policy; Bnext = B_policy; Cnext = C_policy; 

%% Figure to see the adjustment area

figure
surf(par.Agrid, par.Bgrid, squeeze(A_policy(:,:,3)) - par.AAgrid)



%% Two Assets Graves (This one is much faster)


tic
[A_policy, B_policy, V] = vfi_Graves(par);
C_policy = (1+par.r_a)*par.AAA + (1+par.r_b) * par.BBB - par.g(par.AAA, A_policy) + par.SSS*par.w - A_policy - B_policy;
toc

%%
% close all
red = 1; % set to 1 to get the adjustment differences
state = 2;
figure
surf(par.Agrid, par.Bgrid, squeeze(C_policy(:,:,state)))
xlabel('Agrid'); % Label for the x-axis
ylabel('Bgrid'); % Label for the y-axis
zlabel('C_policy'); % Label for the z-axis


figure
if red
    surf(par.Agrid, par.Bgrid, squeeze(A_policy(:,:,state))-par.AAgrid)
else
    surf(par.Agrid, par.Bgrid, squeeze(A_policy(:,:,state)))
end
xlabel('Agrid'); % Label for the x-axis
ylabel('Bgrid'); % Label for the y-axis
zlabel('A_policy'); % Label for the z-axis
% zlim([-4, 4])


figure
if red
    surf(par.Agrid, par.Bgrid, squeeze(B_policy(:,:,state))-par.BBgrid)
else
    surf(par.Agrid, par.Bgrid, squeeze(B_policy(:,:,state)))
end
xlabel('Agrid'); % Label for the x-axis
ylabel('Bgrid'); % Label for the y-axis
zlabel('B_policy'); % Label for the z-axis
% zlim([-4, 4])

