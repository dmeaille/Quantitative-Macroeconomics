function [A, dist, Anext, C, V] = household_ss(r, w, tax, V, N, par)

[Anext, C, V] = egm(r, w, tax, V, N, par);
T=zeros(par.M*par.N);       % Transition matrix

% Filling in the transition matrix
for i=1:par.M
    [~,~,A_int] = histcounts(Anext(i,:)',par.Agrid);
    for j=1:par.N
        Aport=(Anext(i,j)-par.Agrid(A_int(j)))/(par.Agrid(A_int(j)+1)-par.Agrid(A_int(j)));
        for k=1:par.M
            T((j-1)*par.M+i,(A_int(j)-1)*par.M+k)=T((j-1)*par.M+i,(A_int(j)-1)*par.M+k)+(1-Aport)*par.P(i,k);
            T((j-1)*par.M+i,(A_int(j))*par.M+k)=T((j-1)*par.M+i,(A_int(j))*par.M+k)+(Aport)*par.P(i,k);
        end
    end
end
% Initial distribution
dist=zeros(par.N*par.M,1);
% Initialize mass at the middle
dist(ceil(par.N*par.M/2))=1;
error=1;
it=1;
while error>1e-8
    dist1=T'*dist;
    error=max(max(abs(dist1-dist)));
    dist=dist1;
    it=it+1;
end

A=dist'*kron(par.Agrid',ones(5,1));

end

