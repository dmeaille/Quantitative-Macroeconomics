function T = compute_Trans_matrix_transp(Anext, par)
% find the bin in Agrid in which decisions in Anext fall - note
% max(Anext)=par.Amax-eps, for a small number eps
T=zeros(par.M*par.N);

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

end