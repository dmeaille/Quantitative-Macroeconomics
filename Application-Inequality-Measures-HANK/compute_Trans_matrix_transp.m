function T = compute_Trans_matrix_transp(anext, par)
% find the bin in agrid in which decisions in anext fall - note
% max(anext)=par.Amax-eps, for a small number eps
T=zeros(par.nb_states*par.nba);

for i=1:par.nb_states
    [~,~,A_int] = histcounts(anext(i,:)',par.agrid);
    for j=1:par.nba
        Aport=(anext(i,j)-par.agrid(A_int(j)))/(par.agrid(A_int(j)+1)-par.agrid(A_int(j)));
        for k=1:par.nb_states
            T((j-1)*par.nb_states+i,(A_int(j)-1)*par.nb_states+k)=T((j-1)*par.nb_states+i,(A_int(j)-1)*par.nb_states+k)+(1-Aport)*par.P(i,k);
            T((j-1)*par.nb_states+i,(A_int(j))*par.nb_states+k)=T((j-1)*par.nb_states+i,(A_int(j))*par.nb_states+k)+(Aport)*par.P(i,k);
        end
    end
end

end