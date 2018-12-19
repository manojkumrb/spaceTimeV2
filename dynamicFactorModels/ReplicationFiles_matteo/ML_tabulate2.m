% ML_Tabulate2 - Smartest version of tabulateto be used if variables have same range

% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)

function Table=ML_tabulate2(x)

N=size(x,2);
J=unique(x);
Table=[J zeros(length(J),N)];
for nn=1:N;
    T1=tabulate(x(:,nn));
    jj1=ismember(J,T1(:,1));
    jj2=ismember(T1(:,1),J);
    Table(jj1,1+nn)=T1(jj2,3);
end
