function DIVS=ESSdivs_re(DIV_P)
global npas Pdiv
X(:,1)=DIV_P(1:2);
for k=1:npas-1
    if X(1,k)<Pdiv
        X(:,k+1)=ESS_desol(X(:,k));
    else
        break;
    end
end
DIVS=X(:,k);
