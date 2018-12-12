function RE=ESSdist_re(P, div,n, par,tval)
global Rm Rd Pdiv Pdeath re1 re2 A B Altrusm %Igen_speed Dremtim_pdiv Iremtim_pdiv ID_rat
%         val1=par(1);
%         val2=par(2);
        lpar=length(par);
        Int_sp=par(1:lpar/2);
        Dam_sp=par(lpar/2+1:lpar);
        % Anticipated time required for damaged protein to reach Pdiv
        Dgen_speed(div)=mean(Dam_sp);% Mean damage protein speed
        Dremdis_pdiv(div)=Pdeath-P(2);%remaining distance
        Dremtim_pdiv(div)=Dremdis_pdiv(div)/Dgen_speed(div);% Required time
        prevtim=tval;%previous division time (same for I and D)
        Dx(div,n)=div-1+Dremtim_pdiv(div)/prevtim;
        
        % Anticipated time required for intact protein to reach Pdiv
        Igen_speed(div)=mean(Int_sp);% Mean intact protein speed
        Iremdis_pdiv(div)=Pdiv-P(1);%distance from starting point of the current generation
        Iremtim_pdiv(div)=Iremdis_pdiv(div)/Igen_speed(div);% Required time
        if div>1
            Idiv_spdRat(div)=Igen_speed(div)/Igen_speed(div-1);
%         else
%             re2=1;
        end
        ID_rat(div)=Iremtim_pdiv(div)/Dremtim_pdiv(div);
        tim_rat=Iremtim_pdiv(div)/Dremtim_pdiv(div);
        if tim_rat<=100
            re2=1-ID_rat(div);
            B(div)=re2;
%         if (re2<0)
%             re2=0;
%         else
            re2=re2+2*(Altrusm-1/2)*((1-Altrusm)*re2+Altrusm*(1-re2));
%         end
        else
            if re2>0.01
                re2=re2-0.1*re2;
            else
                re2=0;
            end
%             B(div-1)=B(div-1)-0.1;
%             re2=B(div-1);
%             B(div)=B(div-1);
        end
%         re2=((log10(Dx(div,n)-(div-1))/log10(Dx(div,n))));
        if (re2<0)
            re2=0;
        end
        if re2>1
            re2=1;
        end

        re1=re2;
        A(n,div)=re2;
RE=re2;