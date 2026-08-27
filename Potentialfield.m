function [hist2,UUP2,diff221,diff222,diff223,disua,disla,disw,UUU,Sua2,Sla2,Sw2,prevIdxua2,prevIdxla2,prevIdxw2,wuap2,uuap2,wlap2,ulap2,wwp2,uwp2,Ni_prevua2,Ti_prevua2,Pi_prevua2,thetadrua2,thetadrla2,thetadrw2] = Potentialfield(bbuf, abuf,hist,J_max,UUP,DesiredTorque2,t_now,ThetauapdN1, ThetauapdN2 ,ThetalapdN, ThetawpdN,thetadrua,thetadrla,thetadrw,sim,model,winsize,Pi_prevua,Ti_prevua,Ni_prevua,wuap,uuap,wlap,ulap,wwp,uwp,calua,calla,calw,ualua,ualla,ualw,M,C,G,Sua,Sla,Sw,sigma,mkkk,Thetauap,Thetalap,Thetawp,Thetauapdt,Thetalapdt,Thetawpdt,Thetauapp,Thetalapp,Thetawpp,Thetauappdt,Thetalappdt,Thetawppdt,Kua,Kla,Kw,prevIdxua,prevIdxla,prevIdxw) 
        [diff221,xxxua,disua, prevIdxua2,startIdxua,endIdxua] = getClosestInPathWindow(Thetauapp, Thetauap, prevIdxua, winsize);
        [diff222,xxxla,disla, prevIdxla2,startIdxla,endIdxla] = getClosestInPathWindow(Thetalapp, Thetalap, prevIdxla, winsize);
        [diff223,xxxw,disw, prevIdxw2,startIdxw,endIdxw] = getClosestInPathWindow(Thetawpp, Thetawp, prevIdxw, winsize);  
        %disua=getClosestInPath(Thetauapp, Thetauap);
        %disla=getClosestInPath(Thetalapp, Thetalap);
        %disw= getClosestInPath(Thetawpp, Thetawp);
        
        for iii=startIdxua:endIdxua
            wuap (iii-startIdxua+1)= exp(-(1/(sigma^2))*(Thetauap-Thetauapp(:,iii))'*(Thetauap-Thetauapp(:,iii)));
            uuap(iii-startIdxua+1)= ualua(iii)+0.5*(Thetauap-Thetauapp(:,iii))'*Kua(:,3*iii-2:3*iii)*(Thetauap-Thetauapp(:,iii));            
            %diffth = Thetauap-Thetauapp(:,iii)
        end
        for iii=startIdxla:endIdxla
            wlap (iii-startIdxla+1)= exp(-(1/(sigma^2))*(Thetalap-Thetalapp(:,iii))'*(Thetalap-Thetalapp(:,iii)));
            ulap(iii-startIdxla+1)= ualla(iii)+0.5*(Thetalap-Thetalapp(:,iii))'*Kla(:,2*iii-1:2*iii)*(Thetalap-Thetalapp(:,iii));
            %diffth = Thetalap-Thetalapp(:,iii)
        end
        for iii=startIdxw:endIdxw
            wwp (iii-startIdxw+1)= exp(-(1/(sigma^2))*(Thetawp-Thetawpp(:,iii))'*(Thetawp-Thetawpp(:,iii)));  
            uwp(iii-startIdxw+1)= ualw(iii)+0.5*(Thetawp-Thetawpp(:,iii))'*Kw(:,2*iii-1:2*iii)*(Thetawp-Thetawpp(:,iii));           
        end
        
        %{
        for iii=1:mkkk
            wuap (iii)= exp(-(1/(sigma^2))*(Thetauap-Thetauapp(:,iii))'*(Thetauap-Thetauapp(:,iii)));
            wlap (iii)= exp(-(1/(sigma^2))*(Thetalap-Thetalapp(:,iii))'*(Thetalap-Thetalapp(:,iii)));
            wwp (iii)= exp(-(1/(sigma^2))*(Thetawp-Thetawpp(:,iii))'*(Thetawp-Thetawpp(:,iii)));  
            uuap(iii)= ualua(iii)+0.5*(Thetauap-Thetauapp(:,iii))'*Kua(:,3*iii-2:3*iii)*(Thetauap-Thetauapp(:,iii));
            ulap(iii)= ualla(iii)+0.5*(Thetalap-Thetalapp(:,iii))'*Kla(:,2*iii-1:2*iii)*(Thetalap-Thetalapp(:,iii));
            uwp(iii)= ualw(iii)+0.5*(Thetawp-Thetawpp(:,iii))'*Kw(:,2*iii-1:2*iii)*(Thetawp-Thetawpp(:,iii));
        end 
        %}
        %Thetalapp = Thetalapp;
        wuap= wuap/sum(wuap);
        wlap= wlap/sum(wlap);
        assert(~any(isnan(wlap), 'all'), 'NaN detected in variable wlap!');
        wwp= wwp/sum(wwp);
        %lg= size(Thetalapp(:,startIdxla:endIdxla))
        %lg2= size(wlap)
        %ag= size(Thetauapp(:,startIdxua:endIdxua))
        %ag2= size(wuap)
        wjxpua = Thetauapp(:,startIdxua:endIdxua)*wuap;

        wjxpla = Thetalapp(:,startIdxla:endIdxla)*wlap;
        
        wjxpw = Thetawpp(:,startIdxw:endIdxw)*wwp;
        
        wuap2 = wuap;
        uuap2 = uuap;
        wlap2 = wlap;
        ulap2 = ulap;
        wwp2 = wwp;
        uwp2 = uwp;
        uuak=0;
        ulak=0;
        uwk=0;
%{
        for iii=1:mkkk
            wuapv (iii)= exp(-(1/(sigma^2))*(Thetauap-Thetauapp(:,iii))'*(Thetauap-Thetauapp(:,iii)));
            wlapv (iii)= exp(-(1/(sigma^2))*(Thetalap-Thetalapp(:,iii))'*(Thetalap-Thetalapp(:,iii)));
            wwpv (iii)= exp(-(1/(sigma^2))*(Thetawp-Thetawpp(:,iii))'*(Thetawp-Thetawpp(:,iii)));  
        end  
%}
         assert(~any(isnan(Thetalap), 'all'), 'NaN detected in variable Thetalap!');
         %{
        for iii=1:mkkk
            uuak = uuak-((1/sigma^2)*wuap(iii)*uuap(iii)*(Thetauapp(:,iii)-wjxpua)+wuap(iii)*Kua(:,3*iii-2:3*iii)*(Thetauap-Thetauapp(:,iii)));
            assert(~any(isnan(uuak), 'all'), 'NaN detected in variable uuak!');
            %if any(isnan(uuak(:)))
           %     keyboard; % This forces MATLAB to pause right here
           % end
            ulak = ulak-((1/sigma^2)*wlap(iii)*ulap(iii)*(Thetalapp(:,iii)-wjxpla)+wlap(iii)*Kla(:,2*iii-1:2*iii)*(Thetalap-Thetalapp(:,iii)));
            assert(~any(isnan(ulak), 'all'), 'NaN detected in variable ulak!');
            %if any(isnan(ulak(:)))
           %     keyboard; % This forces MATLAB to pause right here
            %end            
       
            uwk = uwk-((1/sigma^2)*wwp(iii)*uwp(iii)*(Thetawpp(:,iii)-wjxpw)+wwp(iii)*Kw(:,2*iii-1:2*iii)*(Thetawp-Thetawpp(:,iii)));
            assert(~any(isnan(uwk), 'all'), 'NaN detected in variable uwk!');
            %if any(isnan(uwk(:)))
            %    keyboard; % This forces MATLAB to pause right here
            %end              
        end
         %}
        for iii=startIdxua:endIdxua
            uuak = uuak-((1/sigma^2)*wuap(iii-startIdxua+1)*uuap(iii-startIdxua+1)*(Thetauapp(:,iii)-wjxpua)+wuap(iii-startIdxua+1)*Kua(:,3*iii-2:3*iii)*(Thetauap-Thetauapp(:,iii)));
            assert(~any(isnan(uuak), 'all'), 'NaN detected in variable uuak!');
            %if any(isnan(uuak(:)))
           %     keyboard; % This forces MATLAB to pause right here
           % end
        end
        for iii=startIdxla:endIdxla
            ulak = ulak-((1/sigma^2)*wlap(iii-startIdxla+1)*ulap(iii-startIdxla+1)*(Thetalapp(:,iii)-wjxpla)+wlap(iii-startIdxla+1)*Kla(:,2*iii-1:2*iii)*(Thetalap-Thetalapp(:,iii)));
            assert(~any(isnan(ulak), 'all'), 'NaN detected in variable ulak!');
            %if any(isnan(ulak(:)))
           %     keyboard; % This forces MATLAB to pause right here
            %end            
        end
        for iii=startIdxw:endIdxw
            uwk = uwk-((1/sigma^2)*wwp(iii-startIdxw+1)*uwp(iii-startIdxw+1)*(Thetawpp(:,iii)-wjxpw)+wwp(iii-startIdxw+1)*Kw(:,2*iii-1:2*iii)*(Thetawp-Thetawpp(:,iii)));
            assert(~any(isnan(uwk), 'all'), 'NaN detected in variable uwk!');
            %if any(isnan(uwk(:)))
            %    keyboard; % This forces MATLAB to pause right here
            %end 
        end    
                 
        %{
        for iii=1:mkkk
            uuak = uuak+(1/sigma^2)*wuap(iii)*uuap(iii)*(Thetauapp(iii)-wjxpua)+wuap(iii)*Kua(:,3*iii-2:3*iii)*((Thetauap-Thetauapp(:,iii)));
            ulak = ulak+(1/sigma^2)*wlap(iii)*ulap(iii)*(Thetalapp(iii)-wjxpla)+wlap(iii)*Kla(:,2*iii-1:2*iii)*((Thetalap-Thetalapp(:,iii)));
            uwk = uwk+(1/sigma^2)*wwp(iii)*uwp(iii)*(Thetawpp(iii)-wjxpw)+wwp(iii)*Kw(:,2*iii-1:2*iii)*((Thetawp-Thetawpp(:,iii)));
        end 
        %}
        UX=[uuak;ulak;uwk];
        
        %brc=0;
        %{
        for iii=1:mkkk
            wuapv (iii)= exp(-(1/(sigma^2))*(Thetauap-Thetauapp(:,iii))'*(Thetauap-Thetauapp(:,iii)));
            wlapv (iii)= exp(-(1/(sigma^2))*(Thetalap-Thetalapp(:,iii))'*(Thetalap-Thetalapp(:,iii)));
            wwpv (iii)= exp(-(1/(sigma^2))*(Thetawp-Thetawpp(:,iii))'*(Thetawp-Thetawpp(:,iii)));  
        end  
  
        wuapv= wuapv/sum(wuapv);
        wlapv= wlapv/sum(wlapv);
        assert(~any(isnan(wlapv), 'all'), 'NaN detected in variable wlap!');
        wwpv= wwpv/sum(wwpv);
        %}
       %[s1,s11]= size(wuapv)
       %[s2,s22]= size(Thetauappdt)
       
        uacp = Thetauappdt(:,startIdxua:endIdxua)*wuap;
        lacp = Thetalappdt(:,startIdxla:endIdxla)*wlap;
        wcp = Thetawppdt(:,startIdxla:endIdxla)*wwp;

        ua3=0.5*(1-tanh(3*(2*disua/model.lvua-1)));
        la3=0.5*(1-tanh(3*(2*disla/model.lvla-1)));
        w3= 0.5*(1-tanh(3*(2*disw/model.lvw-1)));
        
        tii = round(t_now/(sim.rt*sim.ts));
        if tii==0 
            tii=1;
        end
        
        if tii > length(Thetauapp)  
            tii=length(Thetauapp);
        end
        %distmua= norm(Thetauap-Thetauapp(:,tii));
        %distmla= norm(Thetalap-Thetalapp(:,tii));
        %distmw= norm(Thetawp-Thetawpp(:,tii));

        Sd1= Thetauap-Thetauapp(:,prevIdxua2);
        distmua= norm(Sd1);
        Sd2= ThetauapdN1(:,prevIdxua2);
        distmuan1= dot(Sd1,Sd2);              
        Sd2= ThetauapdN2(:,prevIdxua2);
        distmuan2= dot(Sd1,Sd2);
        Sd1= Thetalap-Thetalapp(:,prevIdxla2);
        Sd2= ThetalapdN(:,prevIdxla2);   
        distmla= norm(Sd1);
        distmlan= dot(Sd1,Sd2);
        Sd1= Thetawp-Thetawpp(:,prevIdxw2);
        Sd2= ThetawpdN(:,prevIdxw2);
        distmw= norm(Sd1);
        distmwn= dot(Sd1,Sd2);  
%{        
        ua4 = 0.5*(1-tanh(3*(2*distmua/(2*model.lvua)-2)));
        la4 = 0.5*(1-tanh(3*(2*distmla/(2*model.lvla)-2)));
        w4= 0.5*(1-tanh(3*(2*distmw/(2*model.lvw)-2)));        
%}        

        

        
        muanmv= Thetauapdt-Thetauappdt(:,tii);
        mlanmv= Thetalapdt-Thetalappdt(:,tii);
        mwnmv= Thetawpdt-Thetawppdt(:,tii);
        
        muan1m1= dot(muanmv,ThetauapdN1(:,tii));
        muan1m2= dot(muanmv,ThetauapdN2(:,tii));
        mlan1m= dot(mlanmv,ThetalapdN(:,tii));
        mwn1m= dot(mwnmv,ThetawpdN(:,tii));
        
        ua51 = 0.5*(1+tanh(3*(2*abs(muan1m1)/(model.lvua21/0.2)-1)));
        ua52 = 0.5*(1+tanh(3*(2*abs(muan1m2)/(model.lvua22/0.2)-1)));
        la5 = 0.5*(1+tanh(3*(2*abs(mlan1m)/(model.lvla2/0.2)-1)));
        w5 = 0.5*(1+tanh(3*(2*abs(mwn1m)/(model.lvw2/0.2)-1)));
        
        muan1= -1*ThetauapdN1(:,tii)*muan1m1;
        muan2= -1*ThetauapdN2(:,tii)*muan1m2;
        mlan= -1*ThetalapdN(:,tii)*mlan1m;
        mwn= -1*ThetawpdN(:,tii)*mwn1m;  
        
        thetadrua3 = ua3*uacp;
        thetadrla3 = la3*lacp;
        thetadrw3 = w3*wcp;
        
        thetadrua4 = thetadrua3/norm(thetadrua3);
        thetadrla4 = thetadrla3/norm(thetadrla3);
        thetadrw4 = thetadrw3/norm(thetadrw3);
 %{       
        dua= Thetauap-Thetauapp(:,prevIdxua2);
        dla= Thetalap-Thetalapp(:,prevIdxla2);
        dw= Thetawp-Thetawpp(:,prevIdxw2);
        
        duan2= -dot(dua, ThetauapdN1(:,prevIdxua2))* ThetauapdN1(:,prevIdxua2);
        duan3= -dot(dua, ThetauapdN2(:,prevIdxua2))* ThetauapdN2(:,prevIdxua2);
        dlan= -dot(dla, ThetalapdN(:,prevIdxla2))* ThetalapdN(:,prevIdxla2);
        dwn= -dot(dw, ThetawpdN(:,prevIdxw2))* ThetawpdN(:,prevIdxw2);
        
        duan2nr= abs(dot(dua, ThetauapdN1(:,prevIdxua2)));
        duan3nr= abs(dot(dua, ThetauapdN2(:,prevIdxua2)));
        dlannr= abs(dot(dla, ThetalapdN(:,prevIdxla2)));
        dwnnr= abs(dot(dw, ThetawpdN(:,prevIdxw2)));
        
        uan2=0.5*(1-tanh(-3*(2*duan2nr/model.lvua21-2)));
        uan3=0.5*(1-tanh(-3*(2*duan3nr/model.lvua22-2)));
        lan=0.5*(1-tanh(-3*(2*dlannr/model.lvla2-2)));
        wn=0.5*(1-tanh(-3*(2*dwnnr/model.lvw2-2)));        
%}

        if ~isnan(thetadrua4)
            thetadrua = thetadrua4;
        end
        if ~isnan(thetadrla4)
            thetadrla = thetadrla4;
        end
        if ~isnan(thetadrw4)
            thetadrw = thetadrw4;
        end        

        %{
        if mlknmn == 1 
            Pi_prevua = Thetauap; %only for the first moment  + every turn
            Ti_prevua = thetadrua/norm(thetadrua);
            [Niua, Biua] = transportFrame(P1, T1, Ni_prevua, Pi_prevua, Ti_prevua);  
            mlknmn2 = 2; 
        end
        %}
        Tla = thetadrla;%/norm(thetadrla);
        Nla=Tla;
        ctd = Nla(1);
        Nla(1)= Nla(2);
        Nla(2)= -ctd; 
        Tw = thetadrw;%/norm(thetadrw);
        Nw=Tw;
        ctd = Nw(1);
        Nw(1)= Nw(2);
        Nw(2)= -ctd; 
    
        
        Tiua = thetadrua;%/norm(thetadrua);
        Piua = Thetauap;
        
        [Niua, Biua] = nextBishopFrame(Pi_prevua, Ti_prevua, Ni_prevua, Piua, Tiua);
        Pi_prevua2=Piua;
        Ni_prevua2= Niua;
        Ti_prevua2 = Tiua;
        Qua=[Tiua,Niua,Biua];
        Dua=Qua*diag([calua(1),calua(2),calua(3)])*Qua';
        Qla=[Tla,Nla];
        Dla=Qla*diag([calla(1),calla(2)])*Qla';
        Qw=[Tw,Nw];
        Dw=Qw*diag([calw(1),calw(2)])*Qw';
    
        zua= Thetauapdt'*thetadrua3;
        zua2= Thetauapdt'*(muan1);
        zua3= Thetauapdt'*(muan2);
        zla= Thetalapdt'*thetadrla3;
        zla2= Thetalapdt'*(mlan);
        zw= Thetawpdt'*thetadrw3;
        zw2= Thetawpdt'*(mwn);

        if Sua<60%Sbua 
            asua=1;
        else
            asua=0;
        end
        
        if Sla<45%Sbla 
            asla=1;
        else
            asla=0;
        end
        
        if Sw<30%Sbw 
            asw=1;
        else
            asw=0;
        end
        
        
        if (Sua<=0 && 0<=zua) || (60<=Sua && zua<=0)
            bsua = 0;
        else
            bsua = 1;
        end
        if (Sua<=0 && 0<=zua) 
            ysua = 0;
        else
            ysua = 1;
        end
        

        
        if (Sla<=0 && 0<=zla) || (45<=Sla && zla<=0)
            bsla = 0;
        else
            bsla = 1;
        end
        if (Sla<=0 && 0<=zla) 
            ysla = 0;
        else
            ysla = 1;
        end
    

        
        
        if (Sw<=0 && 0<=zw) || (30<=Sw && zw<=0)
            bsw = 0;
        else
            bsw = 1;
        end
        if (Sw <=0 && 0<=zw) 
            ysw = 0;
        else
            ysw = 1;
        end
        
 
        
        
       
        if (Sua<=0 && 0<=zua2) || (60<=Sua && zua2<=0)
            bsua2 = 0;
        else
            bsua2 = 1;
        end
        if (Sua<=0 && 0<=zua2) 
            ysua2 = 0;
        else
            ysua2 = 1;
        end       

        if (Sua<=0 && 0<=zua3) || (60<=Sua && zua3<=0)
            bsua3 = 0;
        else
            bsua3 = 1;
        end
        if (Sua<=0 && 0<=zua3) 
            ysua3 = 0;
        else
            ysua3 = 1;
        end  
        
        if (Sla<=0 && 0<=zla2) || (45<=Sla && zla2<=0)
            bsla2 = 0;
        else
            bsla2 = 1;
        end
        if (Sla<=0 && 0<=zla2) 
            ysla2 = 0;
        else
            ysla2 = 1;
        end        
        
        if (Sw<=0 && 0<=zw2) || (30<=Sw && zw2<=0)
            bsw2 = 0;
        else
            bsw2 = 1;
        end
        if (Sw <=0 && 0<=zw2) 
            ysw2 = 0;
        else
            ysw2 = 1;
        end         
        
      
        delta_Sua = asua*Thetauapdt'*Dua*Thetauapdt-bsua*calua(1)*zua;%-bsua2*calua(2)*zua2-bsua3*calua(3)*zua3;      %ua4*
        delta_Sla = asla*Thetalapdt'*Dla*Thetalapdt-bsla*calla(1)*zla;%-bsla2*calla(2)*zla2;                %la4*
        delta_Sw = asw*Thetawpdt'*Dw*Thetawpdt-bsw*calw(1)*zw;%-bsw2*calw(2)*zw2;               %w4*
        Sua2 = Sua + delta_Sua*sim.ts;
        Sla2 = Sla + delta_Sla*sim.ts;
        Sw2 = Sw + delta_Sw*sim.ts;
        %;%

%{
        [AA1,AA2]=size(ysua2*uan2*calua(2)*duan2)
        [AA3,AA4]=size(ysua3*uan3*calua(3)*duan3) 
        [BA1,BA2]=size(ysla2*lan*calla(2)*dlan) 
        [CA1,CA2]=size(ysla2*lan*calla(2)*dlan) 
%}
        
        uvua= -1*Dua*Thetauapdt+ysua*calua(1)*thetadrua3;% +1*ysua2*calua(2)*muan1+1*ysua3*calua(3)*muan2; 
        uvla= -1*Dla*Thetalapdt+ysla*calla(1)*thetadrla3;%+1*ysla2*calla(2)*mlan;
        uvw= -1*Dw*Thetawpdt+ysw*calw(1)*thetadrw3;%+1*ysw2*calw(2)*mwn;
        thetadrua2=thetadrua;
        thetadrla2=thetadrla;
        thetadrw2=thetadrw;
        UV= [uvua;uvla;uvw];
        UU=1*UX+1*UV;
        UU= real_time_jerk_limiter(UU, UUP, sim.ts, J_max);
        UUP2=UU;
        [U_CF, hist2] = butter3_exo_step(UU, hist, bbuf, abuf);
        M = M + 0.01*eye(size(M));
        %UUP2=U_CF;
        UUU= M*(1*UU+1*DesiredTorque2(:,tii))+C+G;

end