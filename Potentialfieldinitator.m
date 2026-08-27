function [DesiredTorque3, DesiredTorque4, Kua, Kla, Kw, ualua, ualla, ualw, calua, calla, calw, BBBB, ThetauapdN1, ThetauapdN2 ,ThetalapdN, ThetawpdN] = Potentialfieldinitator(trajectory, trajectory2, model, sigma,DesiredTorque,DesiredTorque2)   
    
    for jj = 1:length(trajectory(1,:))
        M = getM(model, trajectory(:,jj));
        C = getC(model, trajectory(:,jj));
        G = getG(model, trajectory(:,jj));
        DesiredTorque = [DesiredTorque inverseDynamics(trajectory(:,jj), M, C, G)];
        DesiredTorque2 = [DesiredTorque2 inverseDynamics2(trajectory2(:,jj), M, C, G)];
    end
    DesiredTorque3 = DesiredTorque;
    DesiredTorque4 = DesiredTorque2;
    %qtotal = [1 0 0 0];
    mpm=length(trajectory(1,:));
    Kua = zeros (3,3*mpm);
    Kla = zeros (2,2*mpm);
    Kw = zeros (2,2*mpm);
    Torqueua = DesiredTorque2(1:3,:);
    Torqueuat = zeros (1,mpm); 
    Torqueuan1 = zeros (1,mpm); 
    Torqueuan2 = zeros (1,mpm); 
    Torquela = DesiredTorque2(4:5,:);
    Torquelat = zeros (1,mpm); 
    Torquelan = zeros (1,mpm); 
    Torquew = DesiredTorque2(6:7,:);
    Torquewt = zeros (1,mpm); 
    Torquewn = zeros (1,mpm); 
    Thetauap = trajectory (1:3,:);
    Thetalap = trajectory (4:5,:);
    Thetawp = trajectory (6:7,:);
    ThetauapdT=zeros(3,mpm-1);
    ThetalapdT=zeros(2,mpm-1);
    ThetawpdT=zeros(2,mpm-1);
    for im=1:mpm-1
        ThetauapdT(:,im)= Thetauap(:,im+1)-Thetauap(:,im);
        ThetalapdT(:,im)= Thetalap(:,im+1)-Thetalap(:,im);
        ThetawpdT(:,im)= Thetawp(:,im+1)-Thetawp(:,im);
    end
    Thetauapdot = trajectory (8:10,:);
    %Thetauapdotdot = trajectory (15:17,:);
    %ThetauapdT = diff(Thetauap,1,2);
    ThetauapdN1= ThetauapdT;
    ThetauapdN2= ThetauapdT;
    
    Thetalapdot = trajectory (11:12,:);
    %Thetalapdotdot = trajectory (18:19,:);
    %ThetalapdT = diff(Thetalap,1,2);   
    ThetalapdN= ThetalapdT;   
    nmnm = length(ThetauapdT(1,:))-1;
    ctd = ThetalapdN(1,:);
    ThetalapdN(1,:)= ThetalapdN(2,:);
    ThetalapdN(2,:)= -ctd;

    
    %ThetawpdT = diff(Thetawp,1,2);
    Thetawpdot = trajectory (13:14,:);
    %Thetawtdott = zeros (2,mpm);
    ThetawpdN= ThetawpdT;    
    ctd = ThetawpdN(1,:);
    ThetawpdN(1,:)= ThetawpdN(2,:);
    ThetawpdN(2,:)= -ctd; 

    for ii = 1:length(ThetauapdT(1,:))
        ThetauapdT (:,ii) = ThetauapdT (:,ii)/norm(ThetauapdT (:,ii));
        %Thetauapdot(:,ii+1) = Thetauapdot(:,ii+1)/norm(Thetauapdot(:,ii+1));
        ThetalapdT (:,ii) = ThetalapdT (:,ii)/norm(ThetalapdT (:,ii));
        ThetawpdT (:,ii) = ThetawpdT (:,ii)/norm(ThetawpdT (:,ii));
        ThetalapdN (:,ii) = ThetalapdN (:,ii)/norm(ThetalapdN (:,ii));
        ThetawpdN (:,ii) = ThetawpdN (:,ii)/norm(ThetawpdN (:,ii));
    end

   
    Tua1 = ThetauapdT (:,1); % Unit tangent
    B = null(Tua1'); 
    ThetauapdN1(:,1) = B(:,1)/norm(B(:,1)); % First normal
    ThetauapdN2(:,1) = B(:,2)/norm(B(:,2)); % Second normal (binormal)
    
    
    

    
    % 1. Pre-allocate for speed

    

    
    % 3. Initialize first frame

    
    % 4. Main loop: Double Reflection (In-place for speed)
    % Each step computes the frame at i+1 from the frame at i
    for imm = 1:nmnm
        % --- Reflection 1: Transport from P(i) to P(i+1) ---
        v1 = Thetauap(:,imm+1) - Thetauap(:,imm);
        c1 = dot(v1,v1);
        
        
        if c1 > 1e-15
            % Reflect current Tangent and Normal to intermediate state
            Ti_mid = ThetauapdT(:,imm) - (2/c1) * (v1 * ThetauapdT(:,imm)') * v1;
            Ni_mid = ThetauapdN1(:,imm) - (2/c1) * (v1 * ThetauapdN1(:,imm)') * v1;
        else
            Ti_mid = ThetauapdT(:,imm);
            Ni_mid = ThetauapdN1(:,imm);
        end
        
        % --- Reflection 2: Align intermediate tangent with T(i+1) ---
        v2 = ThetauapdT(:,imm+1) - Ti_mid;
        c2 = dot(v2,v2);
        
        
        if c2 > 1e-15
            % This reflection aligns the frame to the new tangent exactly
            ThetauapdN1(:,imm+1) = Ni_mid - (2/c2) * (v2 * Ni_mid') * v2;
        else
            ThetauapdN1(:,imm+1) = Ni_mid;
        end
        
        % Final cross product for orthonormality
        ThetauapdN2(:,imm+1) = cross(ThetauapdT(:,imm+1), ThetauapdN1(:,imm+1));
        ThetauapdN2(:,imm+1) = ThetauapdN2(:,imm+1) / norm(ThetauapdN2(:,imm+1)); 
        ThetauapdN1(:,imm+1) = cross(ThetauapdN2(:,imm+1), ThetauapdT(:,imm+1));
        
        % Periodic normalization to prevent numerical drift over long curves
        if mod(imm, 100) == 0
            ThetauapdN1(:,imm+1) = ThetauapdN1(:,imm+1) / norm(ThetauapdN1(:,imm+1));
            ThetauapdN2(:,imm+1) = ThetauapdN2(:,imm+1) / norm(ThetauapdN2(:,imm+1));
        end
    end
    ThetauapdT = [ThetauapdT,ThetauapdT(:,mpm-1)];
    ThetauapdN1 = [ThetauapdN1,ThetauapdN1(:,mpm-1)];
    ThetauapdN2 = [ThetauapdN2,ThetauapdN2(:,mpm-1)];
    ThetalapdT = [ThetalapdT,ThetalapdT(:,mpm-1)];
    ThetalapdN = [ThetalapdN,ThetalapdN(:,mpm-1)];
    ThetawpdT = [ThetawpdT,ThetawpdT(:,mpm-1)];
    ThetawpdN = [ThetawpdN,ThetawpdN(:,mpm-1)];
    for jjj = 1:length(trajectory(1,:))
        AAA1 = dot(Torqueua(:,jjj),ThetauapdT(:,jjj));
        if AAA1 > 0 
            Torqueuat (jjj) = AAA1;
        end
        Torqueuan1 (jjj) = dot(Torqueua(:,jjj),ThetauapdN1(:,jjj));
        Torqueuan2 (jjj) = dot(Torqueua(:,jjj),ThetauapdN2(:,jjj));
        AAA2 = dot(Torquela(:,jjj),ThetalapdT(:,jjj));
        if AAA2 > 0 
            Torquelat (jjj) = AAA2;
        end        
        Torquelan (jjj) = dot(Torquela(:,jjj),ThetalapdN(:,jjj));
        AAA3 = dot(Torquew(:,jjj),ThetawpdT(:,jjj));
        if AAA3 > 0 
            Torquewt (jjj) = AAA3;
        end        
        Torquewn (jjj) = dot(Torquew(:,jjj),ThetawpdN(:,jjj));
    end
    
    Torqueuan1mx = max (Torqueuan1);
    Torqueuan1mn = min (Torqueuan1);
    Torqueuan2mx = max (Torqueuan2);
    Torqueuan2mn = min (Torqueuan2);
    Torquelanmx = max (Torquelan);
    Torquelanmn = min (Torquelan);
    Torquewnmx = max (Torquewn);
    Torquewnmn = min (Torquewn); 


% Parameters
    Kuamx1= 2.5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kuamx2= 2.5;
    Kuamn1= 1.5;
    Kuamn2= 1.5;
    Klamx= 2.0;
    Kwmx= 1.5;
    Klamn= 1.0;
    Kwmn= 5; 
    Kt=200;

        
    Kuan1=(Torqueuan1-Torqueuan1mn)/(Torqueuan1mx-Torqueuan1mn)*(Kuamx1-Kuamn1)+Kuamn1;
    Kuan2=(Torqueuan2-Torqueuan2mn)/(Torqueuan2mx-Torqueuan2mn)*(Kuamx2-Kuamn2)+Kuamn2;
    Klan=(Torquelan-Torquelanmn)/(Torquelanmx-Torquelanmn)*(Klamx-Klamn)+Klamn;
    Kwn=(Torquewn-Torquewnmn)/(Torquewnmx-Torquewnmn)*(Kwmx-Kwmn)+Kwmn;
    
    if any(isnan(Kuan1), 'all')
        Kuan1=0.5*(Kuamx1+Kuamn1)*ones(1,mpm);
    end    
    if any(isnan(Kuan2), 'all')
        Kuan2=0.5*(Kuamx2+Kuamn2)*ones(1,mpm);
    end 
    if any(isnan(Klan), 'all')
        Klan=0.5*(Klamx+Klamn)*ones(1,mpm);
    end
    if any(isnan(Kwn), 'all')
        Kwn=0.5*(Kwmx+Kwmn)*ones(1,mpm);
    end    
    
   

    for ikk=1:mpm
        Kua(:,3*ikk-2:3*ikk)= [ThetauapdT(:,ikk),ThetauapdN1(:,ikk),ThetauapdN2(:,ikk)]*[Kt,0,0;0,Kuan1(ikk),0;0,0,Kuan2(ikk)]*[ThetauapdT(:,ikk),ThetauapdN1(:,ikk),ThetauapdN2(:,ikk)]';
        Kla(:,2*ikk-1:2*ikk)= [ThetalapdT(:,ikk),ThetalapdN(:,ikk)]*[Kt,0;0,Klan(ikk)]*[ThetalapdT(:,ikk),ThetalapdN(:,ikk)]';
        Kw(:,2*ikk-1:2*ikk)= [ThetawpdT(:,ikk),ThetawpdN(:,ikk)]*[Kt,0;0,Kwn(ikk)]*[ThetawpdT(:,ikk),ThetawpdN(:,ikk)]';    
    end

    Aopt= zeros (nmnm+2, nmnm+2);
    Aopt (nmnm+2,nmnm+2)= -1;
    for kj=1: nmnm+1
        Aopt (kj, kj) = -1;
        Aopt (kj, kj+1) = 1;
    end
    bopt = zeros(nmnm+2,1);
    mkkk=length(Thetauap(1,:));
    wuap=zeros(mkkk,mkkk);
    wlap=zeros(mkkk,mkkk);
    wwp=zeros(mkkk,mkkk);
    %mk= ceil((nmnm+2)/2);
    %Thrauap = Thetauap (:,mk);
    %Thralap = Thetalap (:,mk);
    %Thrawp = Thetawp (:,mk);
   
    for iii=1:mkkk
        for jjj=1:mkkk
            wuap (jjj,iii)= exp(-(1/(sigma^2))*(Thetauap(:,jjj)-Thetauap(:,iii))'*(Thetauap(:,jjj)-Thetauap(:,iii)));
            wlap (jjj,iii)= exp(-(1/(sigma^2))*(Thetalap(:,jjj)-Thetalap(:,iii))'*(Thetalap(:,jjj)-Thetalap(:,iii)));
            wwp (jjj,iii)= exp(-(1/(sigma^2))*(Thetawp(:,jjj)-Thetawp(:,iii))'*(Thetawp(:,jjj)-Thetawp(:,iii)));
        end
    end
    for jjj=1:mkkk
        wuap(jjj,:) = wuap(jjj,:)/sum(wuap(jjj,:));
        wlap(jjj,:) = wlap(jjj,:)/sum(wlap(jjj,:));
        wwp(jjj,:) = wwp(jjj,:)/sum(wwp(jjj,:));
    end

    wuap2=zeros(3,mkkk);
    wlap2=zeros(2,mkkk);
    wwp2=zeros(2,mkkk);
    for ijj=1:mkkk
        for iji=1:mkkk
            wuap2(:,ijj) = wuap2(:,ijj)+wuap(ijj,iji)*Thetauap(:,iji);
            wlap2(:,ijj) = wlap2(:,ijj)+wlap(ijj,iji)*Thetalap(:,iji);
            wwp2(:,ijj) = wwp2(:,ijj)+wwp(ijj,iji)*Thetawp(:,iji);
        end
    end
    fgp=mkkk;
    duua= zeros (3*fgp,1);
    dula= zeros (2*fgp,1);
    duw= zeros (2*fgp,1); 
    cuua= zeros (3*fgp,fgp);
    cula= zeros (2*fgp,fgp);
    cuw= zeros (2*fgp,fgp);
    for ijj=1:fgp
        for ijk=1:fgp
            cuua(3*ijk-2:3*ijk,ijj)=(wuap(ijk,ijj)/sigma^2)* (Thetauap(:,ijj)-wuap2(:,ijk));
            cula(2*ijk-1:2*ijk,ijj)=(wlap(ijk,ijj)/sigma^2)* (Thetalap(:,ijj)-wlap2(:,ijk));
            cuw(2*ijk-1:2*ijk,ijj)=(wwp(ijk,ijj)/sigma^2)* (Thetawp(:,ijj)-wwp2(:,ijk));
            duua(3*ijk-2:3*ijk,1)=duua(3*ijk-2:3*ijk,1)-((wuap(ijk,ijj)/sigma^2)*(0.5*(Thetauap(:,ijk)-Thetauap(:,ijj))'*Kua(:,3*ijj-2:3*ijj)*(Thetauap(:,ijk)-Thetauap(:,ijj))) *(Thetauap(:,ijj)-wuap2(:,ijk))+wuap(ijk,ijj)*Kua(:,3*ijj-2:3*ijj)*(Thetauap(:,ijk)-Thetauap(:,ijj))+Torqueuat(ijk)*ThetauapdT(:,ijk));
            dula(2*ijk-1:2*ijk,1)=dula(2*ijk-1:2*ijk,1)-((wlap(ijk,ijj)/sigma^2)*(0.5*(Thetalap(:,ijk)-Thetalap(:,ijj))'*Kla(:,2*ijj-1:2*ijj)*(Thetalap(:,ijk)-Thetalap(:,ijj))) *(Thetalap(:,ijj)-wlap2(:,ijk))+wlap(ijk,ijj)*Kla(:,2*ijj-1:2*ijj)*(Thetalap(:,ijk)-Thetalap(:,ijj))+Torquelat(ijk)*ThetalapdT(:,ijk));
            duw(2*ijk-1:2*ijk,1)=duw(2*ijk-1:2*ijk,1)-((wwp(ijk,ijj)/sigma^2)*(0.5*(Thetawp(:,ijk)-Thetawp(:,ijj))'*Kw(:,2*ijj-1:2*ijj)*(Thetawp(:,ijk)-Thetawp(:,ijj))) *(Thetawp(:,ijj)-wwp2(:,ijk))+wwp(ijk,ijj)*Kw(:,2*ijj-1:2*ijj)*(Thetawp(:,ijk)-Thetawp(:,ijj))+Torquewt(ijk)*ThetawpdT(:,ijk));
        end
    end

    Aequa = cuua(3*fgp-2:3*fgp,:);
    Aeqla= cula(2*fgp-1:2*fgp,:);
    Aeqw= cuw(2*fgp-1:2*fgp,:); 
    Bequa = duua (3*fgp-2:3*fgp,:);
    Beqla = dula (2*fgp-1:2*fgp,:);
    Beqw = duw (2*fgp-1:2*fgp,:);

    [ualua, result, status] = gurobi_lsqlin(cuua, duua, Aopt, bopt, Aequa, Bequa, [], []);
    [ualla, result, status] = gurobi_lsqlin(cula, dula, Aopt, bopt, Aeqla, Beqla, [], []);
    [ualw, result, status] = gurobi_lsqlin(cuw, duw, Aopt, bopt, Aeqw, Beqw, [], []);
    %ualua = lsqlin (cuua, duua, Aopt, bopt, Aequa, Bequa);
    %ualla = lsqlin (cula, dula, Aopt, bopt, Aeqla, Beqla);
    %ualw = lsqlin (cuw, duw, Aopt, bopt, Aeqw, Beqw);

    Thetauapdott = zeros (3,fgp);
    Thetauapdott2 = Thetauapdott;
    Thetalapdott = zeros (2,fgp);
    Thetalapdott2 = Thetalapdott;
    Thetawpdott = Thetalapdott;
    Thetawpdott2= Thetawpdott;
    ThetauapdottN1 = Thetauapdott;
    ThetauapdottN2 = ThetauapdottN1;
    for ijk=1:fgp
        Thetauapdott(:,ijk)= 0.9975*Thetauapdot*wuap(ijk,:)';
        Thetauapdott2(:,ijk)= Thetauapdott(:,ijk)/norm( Thetauapdott(:,ijk));
        Thetalapdott(:,ijk)= 0.9975*Thetalapdot*wlap(ijk,:)';
        Thetalapdott2(:,ijk)= Thetalapdott(:,ijk)/norm(Thetalapdott(:,ijk));
        Thetawpdott(:,ijk)= 0.9975*Thetawpdot*wwp(ijk,:)';
        Thetawpdott2(:,ijk)= Thetawpdott(:,ijk)/norm(Thetawpdott(:,ijk));
    end
    ThetalapdottN= Thetalapdott2;
    ctd = ThetalapdottN(1,:);
    ThetalapdottN(1,:)= ThetalapdottN(2,:);
    ThetalapdottN(2,:)= -ctd;
    
    ThetawpdottN= Thetawpdott2;
    ctd = ThetawpdottN(1,:);
    ThetawpdottN(1,:)= ThetawpdottN(2,:);
    ThetawpdottN(2,:)= -ctd;
    
    Tua1 = Thetauapdott(:,1); 
    BBBB = null(Tua1'); 
    ThetauapdottN1(:,1) = BBBB(:,1)/norm(BBBB(:,1)); % First normal
    ThetauapdottN2(:,1) = BBBB(:,2)/norm(BBBB(:,2)); % Second normal (binormal)
    
    for ipp=2:fgp
        [ThetauapdottN1(:,ipp),ThetauapdottN2(:,ipp)]= nextBishopFrame(Thetauap(:,ipp-1),Thetauapdott2(:,ipp-1), Thetauapdott2(:,ipp),Thetauap(:,ipp), ThetauapdottN1(:,ipp-1));
    end
    
    for ijj=1:fgp
        cuvua(3*ijj-2:3*ijj,1)= -(Thetauapdott2(:,ijj)*Thetauapdott2(:,ijj)')*Thetauapdot(:,ijj)+Thetauapdott(:,ijj);
        cuvua(3*ijj-2:3*ijj,2)=-(ThetauapdottN1(:,ijj)*ThetauapdottN1(:,ijj)')*Thetauapdot(:,ijj);
        cuvua(3*ijj-2:3*ijj,3)=-(ThetauapdottN2(:,ijj)*ThetauapdottN2(:,ijj)')*Thetauapdot(:,ijj);
        cuvla(2*ijj-1:2*ijj,1)= -(Thetalapdott2(:,ijj)*Thetalapdott2(:,ijj)')*Thetalapdot(:,ijj)+Thetalapdott(:,ijj);
        cuvla(2*ijj-1:2*ijj,2)=-(ThetalapdottN(:,ijj)*ThetalapdottN(:,ijj)')*Thetalapdot(:,ijj);
        cuvw(2*ijj-1:2*ijj,1)= -(Thetawpdott2(:,ijj)*Thetawpdott2(:,ijj)')*Thetawpdot(:,ijj)+Thetawpdott(:,ijj);
        cuvw(2*ijj-1:2*ijj,2)=-(ThetawpdottN(:,ijj)*ThetawpdottN(:,ijj)')*Thetawpdot(:,ijj);
        
        duvua(3*ijj-2:3*ijj,1)= -(Torqueua(:,ijj)+cuua(3*ijj-2:3*ijj,:)*ualua-duua(3*ijj-2:3*ijj,1));
        duvla(2*ijj-1:2*ijj,1)= -(Torquela(:,ijj)+cula(2*ijj-1:2*ijj,:)*ualla-dula(2*ijj-1:2*ijj,1));
        duvw(2*ijj-1:2*ijj,1)= -(Torquew(:,ijj)+cuw(2*ijj-1:2*ijj,:)*ualw-duw(2*ijj-1:2*ijj,1));

    end
    
    bv= [0;0;0]; Av= [-1,0,0;0,-1,0;0,0,-1];
    bv2= [0;0]; Av2= [-1,0;0,-1;];
    
    %calua = lsqlin (cuvua, duvua, Av, bv);
    %calla = lsqlin (cuvla, duvla, Av2, bv2);
    %calw = lsqlin (cuvw, duvw, Av2, bv2);
    Aeqc1= [5 -1 0
            0 1 -1
            0 0 0];
        beqc1= [0;0;0];
        Aeqc2= [5 -1
                0 0];
        beqc2 = [0;0];   
   
    [calua, result, status] = gurobi_lsqlin(cuvua, duvua, Av, bv, Aeqc1, beqc1, [10;10;10], []);
    [calla, result, status] = gurobi_lsqlin(cuvla, duvla, Av2, bv2, Aeqc2, beqc2, [10;10], []);
    [calw, result, status] = gurobi_lsqlin(cuvw, duvw, Av2, bv2, Aeqc2, beqc2, [10;10], []);
   %{ 
    calua(1)=0*300;
    calua(2)=300;
    calua(3)=300;
    calla(1)=0*300;
    calla(2)=300;
    calw(1)=0*150;
    calw(2)=150;
    %}
    disp(calua);
    disp(calla);
    disp(calw);
end