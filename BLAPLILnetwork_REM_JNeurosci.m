%This is a network model for REM
%11/11/21 11/10/2021 updated IDC to int is 3, updated gGie=0.5, gGii=0.1, IDC=3 to int in the BLA
% BLA A and C type have 4 assemblies
% BLA interneurons have 8 assemblies
% PL/IL pyramidal neurons have 2 assemblies 
% STDP on the interneurons in the BLA from pA/pC && STDP btw IL and int
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% REM condition:
% PFC:DS02PYILS, DS02PYPLS, DS02FSPL, DS02FSIL
% BLA: pyrA/C, intG with gNMDA: *(100-10)/100 from PFC to pA & pC
% gSYN_scl is different; gPBs=2, gBPs=1
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%PFC: PL and IL assigining a number of cells
mcell=2;
numAssems=2;
NPLe=6*mcell;NPLi=2*mcell;totNPLe=NPLe*numAssems;
NILe=6*mcell;NILi=2*mcell;totNILe=NILe*numAssems;
%BLA: assigining a number of A/C type cells and interneurons
numAssemsB=4;
NpA=6*mcell;NpC=3*mcell;Nbi=2*mcell;Ncomp=2;
totNpA=NpA*numAssemsB;totNpC=NpC*numAssemsB;
numAssemsBint=8;
totNbi=Nbi*numAssemsBint;
%--------------------------------------------------------------------------
%PFC network:PL(superficial layer)
%--------------------------------------------------------------------------
KPLee=zeros(totNPLe,totNPLe);
KPLii=ones(NPLi,NPLi)-eye(NPLi,NPLi);%Kii=zeros(Ni,Ni);
KPLei=zeros(totNPLe,NPLi); % E->I
KPLie=zeros(NPLi,totNPLe); % I->E
% connectivity probability
p_ee_PL=1;p_ee_ass_PL=.1;
p_ei_PL=1;p_ie_PL=1;
%--------------------------------------------------------------------------
%                            PFC network: IL(superficial)
%--------------------------------------------------------------------------
KILee=zeros(totNILe,totNILe);
KILii=ones(NILi,NILi)-eye(NILi,NILi);%Kii=zeros(Ni,Ni);
KILei=zeros(totNILe,NILi); % A->I
KILie=zeros(NILi,totNILe); % I->A
% connectivity probability
p_ee_IL=1;p_ee_ass_IL=.1;
p_ei_IL=1;p_ie_IL=1;

%--------------------------------------------------------------------------
%                              BLA network
%--------------------------------------------------------------------------
KAee=zeros(totNpA,totNpA);
KCee=zeros(totNpC,totNpC);
Kbii=zeros(totNbi,totNbi);
KACee=zeros(totNpA,totNpC);
KCAee=zeros(totNpC,totNpA);
KAei=zeros(totNpA,totNbi); % A<->I
KAei_stdp=zeros(totNpA,totNbi); % A<->I
KAie=zeros(totNbi,totNpA); % I->A
KCei=zeros(totNpC,totNbi); % C<->I
KCei_stdp=zeros(totNpC,totNbi); % C<->I
KCie=zeros(totNbi,totNpC); % I->C
% connectivity probability
dumi=1;
p_ee_A=1;p_ee_Aass=0.1;
p_ee_C=1;p_ee_Cass=0.1;
p_ee_AC=1;p_ee_ACass=0.1;
p_ee_CA=1;p_ee_CAass=0.1;
p_ii=1;p_ii_ass=.5;

p_ei_A=1*dumi;p_ei_C=1*dumi;
p_ie_A=1*dumi;p_ie_C=1*dumi;
%--------------------------------------------------------------------------
%BLA-PFC(PL) assemblies network
%--------------------------------------------------------------------------
KPLBAee=zeros(totNPLe,totNpA);KBAPLee=zeros(totNpA,totNPLe);
KPLBCee=zeros(totNPLe,totNpC);KBCPLee=zeros(totNpC,totNPLe);

% Assemblies conectivity btw PFC(PL) and BLA_C(BC)
onoffBPL=1;
p_ee_PLBC=onoffBPL*1;p_ee_PLBCass=onoffBPL*0.1;% assemble PFC -> assemble BLA_C
p_ee_BCPL=onoffBPL*1;p_ee_BCPLass=onoffBPL*0.1;% assemble BLA_C -> assemble PFC
% Assemblies conectivity btw PFC(PL) and BLA_A(BA)
p_ee_PLBA=0;p_ee_PLBAass=0;% assemble PFC -> assemble BLA_A
p_ee_BAPL=0;p_ee_BAPLass=0;% assemble BLA_A -> assemble PFC
%--------------------------------------------------------------------------
%BLA-PFC(IL) assemblies network
%--------------------------------------------------------------------------
KILBAee=zeros(totNILe,totNpA);KBAILee=zeros(totNpA,totNILe);
KILBCee=zeros(totNILe,totNpC);KBCILee=zeros(totNpC,totNILe);
KILBIei=zeros(totNILe,totNbi);
onoffBIL=1;
% Assemblies conectivity btw PFC(IL) and BLA_C(BC)
p_ee_ILBC=1*onoffBIL;p_ee_ILBCass=0.1*onoffBIL;% assemble PFC -> assemble BLA_C
p_ee_BCIL=1*onoffBIL;p_ee_BCILass=0.1*onoffBIL;% assemble BLA_C -> assemble PFC
p_ei_ILBI=1*onoffBIL;
% Assemblies conectivity btw PFC(IL) and BLA_A(BA)
p_ee_ILBA=0;p_ee_ILBAass=0;% assemble PFC -> assemble BLA_A
p_ee_BAIL=0;p_ee_BAILass=0;% assemble BLA_A -> assemble PFC

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%PFC network:PL(superficial layer)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
PYtype='DS02PYPLS';PVtype='DS02FSPL';%DS02PYPLS

% Two-compartment Pyramidal cell model ('Es'=soma, 'Ed'=dendrite)
spec=get_LIM_cell(PYtype,totNPLe);
% One-compartment FS model ('FS' = PV+ interneuron inhibiting PY soma)
spec.populations(end+1)=getfield(get_LIM_cell(PVtype,NPLi),'populations');

for ii=1:numAssems
    for jj=1:numAssems
        if ii==jj
            assemPLE=zeros(NPLe,NPLe);
            nodigNe=(NPLe*NPLe-NPLe);dumA=rand(nodigNe,1)< p_ee_PL;
            me=0;
            for iii=1:NPLe
                for jjj=1:NPLe
                    if jjj ~= iii
                        me=me+1;
                        assemPLE(iii,jjj)=dumA(me);
                    end
                end
            end
            inde=(ii-1)*NPLe+(1:NPLe);
            KPLee(inde,inde)=assemPLE;
        else
            assemPLE=rand(NPLe) < p_ee_ass_PL;
            indei=(ii-1)*NPLe+(1:NPLe);indej=(jj-1)*NPLe+(1:NPLe);
            KPLee(indei,indej)=assemPLE;
        end
    end
end

for ii=1:numAssems
    assemei=rand(NPLe,NPLi) < p_ei_PL;
    assemie=rand(NPLi,NPLe) < p_ie_PL;
    inde=(ii-1)*NPLe+(1:NPLe);
    KPLei(inde,:)=assemei;KPLie(:,inde)=assemie;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                            PFC network: IL(superficial)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
PYtype='DS02PYILS';PVtype='DS02FSIL';%DS02PYPLS
% Two-compartment Pyramidal cell model ('Es'=soma, 'Ed'=dendrite)
spec.populations(end+1:end+Ncomp)=getfield(get_LIM_cell(PYtype,totNILe),'populations');
spec.connections(end+1:end+Ncomp)=getfield(get_LIM_cell(PYtype,totNILe),'connections');
% One-compartment FS model ('FS' = PV+ interneuron inhibiting PY soma)
spec.populations(end+1)=getfield(get_LIM_cell(PVtype,NILi),'populations');

for ii=1:numAssems
    for jj=1:numAssems
        if ii==jj
            assemILE=zeros(NILe,NILe);
            nodigNe=(NILe*NILe-NILe);dumA=rand(nodigNe,1)< p_ee_IL;
            me=0;
            for iii=1:NILe
                for jjj=1:NILe
                    if jjj ~= iii
                        me=me+1;
                        assemILE(iii,jjj)=dumA(me);
                    end
                end
            end
            inde=(ii-1)*NILe+(1:NILe);
            KILee(inde,inde)=assemILE;
        else
            assemILE=rand(NILe) < p_ee_ass_IL;
            indei=(ii-1)*NILe+(1:NILe);indej=(jj-1)*NILe+(1:NILe);
            KILee(indei,indej)=assemILE;
        end
    end
end

for ii=1:numAssems
    assemILei=rand(NILe,NILi) < p_ei_IL;
    assemILie=rand(NILi,NILe) < p_ie_IL;
    inde=(ii-1)*NILe+(1:NILe);
    KILei(inde,:)=assemILei;KILie(:,inde)=assemILie;
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                              BLA network
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
PAtype='pyrA';PCtype='pyrC';Intype='intG';

spec.populations(end+1:end+Ncomp)=getfield(get_BLA_cell(PAtype,totNpA),'populations');
spec.connections(end+1:end+Ncomp)=getfield(get_BLA_cell(PAtype,totNpA),'connections');
% Two-comp principle cell type C model('pCs'=soma, 'pBd'=dendrite)
spec.populations(end+1:end+Ncomp)=getfield(get_BLA_cell(PCtype,totNpC),'populations');
spec.connections(end+1:end+Ncomp)=getfield(get_BLA_cell(PCtype,totNpC),'connections');
% Two-comp interneuron model('ints'=soma, 'intd'=dendrite)
spec.populations(end+1:end+Ncomp)=getfield(get_BLA_cell(Intype,totNbi),'populations');
spec.connections(end+1:end+Ncomp)=getfield(get_BLA_cell(Intype,totNbi),'connections');

for ii=1:numAssemsB
    for jj=1:numAssemsB
        if ii==jj
            assemA=zeros(NpA,NpA);assemC=zeros(NpC,NpC);
            nodigNA=(NpA*NpA-NpA);dumA=rand(nodigNA,1)< p_ee_A;
            mA=0;
            for iiiA=1:NpA
                for jjjA=1:NpA
                    if jjjA ~= iiiA
                        mA=mA+1;
                        assemA(iiiA,jjjA)=dumA(mA);
                    end
                end
            end
            nodigNC=(NpC*NpC-NpC);dumC=rand(nodigNC,1)< p_ee_C ;
            mC=0;
            for iiiC=1:NpC
                for jjjC=1:NpC
                    if iiiC ~= jjjC
                        mC=mC+1;
                        assemC(iiiC,jjjC)=dumC(mC);
                    end
                end
            end
            
            indA=(ii-1)*NpA+(1:NpA);
            indC=(ii-1)*NpC+(1:NpC);
            KAee(indA,indA)=assemA;
            KCee(indC,indC)=assemC;
            
        elseif (ii+jj ==3) || (ii+jj==7)   % when number of assem is 4
            assemA=rand(NpA) < p_ee_Aass;
            assemC=rand(NpC) < p_ee_Cass;
            indAi=(ii-1)*NpA+(1:NpA);indAj=(jj-1)*NpA+(1:NpA);
            indCi=(ii-1)*NpC+(1:NpC);indCj=(jj-1)*NpC+(1:NpC);
            KAee(indAi,indAj)=assemA;
            KCee(indCi,indCj)=assemC;
            
        end
    end
end

for ii=1:numAssemsB
    for jj=1:numAssemsB
        if ii==jj
            assemAC=rand(NpA,NpC) < p_ee_AC;
            assemCA=rand(NpC,NpA) < p_ee_CA;
            indA=(ii-1)*NpA+(1:NpA);
            indC=(jj-1)*NpC+(1:NpC);
            KACee(indA,indC)=assemAC;
            KCAee(indC,indA)=assemCA;
        elseif (ii+jj ==3) || (ii+jj==7)
            assemAC=rand(NpA,NpC) < p_ee_ACass;
            assemCA=rand(NpC,NpA) < p_ee_CAass;
            indA=(ii-1)*NpA+(1:NpA);
            indC=(jj-1)*NpC+(1:NpC);
            KACee(indA,indC)=assemAC;
            KCAee(indC,indA)=assemCA;
        end
    end
end

for ii=1:numAssemsBint
    for jj=1:numAssemsBint
        if ii==jj
            assemI=zeros(Nbi,Nbi);
            nodigNI=(Nbi*Nbi-Nbi);dumI=rand(nodigNI,1)< p_ii;
            mI=0;
            for iiiI=1:Nbi
                for jjjI=1:Nbi
                    if jjjI ~= iiiI
                        mI=mI+1;
                        assemI(iiiI,jjjI)=dumI(mI);
                    end
                end
            end
            indI=(ii-1)*Nbi+(1:Nbi);
            Kbii(indI,indI)=assemI;
            
        elseif mod(ii,2) == 1 && (ii-jj) == -1
            assemI=rand(Nbi) < p_ii_ass;
            indi=(ii-1)*Nbi+(1:Nbi);indj=(jj-1)*Nbi+(1:Nbi);
            Kbii(indi,indj)=assemI;
        elseif mod(ii,2) == 0 && (ii-jj) ==1
            assemI=rand(Nbi) < p_ii_ass;
            indi=(ii-1)*Nbi+(1:Nbi);indj=(jj-1)*Nbi+(1:Nbi);
            Kbii(indi,indj)=assemI;
        end
    end
end
% STDP to interneurons in the BLA
for ii=1:numAssemsB  % numAssemsB=4 for A and C
    for jj=1:numAssemsBint    % numAssemsBint=8 for interneurons
        
        if ii == 1 && (jj ==1 || jj ==6)
            if jj ==1
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei(indA,indAi)=assemAei;KCei(indC,indCi)=assemCei;
            else
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei_stdp(indA,indAi)=assemAei;KCei_stdp(indC,indCi)=assemCei;
            end
        elseif ii == 2 && (jj ==3 || jj ==8)
            if jj ==3
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei(indA,indAi)=assemAei;KCei(indC,indCi)=assemCei;
            else
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei_stdp(indA,indAi)=assemAei;KCei_stdp(indC,indCi)=assemCei;
            end
        elseif ii == 3 && (jj ==2 || jj ==5)
            if jj ==5
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei(indA,indAi)=assemAei;KCei(indC,indCi)=assemCei;
            else
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei_stdp(indA,indAi)=assemAei;KCei_stdp(indC,indCi)=assemCei;
            end
        elseif  ii == 4 && (jj ==4 || jj ==7)
            if jj ==7
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei(indA,indAi)=assemAei;KCei(indC,indCi)=assemCei;
            else
                assemAei=rand(NpA,Nbi) < p_ei_A;
                assemCei=rand(NpC,Nbi) < p_ei_C;
                indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
                indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
                KAei_stdp(indA,indAi)=assemAei;KCei_stdp(indC,indCi)=assemCei;
            end
        end
    end
end

for ii=1:numAssemsB  % numAssemsB=4 for A and C
    for jj=1:numAssemsBint    % numAssemsBint=8 for interneurons
        
        if jj == (2*ii -1)
            assemAie=rand(Nbi,NpA) < p_ie_A;
            assemCie=rand(Nbi,NpC) < p_ie_C;
            indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
            indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
            KAie(indAi,indA)=assemAie;
            KCie(indCi,indC)=assemCie;
        elseif jj == 2*ii
            assemAie=rand(Nbi,NpA) < p_ie_A;
            assemCie=rand(Nbi,NpC) < p_ie_C;
            indA=(ii-1)*NpA+(1:NpA);indAi=(jj-1)*Nbi+(1:Nbi);
            indC=(ii-1)*NpC+(1:NpC);indCi=(jj-1)*Nbi+(1:Nbi);
            KAie(indAi,indA)=assemAie;
            KCie(indCi,indC)=assemCie;
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%BLA-PFC(PL) assemblies network
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for jj=1:numAssemsB-2  %% PL-BLA/C assemblies 1 and 2 connectivity
    for ii=1:numAssems
        if jj==ii  % for diagonal component
            % assemblies BLA_C
            assemPLBC=rand(NPLe,NpC) < p_ee_PLBC;
            assemBCPL=rand(NpC,NPLe) < p_ee_BCPL;
            indPL=(ii-1)*NPLe+(1:NPLe);
            indBC=(jj-1)*NpC+(1:NpC);
            KPLBCee(indPL,indBC)=assemPLBC;
            KBCPLee(indBC,indPL)=assemBCPL;
            
            % assemblies BLA_A
            assemPLBA=rand(NPLe,NpA) < p_ee_PLBA;
            assemBAPL=rand(NpA,NPLe) < p_ee_BAPL;
            indPL=(ii-1)*NPLe+(1:NPLe);
            indBA=(jj-1)*NpA+(1:NpA);
            KPLBAee(indPL,indBA)=assemPLBA;
            KBAPLee(indBA,indPL)=assemBAPL;
            
        else % assemblies BLA_C for off-diagonal component
            % assemblies BLA_C
            assemPLBC=rand(NPLe,NpC) < p_ee_PLBCass;
            assemBCPL=rand(NpC,NPLe) < p_ee_BCPLass;
            indPL=(ii-1)*NPLe+(1:NPLe);
            indBC=(jj-1)*NpC+(1:NpC);
            KPLBCee(indPL,indBC)=assemPLBC;
            KBCPLee(indBC,indPL)=assemBCPL;
            
            % assemblies BLA_A
            assemPLBA=rand(NPLe,NpA) < p_ee_PLBAass;
            assemBAPL=rand(NpA,NPLe) < p_ee_BAPLass;
            indPL=(ii-1)*NPLe+(1:NPLe);
            indBA=(jj-1)*NpA+(1:NpA);
            KPLBAee(indPL,indBA)=assemPLBA;
            KBAPLee(indBA,indPL)=assemBAPL;
        end
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%BLA-PFC(IL) assemblies network
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
for jj=3:numAssemsB
    for ii=1:numAssems
        
        if abs(jj-ii) ==2  % for diagonal component
            % assemblies BLA_C
            assemILBC=rand(NILe,NpC) < p_ee_ILBC;
            assemBCIL=rand(NpC,NILe) < p_ee_BCIL;
            indIL=(ii-1)*NILe+(1:NILe);
            indBC=(jj-1)*NpC+(1:NpC);
            KILBCee(indIL,indBC)=assemILBC;
            KBCILee(indBC,indIL)=assemBCIL;
            
            % assemblies BLA_A
            assemILBA=rand(NILe,NpA) < p_ee_ILBA;
            assemBAIL=rand(NpA,NILe) < p_ee_BAIL;
            indIL=(ii-1)*NILe+(1:NILe);
            indBA=(jj-1)*NpA+(1:NpA);
            KILBAee(indIL,indBA)=assemILBA;
            KBAILee(indBA,indIL)=assemBAIL;
            
        else % assemblies BLA_C for off-diagonal component
            % assemblies BLA_C
            assemILBC=rand(NILe,NpC) < p_ee_ILBCass;
            assemBCIL=rand(NpC,NILe) < p_ee_BCILass;
            indIL=(ii-1)*NILe+(1:NILe);
            indBC=(jj-1)*NpC+(1:NpC);
            KILBCee(indIL,indBC)=assemILBC;
            KBCILee(indBC,indIL)=assemBCIL;
            
            % assemblies BLA_A
            assemILBA=rand(NILe,NpA) < p_ee_ILBAass;
            assemBAIL=rand(NpA,NILe) < p_ee_BAILass;
            indIL=(ii-1)*NILe+(1:NILe);
            indBA=(jj-1)*NpA+(1:NpA);
            KILBAee(indIL,indBA)=assemILBA;
            KBAILee(indBA,indIL)=assemBAIL;
        end
    end
end

for jj=1:numAssemsB
    for ii=1:numAssems
        if jj == ii*2
            assemILBI=rand(NILe,Nbi) < p_ei_ILBI;
            indIL=(ii-1)*NILe+(1:NILe);
            indBI=(jj-1)*Nbi+(1:Nbi);
            KILBIei(indIL,indBI)=assemILBI;
        end
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% NETWORK CONNECTIVITY NORMALIZATION by The TOTAL # of PRESYNAPTIC NEURONs
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% for PL and IL
NKPLee=repmat(max(1,sum(KPLee,1)),[size(KPLee,1) 1]);
NKBCPLee=repmat(max(1,sum(KBCPLee,1)),[size(KBCPLee,1) 1]);
NKILee=repmat(max(1,sum(KILee,1)),[size(KILee,1) 1]);
NKBCILee=repmat(max(1,sum(KBCILee,1)),[size(KBCILee,1) 1]);
% for pA in the BLA
NKAee=repmat(max(1,sum(KAee,1)),[size(KAee,1) 1]);
NKCAee=repmat(max(1,sum(KCAee,1)),[size(KCAee,1) 1]);
% for pC in the BLA
NKCee=repmat(max(1,sum(KCee,1)),[size(KCee,1) 1]);
NKACee=repmat(max(1,sum(KACee,1)),[size(KACee,1) 1]);
NKPLBCee=repmat(max(1,sum(KPLBCee,1)),[size(KPLBCee,1) 1]);
NKILBCee=repmat(max(1,sum(KILBCee,1)),[size(KILBCee,1) 1]);

% for int in the BLA
NKAei=repmat(max(1,sum(KAei,1)),[size(KAei,1) 1]);
NKAei_stdp=repmat(max(1,sum(KAei_stdp,1)),[size(KAei_stdp,1) 1]);
NKCei=repmat(max(1,sum(KCei,1)),[size(KCei,1) 1]);
NKCei_stdp=repmat(max(1,sum(KCei_stdp,1)),[size(KCei_stdp,1) 1]);
NKILBIei=repmat(max(1,sum(KILBIei,1)),[size(KILBIei,1) 1]);

% The total number of each cell
totN2PL=NKPLee(1,:)+NKBCPLee(1,:);% PL <- pC/PL
totN2IL=NKILee(1,:)+NKBCILee(1,:);% IL <- pC/IL
totN2A=NKAee(1,:)+NKCAee(1,:);% pA <- pA/pC
totN2C=NKCee(1,:)+NKACee(1,:)+NKPLBCee(1,:)+NKILBCee(1,:);% pC <- pC/pA/IL/PL
totN2int=NKAei(1,:)+NKAei_stdp(1,:)+NKCei(1,:)+NKCei_stdp(1,:)+NKILBIei(1,:);% int <- pA/pC/IL

%PL
KPLee=KPLee./repmat(max(1,totN2PL),[size(KPLee,1) 1]);
KPLii=KPLii./repmat(max(1,sum(KPLii,1)),[size(KPLii,1) 1]);%Kpii=rand(Npi,Npi)<p_ii;
KPLei=KPLei./repmat(max(1,sum(KPLei,1)),[size(KPLei,1) 1]);
KPLie=KPLie./repmat(max(1,sum(KPLie,1)),[size(KPLie,1) 1]);
% IL
KILee=KILee./repmat(max(1,totN2IL),[size(KILee,1) 1]);
KILii=KILii./repmat(max(1,sum(KILii,1)),[size(KILii,1) 1]);%Kpii=rand(Npi,Npi)<p_ii;
KILei=KILei./repmat(max(1,sum(KILei,1)),[size(KILei,1) 1]);
KILie=KILie./repmat(max(1,sum(KILie,1)),[size(KILie,1) 1]);

% BLA
KAee=KAee./repmat(max(1,totN2A),[size(KAee,1) 1]);%%%%% pA <- pA
KCee=KCee./repmat(max(1,totN2C),[size(KCee,1) 1]);%%%%% pC <- pC

KACee=KACee./repmat(max(1,totN2C),[size(KACee,1) 1]);%%%%%%%%% pC <- pA
KCAee=KCAee./repmat(max(1,totN2A),[size(KCAee,1) 1]);%%%%%%%%% pA <- pC
KAie=KAie./repmat(max(1,sum(KAie,1)),[size(KAie,1) 1]);
KCie=KCie./repmat(max(1,sum(KCie,1)),[size(KCie,1) 1]);

Kbii=Kbii./repmat(max(1,sum(Kbii,1)),[size(Kbii,1) 1]);%Kbii=rand(totNbi,totNbi)<p_ii;
KAei=KAei./repmat(max(1,totN2int),[size(KAei,1) 1]);%%%%%% int <- pA
KAei_stdp=KAei_stdp./repmat(max(1,sum(totN2int,1)),[size(KAei_stdp,1) 1]);%%%%%% int <- pA
KCei=KCei./repmat(max(1,totN2int),[size(KCei,1) 1]);%%%%%% int <- pC
KCei_stdp=KCei_stdp./repmat(max(1,sum(totN2int,1)),[size(KCei_stdp,1) 1]);%%%%%% int <- pC

%PL-BLA
KPLBCee=KPLBCee./repmat(max(1,totN2C),[size(KPLBCee,1) 1]);%%%%%% pC <- PL
KBCPLee=KBCPLee./repmat(max(1,totN2PL),[size(KBCPLee,1) 1]);%%%%%% PL <- pC
KPLBAee=KPLBAee./repmat(max(1,sum(KPLBAee,1)),[size(KPLBAee,1) 1]);
KBAPLee=KBAPLee./repmat(max(1,sum(KBAPLee,1)),[size(KBAPLee,1) 1]);

%IL-BLA
KILBCee=KILBCee./repmat(max(1,totN2C),[size(KILBCee,1) 1]);%%%%%% pC <- IL
KBCILee=KBCILee./repmat(max(1,totN2IL),[size(KBCILee,1) 1]);%%%%%% IL <- pC
KILBAee=KILBAee./repmat(max(1,sum(KILBAee,1)),[size(KILBAee,1) 1]);
KBAILee=KBAILee./repmat(max(1,sum(KBAILee,1)),[size(KBAILee,1) 1]);
KILBIei=KILBIei./repmat(max(1,totN2int),[size(KILBIei,1) 1]);%%%%%% int <- IL
%--------------------------------------------------------------------------
%------------------------PL network----------------------------------------
%--------------------------------------------------------------------------
% [DS02,DG07]: AMPA (taur=.2,taud=1,E=0), NMDA (taur=2.3,taud=95,E=0), GABA (taur=.5,taud=5,E=-75)
tauAMPAr=.2;  % ms, AMPA rise time
tauAMPAd=1;   % ms, AMPA decay time
tauNMDAr=2.3; % ms, NMDA rise time
tauNMDAd=95;  % ms, NMDA decay time
tauGABAr=.5;  % ms, GABAa rise time
tauGABAd=5;   % ms, GABAa decay time

onoffP=1;%<---------------Turn on/off PLnetwork----------------------------
gAMPA=onoffP*3e-3;%3     % uS, PY->PY, maximal AMPA conductance
gNMDA=gAMPA/2; % uS, PY->PY, maximal NMDA conductance
gGABA=onoffP*.5e-3;%.2      % uS, FS->PY, maximal GABAa conductance

% recurrent connections between pyramidal cells (add to existing intercompartmental connections)
index=find(strcmp('EPLs->EPLd',{spec.connections.direction}),1,'first');
spec.connections(index).mechanism_list=[{'iAMPA'},{'iNMDA'},spec.connections(index).mechanism_list(:)'];
spec.connections(index).parameters=[{'gAMPA'},{gAMPA},{'gNMDA'},{gNMDA},{'EAMPA'},{0},{'ENMDA'},{0},{'netcon'},{KPLee},...
    {'tauAMPAr'},{tauAMPAr},{'tauAMPA'},{tauAMPAd},{'tauNMDAr'},{tauNMDAr},{'tauNMDA'},{tauNMDAd},spec.connections(index).parameters(:)'];
% pyramidal<->FS connections
spec.connections(end+1).direction='EPLs->FSPL';
spec.connections(end).mechanism_list={'iAMPA','iNMDA'};
spec.connections(end).parameters={'gAMPA',gAMPA,'gNMDA',gNMDA,'netcon',KPLei,'EAMPA',0,'ENMDA',0,...
    'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd,'tauNMDAr',tauNMDAr,'tauNMDA',tauNMDAd};
spec.connections(end+1).direction='FSPL->EPLs';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABA,'netcon',KPLie,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};

% interneuron<->interneuron connections
spec.connections(end+1).direction='FSPL->FSPL';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABA,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75,'netcon',KPLii};
%--------------------------------------------------------------------------
%------------------------IL network----------------------------------------
%--------------------------------------------------------------------------
% [DS02,DG07]: AMPA (taur=.2,taud=1,E=0), NMDA (taur=2.3,taud=95,E=0), GABA (taur=.5,taud=5,E=-75)
tauAMPAr=.2;  % ms, AMPA rise time
tauAMPAd=1;   % ms, AMPA decay time
tauNMDAr=2.3; % ms, NMDA rise time
tauNMDAd=95;  % ms, NMDA decay time
tauGABAr=.5;  % ms, GABAa rise time
tauGABAd=5;   % ms, GABAa decay time

onoffL=1;%<---------------Turn on/off ILnetwork----------------------------
gAMPA=onoffL*3e-3;%3     % uS, PY->PY, maximal AMPA conductance
gNMDA=gAMPA/2; % uS, PY->PY, maximal NMDA conductance
gGABA=onoffL*.5e-3;%.2      % uS, FS->PY, maximal GABAa conductance
% recurrent connections between pyramidal cells (add to existing intercompartmental connections)
index=find(strcmp('EILs->EILd',{spec.connections.direction}),1,'first');
spec.connections(index).mechanism_list=[{'iAMPA'},{'iNMDA'},spec.connections(index).mechanism_list(:)'];
spec.connections(index).parameters=[{'gAMPA'},{gAMPA},{'gNMDA'},{gNMDA},{'EAMPA'},{0},{'ENMDA'},{0},{'netcon'},{KILee},...
    {'tauAMPAr'},{tauAMPAr},{'tauAMPA'},{tauAMPAd},{'tauNMDAr'},{tauNMDAr},{'tauNMDA'},{tauNMDAd},spec.connections(index).parameters(:)'];
% pyramidal<->FS connections
spec.connections(end+1).direction='EILs->FSIL';
spec.connections(end).mechanism_list={'iAMPA','iNMDA'};
spec.connections(end).parameters={'gAMPA',gAMPA,'gNMDA',gNMDA,'netcon',KILei,'EAMPA',0,'ENMDA',0,...
    'tauAMPAr',tauAMPAr,'tauAMPA',tauAMPAd,'tauNMDAr',tauNMDAr,'tauNMDA',tauNMDAd};
spec.connections(end+1).direction='FSIL->EILs';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABA,'netcon',KILie,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75};

% interneuron<->interneuron connections
spec.connections(end+1).direction='FSIL->FSIL';
spec.connections(end).mechanism_list={'iGABA'};
spec.connections(end).parameters={'gGABA',gGABA,'tauGABAr',tauGABAr,'tauGABA',tauGABAd,'EGABA',-75,'netcon',KILii};
spectmpBLA=spec;
%--------------------------------------------------------------------------
%------------------------BLA network---------------------------------------
%--------------------------------------------------------------------------------

onoffB=1;%<---------------Turn on/off BLAnetwork----------------------------
gAmpa = .1*onoffB; gNmda = .05*onoffB;gGaba = .5*onoffB;%mS/cm2

% recurrent connections between pyramidal cells (add to existing intercompartmental connections)
index=find(strcmp('pAs->pAd',{spec.connections.direction}),1,'first');
spec.connections(index).mechanism_list=[{'iAMPAamy'},{'iNMDAamy'},spec.connections(index).mechanism_list(:)'];
spec.connections(index).parameters=[{'gAmpa'},{gAmpa},{'gNmda'},{gNmda},{'netcon'},{KAee},{'eampa'},{0},{'enmda'},{0},...
    spec.connections(index).parameters(:)'];%{'Atau1'},{Atau1},{'Atau2'},{Atau2},{'Ntau1'},{Ntau1},{'Ntau2'},{Ntau2},

index=find(strcmp('pCs->pCd',{spec.connections.direction}),1,'first');
spec.connections(index).mechanism_list=[{'iAMPAamy'},{'iNMDAamy'},spec.connections(index).mechanism_list(:)'];
spec.connections(index).parameters=[{'gAmpa'},{gAmpa},{'gNmda'},{gNmda},{'netcon'},{KCee},{'eampa'},{0},{'enmda'},{0},...
    spec.connections(index).parameters(:)'];%{'Atau1'},{Atau1},{'Atau2'},{Atau2},{'Ntau1'},{Ntau1},{'Ntau2'},{Ntau2},

% % pry A <-> pry C connections
spec.connections(end+1).direction='pCs->pAd';
spec.connections(end).mechanism_list={'iAMPAamy','iNMDAamy'};
spec.connections(end).parameters={'gAmpa',gAmpa,'gNmda',gNmda,'netcon',KCAee,'eampa',0,'enmda',0};
%'Atau1',Atau1,'Atau2',Atau2,'Ntau1',Ntau1,'Ntau2',Ntau2

spec.connections(end+1).direction='pAs->pCd';
spec.connections(end).mechanism_list={'iAMPAamy','iNMDAamy'};
spec.connections(end).parameters={'gAmpa',gAmpa,'gNmda',gNmda,'netcon',KACee,'eampa',0,'enmda',0};
%'Atau1',Atau1,'Atau2',Atau2,'Ntau1',Ntau1,'Ntau2',Ntau2

% pry A <-> int connections
spec.connections(end+1).direction='ints->pAs';
spec.connections(end).mechanism_list={'iGABAamy'};
spec.connections(end).parameters={'gGaba',gGaba,'egaba',-75,'netcon',KAie};%,'Gtau1',Gtau1,'Gtau2',Gtau2

Atau1 = 0.3;% 0.4J;	 %rise time constant: 0.5for pry, 0.3for interneuron
Atau2 = 2.4;     %decay time constant: 7ms for pyramidal, 2.4for interneuron
spec.connections(end+1).direction='pAs->intd';
spec.connections(end).mechanism_list={'iAMPAamy','iNMDAamy'};%,
spec.connections(end).parameters={'gAmpa',gAmpa,'gNmda',gNmda,'netcon',KAei,'eampa',0,'enmda',0,...
    'Atau1',Atau1,'Atau2',Atau2};%,'Ntau1',Ntau1,'Ntau2',Ntau2

% recurrent connections between int <-> int
spec.connections(end+1).direction='ints->ints';
spec.connections(end).mechanism_list={'iGABAamy'};
spec.connections(end).parameters={'gGaba',0.1,'egaba',-75,'netcon',Kbii};

% %pry C <-> int connections
spec.connections(end+1).direction='ints->pCs';
spec.connections(end).mechanism_list={'iGABAamy'};
spec.connections(end).parameters={'gGaba',gGaba,'egaba',-75,'netcon',KCie};

spec.connections(end+1).direction='pCs->intd';
spec.connections(end).mechanism_list={'iAMPAamy','iNMDAamy'};
spec.connections(end).parameters={'gAmpa',gAmpa,'gNmda',gNmda,'netcon',KCei,'eampa',0,'enmda',0,...
    'Atau1',Atau1,'Atau2',Atau2};
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
spectmp=spec;
spectmpORI=spec;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
offW1=20000;%onW1=0;
spec.connections(end+1).direction='pAs->intd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%,
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa,'gNmda',gNmda,'netcon',KAei_stdp,...
    'Atau1',Atau1,'Atau2',Atau2};

spec.connections(end+1).direction='pCs->intd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa,'gNmda',gNmda,'netcon',KCei_stdp,...
    'Atau1',Atau1,'Atau2',Atau2};
%--------------------------------------------------------------------------
%------------------------ PL-BLA network-----------------------------------
%--------------------------------------------------------------------------
% STDP synaptic conductances
% pyramidal soma (length and diameter chosen to match surface area and internal resistance with a sphere of diam 23 microns)
ls=28.618;  % um, length
ds=21.840;  % um, diameter
% pyramidal dendrite
ld=650;     % um, length
dd=6.5;     % um, diameter
As=ds*ls*pi;Ad=dd*ld*pi;

scl_STDP=1;scl_PB=2;% 1 2
parP=-10;scl_REM=(100+parP)/100;%--------REM conditions for pA/pC
gAmpa_STDP = .1*scl_STDP; gNmda_STDP = .05*scl_STDP;%gGaba_STDP = .5*scl_STDP;
% Assemble connection from PFC to BLA_C
spec.connections(end+1).direction='EPLs->pCd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*scl_PB,'gNmda',gNmda_STDP*scl_REM*scl_PB,'netcon',KPLBCee};

% Assemble connection from BLA_C to PFC
spec.connections(end+1).direction='pCs->EPLd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*Ad*1e-5,'gNmda',gNmda_STDP*Ad*1e-5,'netcon',KBCPLee};

% Assemble connection from PFC to BLA_A
spec.connections(end+1).direction='EPLs->pAd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*scl_PB,'gNmda',gNmda_STDP*scl_REM*scl_PB,'netcon',KPLBAee};

% Assemble connection from BLA_A to PFC
spec.connections(end+1).direction='pAs->EPLd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*Ad*1e-5,'gNmda',gNmda_STDP*Ad*1e-5,'netcon',KBAPLee};

%--------------------------------------------------------------------------
%------------------------ IL-BLA network-----------------------------------
%--------------------------------------------------------------------------
% STDP synaptic conductances
% pyramidal soma (length and diameter chosen to match surface area and internal resistance with a sphere of diam 23 microns)
ls=28.618;  % um, length
ds=21.840;  % um, diameter
% pyramidal dendrite
ld=650;     % um, length
dd=6.5;     % um, diameter
As=ds*ls*pi;Ad=dd*ld*pi;

scl_STDP=1;scl_PB=2; % 1 2
parP=-10;scl_REM=(100+parP)/100;%--------REM conditions for pA/pC
gAmpa_STDP = .1*scl_STDP; gNmda_STDP = .05*scl_STDP;%gGaba_STDP = .5*scl_STDP;
% Assemble connection from PFC to BLA_C
spec.connections(end+1).direction='EILs->pCd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*scl_PB,'gNmda',gNmda_STDP*scl_REM*scl_PB,'netcon',KILBCee};

% Assemble connection from BLA_C to PFC
spec.connections(end+1).direction='pCs->EILd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*Ad*1e-5,'gNmda',gNmda_STDP*Ad*1e-5,'netcon',KBCILee};

% Assemble connection from PFC to BLA_A
spec.connections(end+1).direction='EILs->pAd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*scl_PB,'gNmda',gNmda_STDP*scl_REM*scl_PB,'netcon',KILBAee};

% Assemble connection from BLA_A to PFC
spec.connections(end+1).direction='pAs->EILd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*Ad*1e-5,'gNmda',gNmda_STDP*Ad*1e-5,'netcon',KBAILee};

% STDP on interneurons from IL
spec.connections(end+1).direction='EILs->intd';
spec.connections(end).mechanism_list={'STDPupW0_Isyn'};%iAMPA_Hebb .STDP_Isyn
spec.connections(end).parameters={'offW1',offW1,'gAmpa',gAmpa_STDP*scl_PB,'gNmda',gNmda_STDP*scl_REM*scl_PB,'netcon',KILBIei,...
    'Atau1',Atau1,'Atau2',Atau2};

%--------------------------------------------------------------------------
% MAIN SIMULATION
%--------------------------------------------------------------------------
% soma for FS
l=27; d=29; Afs=d*l*pi; % um2, cylinder surface area without the ends
% pyramidal soma (length and diameter chosen to match surface area and internal resistance with a sphere of diam 23 microns)
ls=28.618;ds=21.840;   % um, length,  um, diameter
% pyramidal dendrite
ld=650;  dd=6.5;   % um, length, um, diameter
As=ds*ls*pi;Ad=dd*ld*pi;

noise_ampB=1.5;iexc_avg=0;
noise_amp=noise_ampB*Ad*1e-5; noise_ampFS=noise_ampB*Afs*1e-5;
sin_ampB=0;sin_ampPLon=0;

sin_scl=2;
sin_ampILon=0;%<-------------- sharp peak having a specific frequency
sin_ampPL=sin_ampPLon*Ad*1e-5; sin_ampPLFS=sin_ampPLon*Afs*1e-5;
sin_ampIL=sin_ampILon;sin_ampILFS=sin_ampILon;
sin_NormIL=7.2260e+86/(sin_scl*Ad*1e-5);sin_NormFS=7.2260e+86/(sin_scl*Afs*1e-5);
dum_fB=0;dum_f=4;

targetAssemPL=0;targetAssemIL=4;
targetAssemB=0;targetAssemBint=0;

% Ed:idc=.485:5hz(~7hz w/netw), .6:20hz
vary={'EPLd','Iapp',0.4705;'EPLd','noise_amp',noise_amp;'EPLd','iexc_avg',iexc_avg;'EPLd','sin_amp',sin_ampPL;'EPLd','mf',dum_f;...
    'EPLd','numAssems',numAssems; 'EPLs','numAssems',numAssems;'EPLd','targetAssem',targetAssemPL;...
    'FSPL','Iapp',0.0505;'FSPL','noise_amp',noise_ampFS;'FSPL','iexc_avg',iexc_avg;'FSPL','sin_amp',sin_ampPLFS;'FSPL','mf',dum_f;...
    'EILd','Iapp',0.2405;'EILd','noise_amp',noise_amp;'EILd','iexc_avg',iexc_avg;'EILd','sin_amp',sin_ampIL;'EILd','mf',dum_f;...
    'EILd','numAssems',numAssems;'EILs','numAssems',numAssems;'EILd','targetAssem',targetAssemIL;...
    'FSIL','Iapp',0.0504;'FSIL','noise_amp',noise_ampFS;'FSIL','iexc_avg',iexc_avg;'FSIL','sin_amp',sin_ampILFS;'FSIL','mf',dum_f;...
    'pAd','Iapp',13;'pAd','noise_amp',noise_ampB;'pAd','iexc_avg',iexc_avg;'pAd','sin_amp',sin_ampB;'pAd','mf',dum_fB;
    'pAs','numAssems',numAssemsB;'pAd','numAssems',numAssemsB;'pAd','targetAssem',targetAssemB;...
    'pCd','Iapp',10.8;'pCd','noise_amp',noise_ampB;'pCd','iexc_avg',iexc_avg;'pCd','sin_amp',sin_ampB;'pCd','mf',dum_fB;
    'pCs','numAssems',numAssemsB;'pCd','numAssems',numAssemsB;'pCd','targetAssem',targetAssemB;...
    % intd Iapp was 3.32 replaced by 3 in gGie=0.8, gGii=0.2
    'intd','Iapp',3.;'intd','noise_amp',noise_ampB;'intd','iexc_avg',iexc_avg;'intd','sin_amp',sin_ampB;'intd','mf',dum_fB;
    'ints','numAssems',numAssemsBint;'intd','numAssems',numAssemsBint;'intd','targetAssem',targetAssemBint};

tend=2;
tspan=[0 tend];dt=.01;  % [beg end], ms
solver_options={'tspan',tspan,'solver','euler','dt',dt,'matCompatibility_flg',0,'verbose_flag',1,'parfor_flag',0};

data=dsSimulate(spec,'vary',vary,solver_options{:});

Wdum={'pCd_EPLs_STDPupW0_Isyn_w','EPLd_pCs_STDPupW0_Isyn_w','pAd_EPLs_STDPupW0_Isyn_w',...
    'EPLd_pAs_STDPupW0_Isyn_w','pCd_EILs_STDPupW0_Isyn_w','EILd_pCs_STDPupW0_Isyn_w',...
    'pAd_EILs_STDPupW0_Isyn_w','EILd_pAs_STDPupW0_Isyn_w',...
    'intd_pAs_STDPupW0_Isyn_w','intd_pCs_STDPupW0_Isyn_w','intd_EILs_STDPupW0_Isyn_w'};

ic = [];% only iniW reset
state_variables = data.model.state_variables;
for i = 1:length(state_variables)
    x = data.(state_variables{i});
    ss=find(strcmp(Wdum,state_variables{i}));
    if ndims(x) > 2
        if ~isempty(ss)
            ic.(state_variables{i}) =squeeze(x(1,:,:));
        else
            ic.(state_variables{i}) = squeeze(x(end,:,:));
        end
    else
        if ~isempty(ss)
            ic.(state_variables{i}) = x(1,:);
        else
            ic.(state_variables{i}) = x(end,:);
        end
    end
end

tend=10;
tspan=[0 tend];dt=.01;  % [beg end], ms
solver_options={'tspan',tspan,'solver','euler','dt',dt,'matCompatibility_flg',0,'verbose_flag',1,'parfor_flag',0};
fdatnam=['BLAPLIL_REM_10sec']; % PFC(ILPL) and BLA connected

data=dsSimulate(spec,'ic',ic,'vary',vary,solver_options{:});% PFC(ILPL) and BLA connected

data_pAd_v=data.pAd_V;data_pCd_v=data.pCd_V;data_intd_v=data.intd_V;data_time=data.time;
data_EPLd_v=data.EPLd_V;data_FSPL_v=data.FSPL_V;data_EILd_v=data.EILd_V;data_FSIL_v=data.FSIL_V;

EPLd_pAs_STDP_Isyn_w=data.EPLd_pAs_STDPupW0_Isyn_w;EPLd_pCs_STDP_Isyn_w=data.EPLd_pCs_STDPupW0_Isyn_w;
pAd_EPLs_STDP_Isyn_w=data.pAd_EPLs_STDPupW0_Isyn_w;pCd_EPLs_STDP_Isyn_w=data.pCd_EPLs_STDPupW0_Isyn_w;

EILd_pAs_STDP_Isyn_w=data.EILd_pAs_STDPupW0_Isyn_w;EILd_pCs_STDP_Isyn_w=data.EILd_pCs_STDPupW0_Isyn_w;
pAd_EILs_STDP_Isyn_w=data.pAd_EILs_STDPupW0_Isyn_w;pCd_EILs_STDP_Isyn_w=data.pCd_EILs_STDPupW0_Isyn_w;

intd_EILs_STDP_Isyn_w=data.intd_EILs_STDPupW0_Isyn_w;
intd_pAs_STDP_Isyn_w=data.intd_pAs_STDPupW0_Isyn_w;intd_pCs_STDP_Isyn_w=data.intd_pCs_STDPupW0_Isyn_w;

ic1 = [];
state_variables = data.model.state_variables;
for i = 1:length(state_variables)
    x = data.(state_variables{i});
    ss=find(strcmp(Wdum,state_variables{i}));
    if ndims(x) > 2
        if ~isempty(ss)
            ic1.(state_variables{i}) =squeeze(x(end,:,:));
        else
            ic1.(state_variables{i}) = squeeze(x(end,:,:));
        end
    else
        if ~isempty(ss)
            ic1.(state_variables{i}) = x(end,:);
        else
            ic1.(state_variables{i}) = x(end,:);
        end
    end
end
save(fdatnam,'data_pAd_v','data_pCd_v','data_intd_v','data_time','data_EPLd_v','data_FSPL_v', 'data_EILd_v','data_FSIL_v',...
    'EPLd_pAs_STDP_Isyn_w','EPLd_pCs_STDP_Isyn_w','pAd_EPLs_STDP_Isyn_w','pCd_EPLs_STDP_Isyn_w',...
    'EILd_pAs_STDP_Isyn_w','EILd_pCs_STDP_Isyn_w','pAd_EILs_STDP_Isyn_w','pCd_EILs_STDP_Isyn_w',...
    'intd_EILs_STDP_Isyn_w','intd_pAs_STDP_Isyn_w','intd_pCs_STDP_Isyn_w',...
    'KACee','KAee','KAei','KAie','KBAILee','KBAPLee','KBCILee','KBCPLee','Kbii','KCAee','KCee','KCei','KCie',...
    'KILBAee','KILBCee','KILBIei','KILee','KILei','KILie','KILii','sin_ampB','sin_ampPLon','sin_ampILon',...
    'KPLBAee','KPLBCee','KPLee','KPLei','KPLie','KPLii','vary','spec','spectmp','spectmpBLA','spectmpORI',...
    'numAssems','numAssemsB','numAssemsBint','NPLe','NPLi','NILe','NILi','NpA','NpC','Nbi','KAei_stdp','KCei_stdp',...
    'solver_options','Wdum','ic','ic1','-v7.3');

