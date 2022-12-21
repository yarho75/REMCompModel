
function spec = get_BLA_cell(type,N)
% To retrieve DynaSim specification for single cell LA models.
% Three types of Principle neurons w 2compartments: A, B and C
% One type of Interneuron with 2 compartments

if nargin<1, type='pyrC'; end
if nargin<2, N=1; end % # of cells in population

switch type
    case 'pyrA'
        %-----------------------------------------------
        % Type A:
        %-----------------------------------------------

        %         N=1;
        %shared parameters: gl(.00034:old value used) Updated(6/4/2019)
        ena=45;ek=-80;el=-75;eca=120;eh=-43;cm=1;dptau=1.5;dqtau=569;
        scl=10;Far=96485;Vol=1;carest=.00024;gl=0.034;ahptau=48;
        f1=0.7;tauca1=1;f2=0.024;
        % 0.7 for kca, 0.024 for sAHP % 1 for kca, 1000 for sAHP [ms]

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaamy','iKamy','ileakamy','iMamy','iDamy','ihCaamy','CaDyn2','isAHPamy','CaDyn1','iKCaamy','iHamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        gna=120;gk=12;gca=0.1;gm=0.3;
        gkca=0.0;gahp=0.0;gd=0;gh=0.0;
        sin_amp=0;mf=0;tauca2=1000;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='pAs';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=40;gk=3;gca=0.2;gd=1;
        gkca=0.5;gahp=0.1;gm=0.3;gh=0.1;
        sin_amp=0;mf=0;tauca2=1000;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='pAd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};

        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'pAs' 'pAd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
    case 'pyrAs'
        spec = get_BLA_cell('pyrA',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'pyrAd'
        spec = get_BLA_cell('pyrA',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
    case 'pystneA'
        %-----------------------------------------------
        % Type A w/Neuromoulator(nm) effects
        % (Active ST + NE during wake)
        %-----------------------------------------------
        %         N=1;
        %shared parameters: gl(.00034:old value used) Updated(6/4/2019)
        ena=45;ek=-80;el=-75;eca=120;eh=-43;cm=1;dptau=1.5;dqtau=569;
        scl=10;Far=96485;Vol=1;carest=.00024;gl=0.034;ahptau=48;
        f1=0.7;tauca1=1;f2=0.024;
        % 0.7 for kca, 0.024 for sAHP % 1 for kca, 1000 for sAHP [ms]

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaamy','iKamy','ileakamy','iMamy','iDamy','ihCaamy','CaDyn2','isAHPamy','CaDyn1','iKCaamy','iHamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        per=-20;
        ratio=(100+per)/100;%------------------------------------------------------------------------------------
        %--------------------------------------------------------------------------------------------------------------------------
        gna=120;gk=12*ratio;gca=0.1;gm=0.3;
        gkca=0.0;gahp=0.0*ratio;gd=0;gh=0.0;
        sin_amp=0;mf=0;tauca2=1000;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='pAs';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=40;gk=3*ratio;gca=0.2;gd=1;
        gkca=0.5;gahp=0.1*ratio;gm=0.3;gh=0.1;
        sin_amp=0;mf=0;tauca2=1000;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='pAd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'pAs' 'pAd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
    case 'pystneAs'
        spec = get_BLA_cell('pystneA',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'pystneAd'
        spec = get_BLA_cell('pystneA',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
    case 'pyhneA'
        %-----------------------------------------------
        % Type A w/Neuromoulator(nm) effects
        % (Highest NE during REM for PTSD)
        %-----------------------------------------------
        %         N=1;
        %shared parameters: gl(.00034:old value used) Updated(6/4/2019)
        ena=45;ek=-80;el=-75;eca=120;eh=-43;cm=1;dptau=1.5;dqtau=569;
        scl=10;Far=96485;Vol=1;carest=.00024;gl=0.034;ahptau=48;
        f1=0.7;tauca1=1;f2=0.024;
        % 0.7 for kca, 0.024 for sAHP % 1 for kca, 1000 for sAHP [ms]

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaamy','iKamy','ileakamy','iMamy','iDamy','ihCaamy','CaDyn2','isAHPamy','CaDyn1','iKCaamy','iHamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        per=-30;
        ratio=(100+per)/100;%------------------------------------------------------------------------------------
        %--------------------------------------------------------------------------------------------------------------------------
        gna=120;gk=12*ratio;gca=0.1;gm=0.3;
        gkca=0.0;gahp=0.0*ratio;gd=0;gh=0.0;
        sin_amp=0;mf=0;tauca2=1000;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='pAs';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=40;gk=3*ratio;gca=0.2;gd=1;
        gkca=0.5;gahp=0.1*ratio;gm=0.3;gh=0.1;
        sin_amp=0;mf=0;tauca2=1000;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='pAd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'pAs' 'pAd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
    case 'pyhneAs'
        spec = get_BLA_cell('pyhneA',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'pyhneAd'
        spec = get_BLA_cell('pyhneA',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
        
    case 'pyrC'
        %-----------------------------------------------
        % Type C:
        %-----------------------------------------------
        %         N=1;
        %shared parameters
        ena=45;ek=-80;el=-75;eca=120;eh=-43;cm=1;dptau=1.5;dqtau=569;eh=-43;
        scl=10;Far=96485;Vol=1;carest=.00024;gl=0.034;ahptau=48;
        f1=0.7;tauca1=1;f2=0.024;
        % 0.7 for kca, 0.024 for sAHP % 1 for kca, 1000 for sAHP [ms]

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaamy','iKamy','ileakamy','iMamy','iDamy','ihCaamy','CaDyn2','isAHPamy','CaDyn1','iKCaamy','iHamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        gna=120;gk=12;gca=0.1;gm=0.25;
        gkca=0.0;gahp=0.0;gd=0;gh=0;
        sin_amp=0;mf=0;tauca2=120;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='pCs';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};

        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=40;gk=3;gca=0.2;gd=0.1;
        gkca=0.5;gahp=0.5;gm=0.25;gh=0.1;
        sin_amp=0;mf=0;tauca2=120;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='pCd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'pCs' 'pCd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end

    case 'pyrCs'
        spec = get_BLA_cell('pyrC',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'pyrCd'
        spec = get_BLA_cell('pyrC',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
    case 'pystneC'
        %-----------------------------------------------
        % Type C w/Neuromoulator(nm) effects
        % (Active ST + NE during wake)
        %-----------------------------------------------
        %         N=1;
        %shared parameters
        ena=45;ek=-80;el=-75;eca=120;eh=-43;cm=1;dptau=1.5;dqtau=569;eh=-43;
        scl=10;Far=96485;Vol=1;carest=.00024;gl=0.034;ahptau=48;
        f1=0.7;tauca1=1;f2=0.024;
        % 0.7 for kca, 0.024 for sAHP % 1 for kca, 1000 for sAHP [ms]

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaamy','iKamy','ileakamy','iMamy','iDamy','ihCaamy','CaDyn2','isAHPamy','CaDyn1','iKCaamy','iHamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        per=-20;
        ratio=(100+per)/100;%------------------------------------------------------------------------------------
        %--------------------------------------------------------------------------------------------------------------------------
        gna=120;gk=12*ratio;gca=0.1;gm=0.25;
        gkca=0.0;gahp=0.0*ratio;gd=0;gh=0;
        sin_amp=0;mf=0;tauca2=120;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='pCs';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=40;gk=3*ratio;gca=0.2;gd=0.1;
        gkca=0.5;gahp=0.5*ratio;gm=0.25;gh=0.1;
        sin_amp=0;mf=0;tauca2=120;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='pCd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'pCs' 'pCd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
    case 'pystneCs'
        spec = get_BLA_cell('pystneC',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'pystneCd'
        spec = get_BLA_cell('pystneC',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
    case 'pyhneC'
        %-----------------------------------------------
        % Type C  w/Neuromoulator(nm) effects
        % (Highest NE during REM for PTSD)
        %-----------------------------------------------
        %         N=1;
        %shared parameters
        ena=45;ek=-80;el=-75;eca=120;eh=-43;cm=1;dptau=1.5;dqtau=569;eh=-43;
        scl=10;Far=96485;Vol=1;carest=.00024;gl=0.034;ahptau=48;
        f1=0.7;tauca1=1;f2=0.024;
        % 0.7 for kca, 0.024 for sAHP % 1 for kca, 1000 for sAHP [ms]

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaamy','iKamy','ileakamy','iMamy','iDamy','ihCaamy','CaDyn2','isAHPamy','CaDyn1','iKCaamy','iHamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        per=-30;
        ratio=(100+per)/100;%------------------------------------------------------------------------------------
        %--------------------------------------------------------------------------------------------------------------------------
        gna=120;gk=12*ratio;gca=0.1;gm=0.25;
        gkca=0.0;gahp=0.0*ratio;gd=0;gh=0;
        sin_amp=0;mf=0;tauca2=120;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='pCs';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=40;gk=3*ratio;gca=0.2;gd=0.1;
        gkca=0.5;gahp=0.5*ratio;gm=0.25;gh=0.1;
        sin_amp=0;mf=0;tauca2=120;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='pCd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gm',gm,'gd',gd,...
            'gh',gh,'gahp',gahp,'gkca',gkca,'tauca2',tauca2,'mf',mf','sin_amp',sin_amp,'noise_amp',noise_amp,...
            'numAssems',numAssems,'targetAssem',targetAssem,'iexc_avg',iexc_avg};
        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'pCs' 'pCd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
    case 'pyhneCs'
        spec = get_BLA_cell('pyhneC',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'pyhneCd'
        spec = get_BLA_cell('pyhneC',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
   
    case 'intG'
        %-----------------------------------------------
        % Inter neuron: int
        %-----------------------------------------------
        %shared parameters
        ena=45;ek=-80;el=-70;cm=1;gl=0.05;gSYN=0;

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaint','iKint','ileakamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        gna=35;gk=8;
        sin_amp=0;mf=0;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='ints';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'mf',mf',...
            'numAssems',numAssems,'targetAssem',targetAssem,'sin_amp',sin_amp,'noise_amp',noise_amp,'iexc_avg',iexc_avg};

        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=10;gk=3;
        sin_amp=0;mf=0;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='intd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'mf',mf',...
            'numAssems',numAssems,'targetAssem',targetAssem,'sin_amp',sin_amp,'noise_amp',noise_amp,'iexc_avg',iexc_avg};

        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'ints' 'intd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
        %        Indata=dsSimulate(spec,'tspan',[0 5000],'vary',{'ints','Iapp',1.5;'intd','Iapp',1.5},'solver','rk1','verbose_flag',1);
        %        dsPlot(Indata,'xlim',[0 500],'ylim',[-100 50]);

    case 'intGs'
        spec = get_PFC_cell('intG',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'intGd'
        spec = get_PFC_cell('intG',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
    case 'intstneG'
        %-----------------------------------------------
        % Interneuron w/Neuromoulator(nm) effects
        % (Active ST + NE during wake)
        %-----------------------------------------------
        %shared parameters
        per=-20;
        ratio=(100+per)/100;%------------------------------------------------------------------------------------
        %--------------------------------------------------------------------------------------------------------
        ena=45;ek=-80;el=-70;cm=1;gl=0.05*ratio;gSYN=0;

        % ionic mechanisms and voltage dynamics present in both compartments (see [DS00] Methods for justification)
        mechanism_list={'iNaint','iKint','ileakamy','isin_noise'};
        state_equations={'dV/dt=(@current+Iapp)./cm; Iapp=0; cm=1; V(0)=-70'
            'monitor V.spikes(0)'};
        spec=[];

        % soma
        gna=35;gk=8;
        sin_amp=0;mf=0;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(1).name='ints';
        spec.populations(1).size=N;
        spec.populations(1).equations=state_equations;
        spec.populations(1).mechanism_list=mechanism_list;
        spec.populations(1).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'mf',mf',...
            'numAssems',numAssems,'targetAssem',targetAssem,'sin_amp',sin_amp,'noise_amp',noise_amp,'iexc_avg',iexc_avg};
        % spec.populations(2).parameters={'Iapp',0,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'gca',gca,'gkca',gkca,'gahp',gahp};
        % dend
        gna=10;gk=3;
        sin_amp=0;mf=0;
        noise_amp=0;iexc_avg=0;numAssems=0;targetAssem=0;
        spec.populations(2).name='intd';
        spec.populations(2).size=N;
        spec.populations(2).equations=state_equations;
        spec.populations(2).mechanism_list=mechanism_list;
        spec.populations(2).parameters={'Iapp',0.,'cm',cm,'gna',gna,'gk',gk,'gl',gl,'mf',mf',...
            'numAssems',numAssems,'targetAssem',targetAssem,'sin_amp',sin_amp,'noise_amp',noise_amp,'iexc_avg',iexc_avg};

        %intercompartmental connections
        Ri=150;% Ohm-cm
        compartments={'ints' 'intd'};
        connections={[1 2],[2 1]};
        for c=1:length(connections)
            src=connections{c}(1);
            dst=connections{c}(2);
            spec.connections(c).direction=[compartments{src} '->' compartments{dst}];
            spec.connections(c).mechanism_list={'icomp'};
            gc=0.3;% gc=1/mean(Ri*4*lengths./(pi*diameters.^2));
            spec.connections(c).parameters={'gc',gc};
        end
        %        Indata=dsSimulate(spec,'tspan',[0 5000],'vary',{'ints','Iapp',1.5;'intd','Iapp',1.5},'solver','rk1','verbose_flag',1);
        %        dsPlot(Indata,'xlim',[0 500],'ylim',[-100 50]);

    case 'intstneGs'
        spec = get_PFC_cell('intstneG',N);
        spec.populations=spec.populations(1);
        spec.connections=[];
    case 'intstneGd'
        spec = get_PFC_cell('intstneG',N);
        spec.populations=spec.populations(2);
        spec.connections=[];
end

