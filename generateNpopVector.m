
function upNpop = generateNpopVector(totNpop,numAssems,numtarget,targetAssem1,targetAssem2)
iniKK=zeros(1,totNpop);
if numtarget ==1
    targetAssem=targetAssem1;
    if targetAssem ~= 0
        ncell=totNpop/numAssems;
        ind=((targetAssem-1)*ncell+1):ncell*targetAssem;
        iniKK(1,ind)=1;
        upNpop=iniKK;
    elseif numAssems ==0
        iniKK=ones(1,totNpop);
        upNpop=iniKK;
    else
        upNpop=iniKK;
    end
else
    if targetAssem1 ~= 0
        ncell=totNpop/numAssems;
        ind1=((targetAssem1-1)*ncell+1):ncell*targetAssem1;
        iniKK(1,ind1)=1;
        ind2=((targetAssem2-1)*ncell+1):ncell*targetAssem2;
        iniKK(1,ind2)=1;
        upNpop=iniKK;
    elseif numAssems ==0
        iniKK=ones(1,totNpop);
        upNpop=iniKK;
    else
        upNpop=iniKK;
    end
end
