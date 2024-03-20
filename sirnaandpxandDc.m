clc
close all
clear all
%%parametrs
gamaa=0.000038; %%γ %%Effect of the vaccine on adenosine
effecthypoxionadenosine=0.9412;  %%a %%Effect of hypoxia on adenosine
r=0.07;  %%r %%Adenosine elimination rate
y1=12; %%β  %%Constant rate of adenosine production
effectofadenosineonTumor=0.000025;  %%rAt  %%Effect of adenosine on tumor
effectofhypoxiaonTumor=0.0021;   %%rHt   %%Effect of hypoxia on tumor
effectadenoonTreg=0.0015;     %%rAr  %%Adenosine effect on Treg
effecthypoxiaonTreg=0.0047;  %%rHr   %%Hypoxia effect on Treg
landa = 0.002;  %%λ   %%The elimination rate of adenosine inhibitor vaccine
activeDc=0.003; %%Dact  %%DC activation rate
activeCTL=0.24; %%Cact  %%T lymphocyte activation rate
miu=0.45; %%µ   %%Effector cell cytotoxicity rate
effectadenoonmiu=0.001;  %%rAe  %%Effect of adenosine on effector cell cytotoxicity rate
effecthypoxionmiu=0.014;  %%rHe  %%Effect of hypoxia on effector cell cytotoxicity rate
inactiveDc=0.01; %%Dinact  %%DCs inactivation rate by Treg
landahypoxi=0.014; %%η   %%Elimination rate of hypoxia inhibitory vaccine
effectpxonoxgyen=0.00004; %%K   %%The effect of the vaccine on oxygen
pTreg=0.3;  %%PTreg    %%Treg cell production rate
pTumor=0.001;  %%Ptumor  %%Tumor growth rate
fh=2.3*10^(-3); %%¯f   %%Rate of oxygen consumption by the tumor cell
parametrs=[gamaa,effecthypoxionadenosine,r,y1,effectofadenosineonTumor,effectofhypoxiaonTumor,effectadenoonTreg,effecthypoxiaonTreg,landa,activeDc,activeCTL,miu,effectadenoonmiu,effecthypoxionmiu,inactiveDc,landahypoxi,effectpxonoxgyen,pTreg,pTumor,fh];
N=210;
x = [0, 1000, 2135, 2045, 792, 635, 0]; % ex. data
t1 = [0, 4, 16, 36, 48, 72 96]*30+(6*24)*30; % time in hours
dT = 1; % time step
t = 0:dT:18*24*30; % time vector
n = length(t);
Adeno = zeros(1,n);
vv = zeros(1,n);
days=2;
IT = [ 6,9,12,15]*30*24; % injection times for siRNA
vv(round(IT/dT)+1) =4200;
IT2= [ 6*24*30,8*24*30,10*24*30,13*24*30,15*24*30,17*24*30 ]; % injection times for px-478
v = zeros(1,n);
v(round(IT2/dT)+1) =40;
HD=zeros(1,n);
M=800;
landaw=1/6;
z=zeros(1,M);
gama=0.1;
x2(6*24*30)=0.001;
x3(6*24*30)=0.001;
y(1)=0;
Tumor=zeros(N+2,N+2);
xi=1:N+2;
yi=1:N+2;
A=.41;
B=.11;
center=(N+2)/2;
Tumor(center,center)=1;
gamaa=parametrs(1);
effecthypoxionadenosine=parametrs(2);
r=parametrs(3);
y1=parametrs(4);
effectofadenosineonTumor=parametrs(5);
effectofhypoxiaonTumor=parametrs(6);
effectadenoonTreg=parametrs(7);
effecthypoxiaonTreg=parametrs(8);
landa =parametrs(9);
activeDc=parametrs(10);
activeCTL=parametrs(11);
miu=parametrs(12);
miuu=miu;
effectadenoonmiu=parametrs(13);
effecthypoxionmiu=parametrs(14);
inactiveDc=parametrs(15);
landahypoxi=parametrs(16);
effectpxonoxgyen=parametrs(17);
pTreg=parametrs(18);
pTumor=parametrs(19);
fh=parametrs(20);
mm=randperm((N-1)*(N-1));
for i=1:N+1
    mm(mm==i)=[];
end
for i=1:N+1
    mm(mm==i*(N+2)+1)=[];
end
for i=1:N+2
    mm(mm==i*(N+2)-1)=[];
end
for j=1:100
    Tumor(mm(j))=4;%%100 inactive DC cells
end
for j=101:850
    Tumor(mm(j))=5; %%750 inactive CTL cells
end
Tumor(center,center)=1;
mov1=0;mov2=0;
Gam=10;
ns=1000;
for i=1:N+2
    for j=1:N+2
        savematrix4{i,j}=[ns,mov1,mov2,Gam];
        savematrix5{i,j}=[ns,mov1,mov2,Gam];
    end
end
%%Initial conditions for hypoxia
wo=0.8*ones(N+2,N+2);
xi=1:N+2;
yi=1:N+2;
cc=0;
l=2;
ao=2;
bo=6;
c=8;
d=2;
wo(1,:)=ao;
wo(2,:)=ao;
wo(N+1,:)=ao;
wo(N+2,:)=ao;
wo(:,1)=ao;
wo(:,2)=ao;
wo(:,N+1)=ao;
wo(:,N+2)=ao;
alpha=1;
co=2;
for i=1:n-1
    HD(i+1) = HD(i) + dT*( v(i)/dT - landahypoxi*HD(i) );
end
tic
for ii=6*24*30:18*24*30
    eff=0;
     %% DC Vaccine
     if (ii==7*24*30)
        vdc=find(Tumor==0); %%injection DC vaccine
  for i=1:N+1
    vdc(vdc==i)=[];
end
for i=1:N+1
    vdc(vdc==i*(N+2)+1)=[];
end
for i=1:N+2
    vdc(vdc==i*(N+2)-1)=[];
end
vdc=randsample(vdc,(70));
for j=1:length(vdc)
    Tumor(vdc(j))=6;
end
    end
    %% 
    % hypoxia
      co=co+1; % center of blocks = (co,co)
            if co==N+1
                co=2;
            end
            L1=co:-3:2;
            R1=co:3:N+1;R1(1)=[];
            ind=[L1(end:-1:1) R1]; %center of block = (ind(j), ind(k))
            %     Update oxygen
            for j=1:length(ind)
                for k=1:length(ind)
                    s1=ind(j);s2=ind(k);
                    if Tumor(s1,s2)==1
                        fh=parametrs(20);
                    else
                        fh=0;
                    end
                    wo(s1,s2)=alpha*1/9*(wo(s1+1,s2+1)+wo(s1+1,s2)+wo(s1+1,s2-1)+wo(s1,s2+1)+wo(s1,s2)+wo(s1,s2-1)+wo(s1-1,s2+1)+wo(s1-1,s2)+wo(s1-1,s2-1))+(1-alpha)*wo(s1,s2)-fh+(effectpxonoxgyen*HD(ii));
                    wo(s1-1,s2-1)=wo(s1,s2);wo(s1-1,s2)=wo(s1,s2);wo(s1-1,s2+1)=wo(s1,s2);wo(s1+1,s2)=wo(s1,s2);wo(s1+1,s2-1)=wo(s1,s2);wo(s1+1,s2+1)=wo(s1,s2);wo(s1,s2-1)=wo(s1,s2);wo(s1,s2+1)=wo(s1,s2);
                    wo(1,:)=ao;wo(2,:)=ao;wo(N+1,:)=ao;wo(N+2,:)=ao;wo(:,1)=ao;wo(:,2)=ao;wo(:,N+1)=ao;wo(:,N+2)=ao;
                end
            end
             hypoxi=1/(mean(mean(wo)));
             nhypoxi(ii)=1/(mean(mean(wo)));
            %%
            % adenosine
            t3=0:t(ii+1);
        Adeno(ii+1) = Adeno(ii) + dT*( vv(ii)/dT - landa*Adeno(ii) );
        Adeno=[z,Adeno];
        nn=M+ii;
        y(ii+1)=0;
        for i=0:M-1
            y(ii+1)=(Adeno(nn-i))+y(ii+1);
        end
        y(ii+1)=(y(ii+1)/M);
        Adeno(1:M)=[];
        x3(ii+1)=x3(ii)+0.1*(-r*x3(ii)+y1)+(effecthypoxionadenosine*hypoxi );
       x2(ii+1)=x2(ii)+0.1*((-r*x2(ii)+y1)+(effecthypoxionadenosine*hypoxi )-(gamaa*y(ii)*x2(ii)));
        %%
        % Tumor growth
        TN=Tumor==1;
        T11=(diff(TN,1,1)>0);T11=[0*Tumor(1,:);T11];
        T12=(diff(1-TN,1,1)>0);T12=[T12;0*Tumor(1,:)];
        T21=(diff(TN,1,2)>0);T21=[0*Tumor(:,1) T21];
        T22=(diff(1-TN,1,2)>0);T22=[T22 0*Tumor(:,1)];
        TN=(T11+T12+T21+T22)>0;
    [T1,T2]=find(TN==1);
    lT=randperm(length(T1));
    for j=1:length(lT)
        e=T1(lT(j));s=T2(lT(j));
        k=[Tumor(e-1,s-1),Tumor(e-1,s),Tumor(e-1,s+1),Tumor(e,s-1),Tumor(e,s+1),Tumor(e+1,s-1),Tumor(e+1,s+1),Tumor(e+1,s)];
        l=find(k==0);
          if (~isempty(l)) 
                    po=rand;
                    if po<(pTumor+effectofadenosineonTumor*x2(ii)+effectofhypoxiaonTumor*(1/wo(e,s)))
                        o=length(l);
                        u=randi(o);
                        k(l(u))=1;
                        Tumor(e-1,s-1)=k(1);
                        Tumor(e-1,s)=k(2);
                        Tumor(e-1,s+1)=k(3);
                        Tumor(e,s-1)=k(4);
                        Tumor(e,s+1)=k(5);
                        Tumor(e+1,s-1)=k(6);
                        Tumor(e+1,s+1)=k(7);
                        Tumor(e+1,s)=k(8);
                    end
          end
    end
    
            %%
            %immune cells
            %%%inactive DC
    [T1,T2]=find(Tumor==4);
     lo=randperm(length(T1));
    for j=1:length(lo)
        e=T1(lo(j));s=T2(lo(j));
          variable=savematrix4{e,s};
        if variable(1)<variable(4)
              if (2<e+variable(2))&&(e+variable(2)<N )&& (2<s+variable(3))&&(s+variable(3)<N )
                        if Tumor(e+variable(2),s+variable(3))==0
                            p1=[e-1,e,e+1];
                            p2=[s-1,s,s+1];
                            k1=[];
                            k=[];
                            h=1;
                            for i=1:3
                                for j=1:3
                                    if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                        k1{h}=[p1(i),p2(j)];
                                        k(h)=Tumor(p1(i),p2(j));
                                        h=h+1;
                                    end
                                end
                            end
                            h=find(k==1);
                            if isempty(h)
                                Tumor(e+variable(2),s+variable(3))=4;
                                Tumor(e,s)=0;
                                savematrix4{e+variable(2),s+variable(3)}=[variable(1)+1,variable(2),variable(3),variable(4)];
                            else
                                a=0;
                                b=1;
                                tr=a+(b-a)*rand;
                                if tr<(activeDc)
                                    Tumor(e+variable(2),s+variable(3))=6;
                                    Tumor(e,s)=0;
                                else
                                    Tumor(e+variable(2),s+variable(3))=4;
                                    Tumor(e,s)=0;
                                    savematrix4{e+variable(2),s+variable(3)}=[variable(1)+1,variable(2),variable(3),variable(4)];
                                end
                            end
                        else
                            p1=[e-1,e,e+1];
                            p2=[s-1,s,s+1];
                            k1=[];
                            k=[];
                            h=1;
                            for i=1:3
                                for j=1:3
                                    if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                        k1{h}=[p1(i),p2(j)];
                                        k(h)=Tumor(p1(i),p2(j));
                                        h=h+1;
                                    end
                                end
                            end
                            l=find(k==0);
                            h=find(k==1);
                            if (~isempty(l)) &&( isempty(h))
                                o=length(l);
                                u=randi(o);
                                g=l(u);
                                mb=k1{g};
                                if mb(1)>e
                                    mov1=+1;
                                elseif mb(1)==e
                                    mov1=0;
                                else
                                    mov1=-1;
                                end
                                if mb(2)>s
                                    mov2=+1;
                                elseif mb(2)==s
                                    mov2=0;
                                else
                                    mov2=-1;
                                end
                                Tumor(e+mov1,s+mov2)=4;
                                Tumor(e,s)=0;
                                v=e+mov1;
                                w=s+mov2;
                                savematrix4{v,w}=[variable(1)+1,mov1,mov2,variable(4)];
                            elseif  (~isempty(l)) && (~isempty(h))
                                a=0;
                                b=1;
                                tr=a+(b-a)*rand;
                                    o=length(l);
                                    u=randi(o);
                                    g=l(u);
                                    mb=k1{g};
                                    if mb(1)>e
                                        mov1=+1;
                                    elseif mb(1)==e
                                        mov1=0;
                                    else
                                        mov1=-1;
                                    end
                                    if mb(2)>s
                                        mov2=+1;
                                    elseif mb(2)==s
                                        mov2=0;
                                    else
                                        mov2=-1;
                                    end
                                       if tr<(activeDc)
                                    Tumor(e+mov1,s+mov2)=6;
                                    Tumor(e,s)=0;
                                       else
                                    Tumor(e+mov1,s+mov2)=4;
                                    Tumor(e,s)=0;
                                    v=e+mov1;
                                    w=s+mov2;
                                    savematrix4{v,w}=[variable(1)+1,mov1,mov2,variable(4)];
                                end
                            elseif  (isempty(l)) && (~isempty(h))
                                a=0;
                                b=1;
                                tr=a+(b-a)*rand;
                                if  tr<(activeDc)
                                    Tumor(e,s)=6;
                                end
                            end
                        end
                    else
                        p1=[e-1,e,e+1];
                        p2=[s-1,s,s+1];
                        k1=[];
                        k=[];
                        h=1;
                        for i=1:3
                            for j=1:3
                                if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                    k1{h}=[p1(i),p2(j)];
                                    k(h)=Tumor(p1(i),p2(j));
                                    h=h+1;
                                end
                            end
                        end
                        l=find(k==0);
                        h=find(k==1);
                        if (~isempty(l)) &&( isempty(h))
                            o=length(l);
                            u=randi(o);
                            g=l(u);
                            mb=k1{g};
                            if mb(1)>e
                                mov1=+1;
                            elseif mb(1)==e
                                mov1=0;
                            else
                                mov1=-1;
                            end
                            if mb(2)>s
                                mov2=+1;
                            elseif mb(2)==s
                                mov2=0;
                            else
                                mov2=-1;
                            end
                            Tumor(e+mov1,s+mov2)=4;
                            Tumor(e,s)=0;
                            variable(1)=variable(1)+1;
                            v=e+mov1;
                            w=s+mov2;
                            savematrix4{v,w}=[variable(1),mov1,mov2,variable(4)];
                        elseif  (~isempty(l)) && (~isempty(h))
                            a=0;
                            b=1;
                            tr=a+(b-a)*rand;
                                o=length(l);
                                u=randi(o);
                                g=l(u);
                                mb=k1{g};
                                if mb(1)>e
                                    mov1=+1;
                                elseif mb(1)==e
                                    mb(1)=0;
                                else
                                    mov1=-1;
                                end
                                if mb(2)>s
                                    mov2=+1;
                                elseif mb(2)==s
                                    mov2=0;
                                else
                                    mov2=-1;
                                end
                                 if tr<(activeDc)
                                Tumor(e+mov1,s+mov2)=6;
                                Tumor(e,s)=0;
                                 else
                                Tumor(e+mov1,s+mov2)=4;
                                Tumor(e,s)=0;
                                variable(1)=variable(1)+1;
                                v=e+mov1;
                                w=s+mov2;
                                savematrix4{v,w}=[variable(1),mov1,mov2,variable(4)];
                            end
                        elseif  (isempty(l)) && (~isempty(h))
                            a=0;
                            b=1;
                            tr=a+(b-a)*rand;
                            if  tr<(activeDc)
                                Tumor(e,s)=6;
                            end
                        end
                    end
                else
                    F=rand;
                    a=1/1.15;
                    T=1/((1-F)^a);
                    Gam=round(T);
                    p1=[e-1,e,e+1];
                    p2=[s-1,s,s+1];
                    k1=[];
                    k=[];
                    h=1;
                    for i=1:3
                        for j=1:3
                            if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                k1{h}=[p1(i),p2(j)];
                                k(h)=Tumor(p1(i),p2(j));
                                h=h+1;
                            end
                        end
                    end
                    l=find(k==0);
                    h=find(k==1);
                    if (~isempty(l)) &&( isempty(h))
                        o=length(l);
                        u=randi(o);
                        g=l(u);
                        mb=k1{g};
                        if mb(1)>e
                            mov1=+1;
                        elseif mb(1)==e
                            mov1=0;
                        else
                            mov1=-1;
                        end
                        if mb(2)>s
                            mov2=+1;
                        elseif mb(2)==s
                            mov2=0;
                        else
                            mov2=-1;
                        end
                        Tumor(e+mov1,s+mov2)=4;
                        Tumor(e,s)=0;
                        v=e+mov1;
                        w=s+mov2;
                        savematrix4{v,w}=[1,mov1,mov2,Gam];
                    elseif  (~isempty(l)) && (~isempty(h))
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                            o=length(l);
                            u=randi(o);
                            g=l(u);
                            mb=k1{g};
                            if mb(1)>e
                                mov1=+1;
                            elseif mb(1)==e
                                mov1=0;
                            else
                                mov1=-1;
                            end
                            if mb(2)>s
                                mov2=+1;
                            elseif mb(2)==s
                                mov2=0;
                            else
                                mov2=-1;
                            end
                             if tr<(activeDc)
                            Tumor(e+mov1,s+mov2)=6;
                            Tumor(e,s)=0;
                             else
                            Tumor(e+mov1,s+mov2)=4;
                            Tumor(e,s)=0;
                            v=e+mov1;
                            w=s+mov2;
                            savematrix4{v,w}=[1,mov1,mov2,Gam];
                        end
                    elseif  (isempty(l)) && (~isempty(h))
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                        if  tr<(activeDc)
                            Tumor(e,s)=6;
                        end
                    end
                end
    end
             %%
             %%inactive T cell
              [T1,T2]=find(Tumor==5);
     l5=randperm(length(T1));
    for j=1:length(l5)
        e=T1(l5(j));s=T2(l5(j));
           variable=savematrix5{e,s};
                if variable(1)<variable(4)
                    if (2<e+variable(2))&&(e+variable(2)<N )&& (2<s+variable(3))&&(s+variable(3)<N )
                        if Tumor(e+variable(2),s+variable(3))==0
                            p1=[e-1,e,e+1];
                            p2=[s-1,s,s+1];
                            k1=[];
                            k=[];
                            h=1;
                            for i=1:3
                                for j=1:3
                                    if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                        k1{h}=[p1(i),p2(j)];
                                        k(h)=Tumor(p1(i),p2(j));
                                        h=h+1;
                                    end
                                end
                            end
                            h=find(k==6);
                            if isempty(h)
                                Tumor(e+variable(2),s+variable(3))=5;
                                Tumor(e,s)=0;
                                savematrix5{e+variable(2),s+variable(3)}=[variable(1)+1,variable(2),variable(3),variable(4)];
                            else
                                a=0;
                                b=1;
                                tr=a+(b-a)*rand;
                                if tr<(activeCTL)
                                    a=0;
                                    b=1;
                                    tr1=a+(b-a)*rand;
                                    if tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                        Tumor(e+variable(2),s+variable(3))=8;
                                        Tumor(e,s)=0;
                                        savematrix7{e+variable(2),s+variable(3)}=[1];
                                    else
                                        Tumor(e+variable(2),s+variable(3))=7;
                                        Tumor(e,s)=0;
                                        savematrix8{e+variable(2),s+variable(3)}=[1];
                                    end
                                else
                                    Tumor(e+variable(2),s+variable(3))=5;
                                    Tumor(e,s)=0;
                                    savematrix5{e+variable(2),s+variable(3)}=[variable(1)+1,variable(2),variable(3),variable(4)];
                                end
                            end
                        else
                            p1=[e-1,e,e+1];
                            p2=[s-1,s,s+1];
                            k1=[];
                            k=[];
                            h=1;
                            for i=1:3
                                for j=1:3
                                    if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                        k1{h}=[p1(i),p2(j)];
                                        k(h)=Tumor(p1(i),p2(j));
                                        h=h+1;
                                    end
                                end
                            end
                            l=find(k==0);
                            h=find(k==6);
                            if (~isempty(l)) &&( isempty(h))
                                o=length(l);
                                u=randi(o);
                                g=l(u);
                                mb=k1{g};
                                if mb(1)>e
                                    mov1=+1;
                                elseif mb(1)==e
                                    mov1=0;
                                else
                                    mov1=-1;
                                end
                                if mb(2)>s
                                    mov2=+1;
                                elseif mb(2)==s
                                    mov2=0;
                                else
                                    mov2=-1;
                                end
                                Tumor(e+mov1,s+mov2)=5;
                                Tumor(e,s)=0;
                                v=e+mov1;
                                w=s+mov2;
                                savematrix5{v,w}=[variable(1),mov1,mov2,variable(4)];
                            elseif  (~isempty(l)) && (~isempty(h))
                                a=0;
                                b=1;
                                tr=a+(b-a)*rand;
                                    a=0;
                                    b=1;
                                    tr1=a+(b-a)*rand;
                                        o=length(l);
                                        u=randi(o);
                                        g=l(u);
                                        mb=k1{g};
                                        if mb(1)>e
                                            mov1=+1;
                                        elseif mb(1)==e
                                            mov1=0;
                                        else
                                            mov1=-1;
                                        end
                                        if mb(2)>s
                                            mov2=+1;
                                        elseif mb(2)==s
                                            mov2=0;
                                        else
                                            mov2=-1;
                                        end
                                            if tr<(activeCTL)
                                                 if tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                        Tumor(e+mov1,s+mov2)=8;
                                        Tumor(e,s)=0;
                                    else
                                        Tumor(e+mov1,s+mov2)=7;
                                        Tumor(e,s)=0;
                                    end
                                            else
                                    Tumor(e+mov1,s+mov2)=5;
                                    Tumor(e,s)=0;
                                    v=e+mov1;
                                    w=s+mov2;
                                    savematrix5{v,w}=[variable(1),mov1,mov2,variable(4)];
                                end
                            elseif  (isempty(l)) && (~isempty(h))
                                a=0;
                                b=1;
                                tr=a+(b-a)*rand;
                                if  tr<(activeCTL)
                                    a=0;
                                    b=1;
                                    tr1=a+(b-a)*rand;
                                    if tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                        Tumor(e,s)=8;
                                    else
                                        Tumor(e,s)=7;
                                    end
                                end
                            end
                        end
                    else
                        p1=[e-1,e,e+1];
                        p2=[s-1,s,s+1];
                        k1=[];
                        k=[];
                        h=1;
                        for i=1:3
                            for j=1:3
                                if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                    k1{h}=[p1(i),p2(j)];
                                    k(h)=Tumor(p1(i),p2(j));
                                    h=h+1;
                                end
                            end
                        end
                        l=find(k==0);
                        h=find(k==6);
                        if (~isempty(l)) &&( isempty(h))
                            o=length(l);
                            u=randi(o);
                            g=l(u);
                            mb=k1{g};
                            if mb(1)>e
                                mov1=+1;
                            elseif mb(1)==e
                                mov1=0;
                            else
                                mov1=-1;
                            end
                            if mb(2)>s
                                mov2=+1;
                            elseif mb(2)==s
                                mov2=0;
                            else
                                mov2=-1;
                            end
                            Tumor(e+mov1,s+mov2)=5;
                            Tumor(e,s)=0;
                            variable(1)=variable(1)+1;
                            v=e+mov1;
                            w=s+mov2;
                            savematrix5{v,w}=[variable(1),mov1,mov2,variable(4)];
                        elseif  (~isempty(l)) && (~isempty(h))
                            a=0;
                            b=1;
                            tr=a+(b-a)*rand;
                                a=0;
                                b=1;
                                tr1=a+(b-a)*rand;
                                    o=length(l);
                                    u=randi(o);
                                    g=l(u);
                                    mb=k1{g};
                                    if mb(1)>e
                                        mov1=+1;
                                    elseif mb(1)==e
                                        mov1=0;
                                    else
                                        mov1=-1;
                                    end
                                    if mb(2)>s
                                        mov2=+1;
                                    elseif mb(2)==s
                                        mov2=0;
                                    else
                                        mov2=-1;
                                    end
                                     if tr<(activeCTL)
                                           if tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                    Tumor(e+mov1,s+mov2)=8;
                                    Tumor(e,s)=0;
                                           else
                                    Tumor(e+mov1,s+mov2)=7;
                                    Tumor(e,s)=0;
                                end
                                     else
                                Tumor(e+mov1,s+mov2)=5;
                                Tumor(e,s)=0;
                                variable(1)=variable(1)+1;
                                v=e+mov1;
                                w=s+mov2;
                                savematrix5{v,w}=[variable(1),mov1,mov2,variable(4)];
                            end
                        elseif  (isempty(l)) && (~isempty(h))
                            a=0;
                            b=1;
                            tr=a+(b-a)*rand;
                            a=0;
                            b=1;
                            tr1=a+(b-a)*rand;
                            if tr<(activeCTL)
                            if  tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                Tumor(e,s)=8;
                            else
                                Tumor(e,s)=7;
                            end
                            end
                        end
                    end
                else
                    F=rand;
                    a=1/1.15;
                    T=1/((1-F)^a);
                    Gam=round(T);
                    p1=[e-1,e,e+1];
                    p2=[s-1,s,s+1];
                    k1=[];
                    k=[];
                    h=1;
                    for i=1:3
                        for j=1:3
                            if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                k1{h}=[p1(i),p2(j)];
                                k(h)=Tumor(p1(i),p2(j));
                                h=h+1;
                            end
                        end
                    end
                    l=find(k==0);
                    h=find(k==6);
                    if (~isempty(l)) &&( isempty(h))
                        o=length(l);
                        u=randi(o);
                        g=l(u);
                        mb=k1{g};
                        if mb(1)>e
                            mov1=+1;
                        elseif mb(1)==e
                            mov1=0;
                        else
                            mov1=-1;
                        end
                        if mb(2)>s
                            mov2=+1;
                        elseif mb(2)==s
                            mov2=0;
                        else
                            mov2=-1;
                        end
                        Tumor(e,s)=0;
                        Tumor(e+mov1,s+mov2)=5;
                        v=e+mov1;
                        w=s+mov2;
                        savematrix5{v,w}=[1,mov1,mov2,Gam];
                    elseif  (~isempty(l)) && (~isempty(h))
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                            a=0;
                            b=1;
                            tr1=a+(b-a)*rand;
                                o=length(l);
                                u=randi(o);
                                g=l(u);
                                mb=k1{g};
                                if mb(1)>e
                                    mov1=+1;
                                elseif mb(1)==e
                                    mov1=0;
                                else
                                    mov1=-1;
                                end
                                if mb(2)>s
                                    mov2=+1;
                                elseif mb(2)==s
                                    mov2=0;
                                else
                                    mov2=-1;
                                end
                                  if tr<(activeCTL)
                                      if tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                Tumor(e+mov1,s+mov2)=8;
                                Tumor(e,s)=0;
                                      else
                                Tumor(e+mov1,s+mov2)=7;
                                Tumor(e,s)=0;
                            end
                                  else
                            Tumor(e,s)=0;
                            Tumor(e+mov1,s+mov2)=5;
                            v=e+mov1;
                            w=s+mov2;
                            savematrix5{v,w}=[1,mov1,mov2,Gam];
                        end
                    elseif  (isempty(l)) && (~isempty(h))
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                        if  tr<(activeCTL)
                            a=0;
                            b=1;
                            tr1=a+(b-a)*rand;
                            if tr1<(pTreg+effectadenoonTreg*x2(ii)+effecthypoxiaonTreg*hypoxi)
                                Tumor(e,s)=8;
                            else
                                Tumor(e,s)=7;
                            end
                        end
                    end
                end
    end
            %%
            %%active DC
             [T1,T2]=find(Tumor==6);
     l6=randperm(length(T1));
    for j=1:length(l6)
        e=T1(l6(j));s=T2(l6(j));
        p1=[e-1,e,e+1];
                    p2=[s-1,s,s+1];
                    k1=[];
                    k=[];
                    h=1;
                    for i=1:3
                        for j=1:3
                            if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                k1{h}=[p1(i),p2(j)];
                                k(h)=Tumor(p1(i),p2(j));
                                h=h+1;
                            end
                        end
                    end
                    l2=find(k==0);
                    if ~isempty(l2)
                        o=length(l2);
                        u=randi(o);
                        k(l2(u))=6;
                        for i=1:length(k)
                            sh=k1{i};
                            Tumor(sh(1),sh(2))=k(i);
                        end
                        Tumor(e,s)=0;
                    end
    end
              %%
              %Effecror Cell
              [T1,T2]=find(Tumor==7);
     l7=randperm(length(T1));
    for j=1:length(l7)
        e=T1(l7(j));s=T2(l7(j));
          p1=[e-1,e,e+1];
                p2=[s-1,s,s+1];
                k1=[];
                k=[];
                h=1;
                for i=1:3
                    for j=1:3
                        if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                            k1{h}=[p1(i),p2(j)];
                            k(h)=Tumor(p1(i),p2(j));
                            h=h+1;
                        end
                    end
                end
                l1=find(k==1);
                  l2=find(k==0);
                if ~isempty(l1) && ~isempty(l2)
                    eff=eff+1;
                    a=0;
                    b=1;
                    tr=a+(b-a)*rand;
                    miuu=((miu)-(effectadenoonmiu*x2(ii))-effecthypoxionmiu*hypoxi)
                    if tr<miuu
                        o=length(l1);
                        u=randi(o);
                        k(l1(u))=0;
                        for i=1:length(k)
                            sh=k1{i};
                            Tumor(sh(1),sh(2))=k(i);
                        end
                    end
                            o=length(l2);
                            u=randi(o);
                            k(l2(u))=7;
                            for i=1:length(k)
                                sh=k1{i};
                                Tumor(sh(1),sh(2))=k(i);
                            end
                            Tumor(e,s)=0;
                elseif isempty(l1) && (~isempty(l2))
                        o=length(l2);
                        u=randi(o);
                        k(l2(u))=7;
                        for i=1:length(k)
                            sh=k1{i};
                            Tumor(sh(1),sh(2))=k(i);
                        end
                        Tumor(e,s)=0;
                elseif ~isempty(l1) && isempty(l2)
                    eff=eff+1;
                    a=0;
                    b=1;
                    tr=a+(b-a)*rand;
                    miuu=((miu)-(effectadenoonmiu*x2(ii))-effecthypoxionmiu*hypoxi)
                    if tr<miuu
                        o=length(l1);
                        u=randi(o);
                        k(l1(u))=0;
                        for i=1:length(k)
                            sh=k1{i};
                            Tumor(sh(1),sh(2))=k(i);
                        end
                    end
                end
    end
    %%
    % Treg
     [T1,T2]=find(Tumor==8);
     l8=randperm(length(T1));
    for j=1:length(l8)
        e=T1(l8(j));s=T2(l8(j));
        if e>2 && e<N && s>2 && s<N
                    k=[Tumor(e-1,s-1),Tumor(e-1,s),Tumor(e-1,s+1),Tumor(e,s-1),Tumor(e,s+1),Tumor(e+1,s-1),Tumor(e+1,s+1),Tumor(e+1,s)];
                    l2=find(k==6);
                    l1=find(k==0);
                    if ~isempty(l2) && isempty(l1)
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                        if tr<inactiveDc
                            o=length(l2);
                            u=randi(o);
                            k(l2(u))=4;
                            Tumor(e-1,s-1)=k(1);
                            Tumor(e-1,s)=k(2);
                            Tumor(e-1,s+1)=k(3);
                            Tumor(e,s-1)=k(4);
                            Tumor(e,s+1)=k(5);
                            Tumor(e+1,s-1)=k(6);
                            Tumor(e+1,s+1)=k(7);
                            Tumor(e+1,s)=k(8);
                        end
                    elseif ~isempty(l2) && ~isempty(l1)
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                                o=length(l1);
                                u=randi(o);
                                k(l1(u))=8;
                                Tumor(e-1,s-1)=k(1);
                                Tumor(e-1,s)=k(2);
                                Tumor(e-1,s+1)=k(3);
                                Tumor(e,s-1)=k(4);
                                Tumor(e,s+1)=k(5);
                                Tumor(e+1,s-1)=k(6);
                                Tumor(e+1,s+1)=k(7);
                                Tumor(e+1,s)=k(8);
                                Tumor(e,s)=0;
                                if tr<inactiveDc
                            o=length(l2);
                            u=randi(o);
                            k(l2(u))=4;
                                Tumor(e-1,s-1)=k(1);
                                Tumor(e-1,s)=k(2);
                                Tumor(e-1,s+1)=k(3);
                                Tumor(e,s-1)=k(4);
                                Tumor(e,s+1)=k(5);
                                Tumor(e+1,s-1)=k(6);
                                Tumor(e+1,s+1)=k(7);
                                Tumor(e+1,s)=k(8);
                        end
                    elseif isempty(l2) && ~isempty(l1)
                            o=length(l1);
                            u=randi(o);
                            k(l1(u))=8;
                            Tumor(e,s)=0;
                            Tumor(e-1,s-1)=k(1);
                            Tumor(e-1,s)=k(2);
                            Tumor(e-1,s+1)=k(3);
                            Tumor(e,s-1)=k(4);
                            Tumor(e,s+1)=k(5);
                            Tumor(e+1,s-1)=k(6);
                            Tumor(e+1,s+1)=k(7);
                            Tumor(e+1,s)=k(8);
                    end
                else
                    p1=[e-1,e,e+1];
                    p2=[s-1,s,s+1];
                    k1=[];
                    k=[];
                    h=1;
                    for i=1:3
                        for j=1:3
                            if (p1(i) >2) && (p1(i)<N)&& (p2(j) >2) &&(p2(j)<N)
                                k1{h}=[p1(i),p2(j)];
                                k(h)=Tumor(p1(i),p2(j));
                                h=h+1;
                            end
                        end
                    end
                    l1=find(k==0);
                    l2=find(k==6);
                    if ~isempty(l2) &&  ~isempty(l1)
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                                o=length(l1);
                                u=randi(o);
                                k(l1(u))=8;
                                for i=1:length(k)
                                    sh=k1{i};
                                    Tumor(sh(1),sh(2))=k(i);
                                end
                                Tumor(e,s)=0;
                                if tr<inactiveDc
                            o=length(l2);
                            u=randi(o);
                            k(l2(u))=4;
                            for i=1:length(k)
                                sh=k1{i};
                                Tumor(sh(1),sh(2))=k(i);
                            end
                                end
                    elseif ~isempty(l2) && isempty(l1)
                        a=0;
                        b=1;
                        tr=a+(b-a)*rand;
                        if tr<inactiveDc
                            o=length(l2);
                            u=randi(o);
                            k(l2(u))=4;
                            for i=1:length(k)
                                sh=k1{i};
                                Tumor(sh(1),sh(2))=k(i);
                            end
                        end
                    elseif isempty(l2) && ~isempty(l1)
                            o=length(l1);
                            u=randi(o);
                            k(l1(u))=8;
                            for i=1:length(k)
                                sh=k1{i};
                                Tumor(sh(1),sh(2))=k(i);
                            end
                            Tumor(e,s)=0;
                    end
                end
    end
    numberofTumor(ii)=length(find(Tumor==1));
     volumofTumor(ii)=(numberofTumor(ii)*0.01)*2.2;
     zaman=ii
         if mod(ii,200)==0
%         subplot(2,2,1)
        pcolor(xi,yi,Tumor)
        title('Tumor && immune system')
        colormap jet
        shading faceted
        drawnow
         end
end
ii=1:length(t)-1;
ii=ii/(24*30);
t1 = [0, 4, 16, 36, 48, 72 96]*30+(6*24)*30; % time in hours
t1=t1/(24*30);
figure,
plot(ii,volumofTumor)


