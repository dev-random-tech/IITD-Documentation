clc
clear all
S=xlsread('S_integrated','sp')% INPUT PARAMETERS OF SCHEDULING 
M=xlsread('M_integrated','mp')% INPUT PARAMETERS OF MAINTENANCE 
n=length(S(:,1));c=length(M(:,1));Rf=0.7;%OTHER INPUT PARAMETERS
PM=M(:,6)';icc=1.5;bs=500;ttr_cm=12;
pop_size=10;gen=10;pz=c*n+2*n+1;

% Initial Job Priority Vector (A1)
A1= 1+(n-1)*rand(pop_size,n);
ps=length(A1(:,1));
py=length(A1(1,:));
A2 = sort(A1,2,'descend');% arrange in non increasing order (A2)

for ri=1:ps 
    for rj=1:py
    A3(ri,rj)=find(A2(ri,rj)==A1(ri,:));
    end
end
A3;
% CALCULATION OF OBJECTIVE FUNCTION (INITIAL SOLUTION)
for ci=1:ps
PM_int=zeros(n,c);ai_1=zeros(n,c);ai=ai_1;Ai=ai_1;Aii=ai_1;PFC=zeros(n,c);
CT=zeros(n,1);PM_temp=PM;
for cj=1:py
        temp=A3(ci,cj);
        newP(cj,:)= S(temp,:);
end
newP;
CT(1,1)=newP(1,4); %INITIAL AGE OF COMPONENTS
ai(1,:)=ai_1(1,:)+newP(1,4)';Ai(1,:)=(ai(1,:)./M(:,2)').^(M(:,3)');
Aii(1,:)=(ai_1(1,:)./M(:,2)').^(M(:,3)');PFC(1,:)=exp(Aii(1,:)-Ai(1,:));    
for i=1:n %CALCULATION OF PM INTERVALS OF GROUPED COMPONENTS //START//
    if(i>1)
    CT(i,1)=CT(i-1,1)+newP(i,4);
    end
    for j=1:c
     if(PM_temp(1,j)<=CT(i,1)||0.9*PM_temp(1,j)<=CT(i,1))
         PM_int(i,j)=1;
         PM_temp(1,j)=CT(i,1)*PM_int(i,j)+PM(1,j);
     end
    end 
    MTR_PM=PM_int(i,:).*M(:,4)';
    CT(i,1)=CT(i,1)+max(MTR_PM);
     PM_temp;%//END//
end
PM_int;
CT;
%MAKING PM INTERVALS AND COMPLETION TIME IN ONE ROW
COM=horzcat(CT,PM_int);
count=c+1;temp1=1;temp2=count;
for i=1:n
    PMCT(1,temp1:temp2)=COM(i,:);
    temp1=temp1+count;temp2=temp2+count;  
end
PMCT;
for i=2:n %CALCULATION OF RELABILBITY OF EACH COMPONENT WHILE JOBS ARE PROCESSED (PFC)//START//
  for j=1:c
      ai_1(i,j)=ai(i-1,j)*(1-Rf*PM_int(i,j));
      ai(i,j)=ai_1(i,j)+newP(i,4);
      Ai(i,j)=(ai(i,j)/M(j,2))^M(j,3);
      Aii(i,j)= (ai_1(i,j)/M(j,2)')^(M(j,3)');
      PFC(i,j)=exp(Aii(i,j)-Ai(i,j));
  end
end      %//END//
PFC;   
PFM=(1-prod(PFC,2))';% PROBABILITY OF FAILUE OF MACHINE WHILE THE JOB IS PROCESSED (PFM) //START//END//
ft = PFM;JCT=zeros(n,n+1);%CALCULATION OF PROBABILITY MASS FUNCTION (PMF)//START//
f_t = 1-ft;
result = zeros(length(ft),length(ft)+1);
index = 1:length(ft);
for i = index
   for j = 1:i
       comb = nchoosek(index(1:i),length(index(1:i))+1-j);
       s = 0;
       for k = 1:size(comb,1)
           newindex = index(1:i);
           newindex([comb(k,:)]) = [];
           s = s+prod(ft(comb(k,:)))*prod(f_t(newindex));          
       end
       result(i,j) = s;
   end
   result(i,i+1) = prod(f_t(1:i));
end %//END//
PMF=result;
for i=1:n% CALCULATION OF TOTAL PENALY COST DUE TO BATCH  DELAY//START//
    for j=1:i
        JCT(i,j)=ttr_cm*(i-(j-1));
    end
end
JCT;
for i=1:n
    for j=1:i+1
        JCTF(i,j)=JCT(i,j)+CT(i,1);
    end
end
JCTF;
for i=1:n
    for j=1:n+1
        JT(i,j)=max(JCTF(i,j)-newP(i,2),0);
    end
end
JT;
for i=1:n
    for j=1:n+1
        JPC(i,j)=PMF(i,j)*JT(i,j)*newP(i,3);
    end
end
JPC;
TPC=sum(JPC(:));
% fprintf('Total penalty cost due to batch delay = %0.2f \n',TPC/CT(n,1))%//END//
for i=1:n%CALCULATION OF TOTAL INVENTORY CARRYING COST OF THE SYSTEM //START//
   if(CT(i,1)<=newP(i,2))
       ICC(i,1)=0;
       continue
   end
   if(CT(i,1)>newP(i,2))
     gr=find(CT(:,1)>newP(i,2),1,'first');
   end
   if(gr==i)
       ICC(i,1)=icc*(bs/2)*(CT(i,1)-newP(i,2));
   else
     ICC(i,1)=icc*bs*(CT(i-1,1)-newP(i,2))+ icc*(bs/2)*(CT(i,1)-CT(i-1,1));
   end
end
TICC=sum(ICC(:));
% fprintf('Total inventory carrying cost of the system= %0.2f \n',TICC/CT(n,1))%//END//
ECPUT = (TPC+TICC)/CT(n,1); %CALCULATION OF ECPUT OF THE SYSTEM//START
% fprintf('Total Expected Cost Per Unit Time of the system E[CPUT]s*m= %0.2f \n',ECPUT)%//END//
A3(ci,py+1:(pz-1))=PMCT;A3(ci,pz)=ECPUT;
end
A3;
A1(:,py+1:pz)=A3(:,py+1:pz);
A1;
[min_value, min_row]=(min(A1(:,pz)));
[max_value, max_row]=(max(A1(:,pz)));
C1=A1;
C3=A3;
for it=1:gen

for i=1:ps
    B1(i,1:py)= C1(i,1:py)+ rand(1,py).*(C1(min_row,1:py)- C1(i,1:py))...
    - rand(1,py).*(C1(max_row,1:py)-C1(i,1:py));
end
B2 = sort(B1,2,'descend'); %NON INCREASING SORTED DIFFERENCE MEAN MATRIX (B2)
for ri=1:ps %CONVERSION OF SORTED MATRIX(B2) TO JOB PERMUTATION MATRIX(B3)
    for rj=1:py
    B3(ri,rj)=find(B2(ri,rj)==B1(ri,:));
    end
end
B1;
B3;
%CALCULATION OF OBJECTIVE FUNCTION (SOLUTION COMPARISON)
for ci=1:ps
PM_int=zeros(n,c);ai_1=zeros(n,c);ai=ai_1;Ai=ai_1;Aii=ai_1;PFC=zeros(n,c);
CT=zeros(n,1);PM_temp=PM;
for cj=1:py
        temp=B3(ci,cj);
        newP(cj,:)= S(temp,:);
end
newP;
CT(1,1)=newP(1,4); %INITIAL AGE OF COMPONENTS
ai(1,:)=ai_1(1,:)+newP(1,4)';Ai(1,:)=(ai(1,:)./M(:,2)').^(M(:,3)');
Aii(1,:)=(ai_1(1,:)./M(:,2)').^(M(:,3)');PFC(1,:)=exp(Aii(1,:)-Ai(1,:));    
for i=1:n %CALCULATION OF PM INTERVALS OF GROUPED COMPONENTS //START//
    if(i>1)
    CT(i,1)=CT(i-1,1)+newP(i,4);
    end
    for j=1:c
     if(PM_temp(1,j)<=CT(i,1)||0.9*PM_temp(1,j)<=CT(i,1))
         PM_int(i,j)=1;
         PM_temp(1,j)=CT(i,1)*PM_int(i,j)+PM(1,j);
     end
    end 
    MTR_PM=PM_int(i,:).*M(:,4)';
    CT(i,1)=CT(i,1)+max(MTR_PM);
     PM_temp;%//END//
end
PM_int;
CT;
%MAKING PM INTERVALS AND COMPLETION TIME IN ONE ROW
COM=horzcat(CT,PM_int);count=c+1;temp1=1;temp2=count;
for i=1:n
    PMCT(1,temp1:temp2)=COM(i,:);
    temp1=temp1+count;temp2=temp2+count;  
end
PMCT;
for i=2:n %CALCULATION OF RELABILBITY OF EACH COMPONENT WHILE JOBS ARE PROCESSED (PFC)//START//
  for j=1:c
      ai_1(i,j)=ai(i-1,j)*(1-Rf*PM_int(i,j));
      ai(i,j)=ai_1(i,j)+newP(i,4);
      Ai(i,j)=(ai(i,j)/M(j,2))^M(j,3);
      Aii(i,j)= (ai_1(i,j)/M(j,2)')^(M(j,3)');
      PFC(i,j)=exp(Aii(i,j)-Ai(i,j));
  end
end      %//END//
PFC;   
PFM=(1-prod(PFC,2))';% PROBABILITY OF FAILUE OF MACHINE WHILE THE JOB IS PROCESSED (PFM) //START//END//
ft = PFM;JCT=zeros(n,n+1);%CALCULATION OF PROBABILITY MASS FUNCTION (PMF)//START//
f_t = 1-ft;
result = zeros(length(ft),length(ft)+1);
index = 1:length(ft);
for i = index
   for j = 1:i
       comb = nchoosek(index(1:i),length(index(1:i))+1-j);
       s = 0;
       for k = 1:size(comb,1)
           newindex = index(1:i);
           newindex([comb(k,:)]) = [];
           s = s+prod(ft(comb(k,:)))*prod(f_t(newindex));          
       end
       result(i,j) = s;
   end
   result(i,i+1) = prod(f_t(1:i));
end %//END//
PMF=result;
for i=1:n% CALCULATION OF TOTAL PENALY COST DUE TO BATCH  DELAY//START//
    for j=1:i
        JCT(i,j)=ttr_cm*(i-(j-1));
    end
end
JCT;
for i=1:n
    for j=1:i+1
        JCTF(i,j)=JCT(i,j)+CT(i,1);
    end
end
JCTF;
for i=1:n
    for j=1:n+1
        JT(i,j)=max(JCTF(i,j)-newP(i,2),0);
    end
end
JT;
for i=1:n
    for j=1:n+1
        JPC(i,j)=PMF(i,j)*JT(i,j)*newP(i,3);
    end
end
JPC;
TPC=sum(JPC(:));
% fprintf('Total penalty cost due to batch delay = %0.2f \n',TPC/CT(n,1))%//END//
for i=1:n%CALCULATION OF TOTAL INVENTORY CARRYING COST OF THE SYSTEM //START//
   if(CT(i,1)<=newP(i,2))
       ICC(i,1)=0;
       continue
   end
   if(CT(i,1)>newP(i,2))
     gr=find(CT(:,1)>newP(i,2),1,'first');
   end
   if(gr==i)
       ICC(i,1)=icc*(bs/2)*(CT(i,1)-newP(i,2));
   else
     ICC(i,1)=icc*bs*(CT(i-1,1)-newP(i,2))+ icc*(bs/2)*(CT(i,1)-CT(i-1,1));
   end
end
TICC=sum(ICC(:));
% fprintf('Total inventory carrying cost of the system= %0.2f \n',TICC/CT(n,1))%//END//
ECPUT = (TPC+TICC)/CT(n,1); %CALCULATION OF ECPUT OF THE SYSTEM//START
% fprintf('Total Expected Cost Per Unit Time of the system E[CPUT]s*m= %0.2f \n',ECPUT)%//END//
B3(ci,py+1:(pz-1))=PMCT;B3(ci,pz)=ECPUT;
end
B1(:,py+1:pz)=B3(:,py+1:pz);
B3;
A3;
for i =1:ps % COMPARISON OF MATRIX B3/B2 WITH A3/A2 TO GET (C1/C2)
    if(B3(i,pz)<=C3(i,pz))
        C1(i,:)= B1(i,:);
        C3(i,:)=B3(i,:);
    end
end
A3
B3
C3
[min_value, min_row]=(min(C1(:,py+1)));
[max_value, max_row]=(max(C1(:,py+1)));
B1(:,py+1:pz)=[];
B3(:,py+1:pz)=[];
end
[Min_cput, row]=(min(C3(:,pz)))
Optimum_Solution= C3(row,:)
