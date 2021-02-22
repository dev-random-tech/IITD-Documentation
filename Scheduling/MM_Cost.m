function [cost,completion_time,PM_Intervals] = MM_Cost()
    S=xlsread('S_integrated'); %Input Parameters of Scheduling 
    M=xlsread('M_integrated','mp'); %Input Parameters of Maintenance

    m=5; %Number of Machines
    n=length(S(:,1)); %Number of jobs
    c=length(M(:,1)); %Total Number of Components
    Rf=0.7; %Restoring Factor
    PM=M(:,8)'; %Optimal PM 
    PM_temp=PM;
    MTTR=M(:,4); %Time to repair the machine
    Due=S(:,2);
    eta=M(:,2);
    beta=M(:,3);
    PM_int=zeros(n,c); %The matrix to indicated whether PM has been done or not
    comp_ind=[0,3,6,9,12,15]; %To determine the indices of the starting components of machines
    PFC=zeros(n,c); %Probability of failure of components while job is being processed
    i_PFcheck=m+abs(n-m); %To determine the no. of parallel checks to be done,
    %while the various jobs pass through the machines
    PFM=zeros(n,m); %Probability of M/C Failure while job is processed
    PMF=zeros(n,n+1,m); %Probability Mass Function, each sheet is for a machine
    JCT=zeros(n,n+1,m); %JCT(i,j,k)=failure time when ith job is performed and 
    %j is the number of times the machine has failed for kth machine 

    %Other parameters
    ttr_cm=12;PR=20;Clp=40;LC=500;
    icc=1.5; %Inventory Carrying Cost
    bs=500;Cd=PR*Clp+LC;F_cm=19645;DCMC=ttr_cm*Cd;

    mttr_cm_1=round((12+8+18)/3);
    mttr_cm_2=round((11+10+25)/3);
    mttr_cm_3=round((22+20+18)/3);
    mttr_cm_4=round((12+15+11)/3);
    mttr_cm_5=round((21+22+20)/3);

    PT=S(:,3+1:3+m); %Processing times
    mrt=PT; %Machine Running Time after ith job is completed in jth machine [i,j]
    for i=2:size(mrt,1)
        mrt(i,:)=sum(PT(1:i,:),1);
    end

    mul=zeros(n,n+1);
    for i=1:n
        for j=1:n
            if(j>i)
                break
            else
                mul(i,j)=i+1-j;
            end
        end
    end

    for i=1:m
        eval(['JCT(:,:,i)=' 'mttr_cm_' num2str(i) '*mul']);
    end

    JT=zeros(n,n+1,m); %Job Tardiness
    JCTF=zeros(n,n+1,m); %Job Completion Time With Failure
    PC=S(:,3); %(Penalty Cost)/hr.
    JPC=zeros(n,n+1,m);

    ai_1=zeros(n,c); %Virtual Age
    ai=zeros(n,c); %Actual Age  
    Ai=zeros(n,c); %Part Probability of Actual Age
    Aii=zeros(n,c); %Part Probability of Virtual Age
    CT=zeros(n,m); %Completion Time of various jobs in various machines
    CT(1,1)=S(1,4); %As first job in first machine will have no other PM


    PM_int=zeros(n,c);
    i_n=1;i_m=1;
    %Direct calculations without PM for the first job in first machine
    ai(i_n,(comp_ind(i_m)+1):comp_ind(i_m+1))=ai_1(i_n,(comp_ind(i_m)+1):comp_ind(i_m+1))+S(1,3+i_m);
    Ai(i_n,(comp_ind(i_m)+1):comp_ind(i_m+1))=(ai(i_n,comp_ind(i_m)+1:comp_ind(i_m+1))./M(comp_ind(i_m)+1:comp_ind(i_m+1),2)').^(M((comp_ind(i_m)+1):comp_ind(i_m+1),3)');
    Aii(i_n,(comp_ind(i_m)+1):comp_ind(i_m+1))=(ai_1(i_n,comp_ind(i_m)+1:comp_ind(i_m+1))./M(comp_ind(i_m)+1:comp_ind(i_m+1),2)').^(M((comp_ind(i_m)+1):comp_ind(i_m+1),3)');

    PFC(i_n,comp_ind(i_m)+1:comp_ind(i_m+1))=exp(Aii(1,comp_ind(i_m)+1:comp_ind(i_m+1))-Ai(1,comp_ind(i_m)+1:comp_ind(i_m+1)));

    i_n=1; %index for jobs
    %For first machine, iterating over various jobs
    i_m=1; %index for machines

    PMTimeScale=zeros(n,c);

    for i=1:n 
        %'CT' [n,m]. Thus this iterates over rows in CT
        if i>1
            CT(i,i_m)=CT(i-1,i_m)+S(i,3+i_m);%Machine completion time starts after 3rd column
        end
        for j=comp_ind(i_m)+1:comp_ind(i_m+1)
            if(PM_temp(1,j)<=mrt(i,i_m)||0.8*PM_temp(1,j)<=mrt(i,i_m))
                PM_int(i,j)=1;
                PM_temp(1,j)=mrt(i,i_m)+PM_int(i,j)*PM(1,j);
            end
            %Calculation of Reliability of each component while jobs are processed
            if i>1 
                ai_1(i,j)=ai(i-1,j)*(1-Rf*PM_int(i,j));
                ai(i,j)=ai_1(i,j)+S(i,3+i_m);
                Ai(i,j)=(ai(i,j)/M(j,2))^M(j,3);
                Aii(i,j)= (ai_1(i,j)/M(j,2)')^(M(j,3)');
                PFC(i,j)=exp(Aii(i,j)-Ai(i,j));
            end
        end 
        MTR_PM=PM_int(i,:).*MTTR';
        CT(i,i_m)=CT(i,i_m)+max(MTR_PM(comp_ind(i_m)+1:comp_ind(i_m+1)));   

    end

    %For the first row: One job over various machines
    i_n=1;
    for i=2:m
        %'CT' [n,m]. Thus this iterates over columns in CT
        if i>1
            CT(i_n,i)=CT(i_n,i-1)+S(i_n,3+i);
        end
        for j=(comp_ind(i)+1):comp_ind(i+1)
            if(PM_temp(1,j)<=mrt(i_n,i)||0.8*PM_temp(1,j)<=mrt(i_n,i))
                PM_int(1,j)=1;
                PM_temp(1,j)=mrt(i_n,i)+PM(i_n,j)*PM_int(i_n,j);
            end
            ai(i_n,(comp_ind(i)+1):comp_ind(i+1))=ai_1(i_n,(comp_ind(i)+1):comp_ind(i+1))+S(1,3+i);
            Ai(i_n,(comp_ind(i)+1):comp_ind(i+1))=(ai(i_n,comp_ind(i)+1:comp_ind(i+1))./M(comp_ind(i)+1:comp_ind(i+1),2)').^(M((comp_ind(i)+1):comp_ind(i+1),3)');
            Aii(i_n,(comp_ind(i)+1):comp_ind(i+1))=(ai_1(i_n,comp_ind(i)+1:comp_ind(i+1))./M(comp_ind(i)+1:comp_ind(i+1),2)').^(M((comp_ind(i)+1):comp_ind(i+1),3)');

            PFC(i_n,comp_ind(i)+1:comp_ind(i+1))=exp(Aii(1,comp_ind(i)+1:comp_ind(i+1))-Ai(1,comp_ind(i)+1:comp_ind(i+1)));
        end 

        MTR_PM=PM_int(i_n,:).*MTTR';
        CT(i_n,i)=CT(i_n,i)+max(MTR_PM(comp_ind(i)+1:comp_ind(i+1)));

    end

    %Generalised CT
    for i_n=2:n %for jobs
        for i_m=2:m %for machines

            CT(i_n,i_m)=max(CT(i_n-1,i_m),CT(i_n,i_m-1));
            CT(i_n,i_m)=CT(i_n,i_m)+S(i_n,3+i_m);


            for j=(comp_ind(i_m)+1):comp_ind(i_m+1)
                if(PM_temp(1,j)<=mrt(i_n,i_m)||0.8*PM_temp(1,j)<=mrt(i_n,i_m))
                    PM_int(i_n,j)=1;
                    PM_temp(1,j)=mrt(i_n,i_m)+PM(1,j)*PM_int(i_n,j);
                end
                ai_1(i_n,j)=ai(i_n-1,j)*(1-Rf*PM_int(i_n,j));
                ai(i_n,j)=ai_1(i_n,j)+S(i_n,3+i_m);
                Ai(i_n,j)=(ai(i_n,j)/M(j,2))^M(j,3);
                Aii(i_n,j)= (ai_1(i_n,j)/M(j,2)')^(M(j,3)');
                PFC(i_n,j)=exp(Aii(i_n,j)-Ai(i_n,j));
            end

            MTR_PM=PM_int(i_n,:).*MTTR';
            CT(i_n,i_m)=CT(i_n,i_m)+max(MTR_PM(comp_ind(i_m)+1:comp_ind(i_m+1)));        

        end
    end

    %To find out PMF
    for i_m=1:m
        PFCtemp=PFC(:,comp_ind(i_m)+1:comp_ind(i_m+1));
        PFM(:,i_m)=(1-prod(PFCtemp,2));
        PMF(:,:,i_m)=pmcfail(PFM);
    end

    %To find JCTF
    for k=1:m
        for i=1:n
            for j=1:i+1
                JCTF(i,j,k)=JCT(i,j)+CT(i,k);
            end
        end
    end

    %To find JT, JPC
    for k=1:m
        for i=1:n
            for j=1:n+1
                JT(i,j,k)=max(JCTF(i,j,k)-Due(i,1),0);
                JPC(i,j,k)=PMF(i,j,k)*JT(i,j,k)*PC(i,1);
            end
        end
    end

    for i=1:m
        for j=1:n
            PMTimeScale(j,comp_ind(i)+1:comp_ind(i+1))=PM_int(j,comp_ind(i)+1:comp_ind(i+1)).*CT(j,i);
        end
    end

    TPC=sum(JPC(:)); %Total Penalty Cost per hour
    BWC=0; %Batch Waiting Cost
    DueTemp=Due;
    for j=1:n
        for i=1:m
            if DueTemp(j,1)<CT(j,i)
                BWC=BWC+(CT(j,i)-DueTemp(j,1))*icc*bs/2;
                if i<m 
                    DueTemp(j,1)=CT(j,i);
                end
                if i<m && j>1
                    if CT(j,i)<CT(j-1,i+1)
                        DueTemp(j,1)=CT(j-1,i+1);
                        BWC=BWC+(CT(j-1,i+1)-CT(j,i))*icc*bs;

                    end
                end
            end

        end
    end
    TotalCost=(TPC+BWC)/CT(n,m);
    cost = TotalCost;
    completion_time = CT;
    PM_Intervals = PM_int;
    
end