ncandidates = 15; %Population size
njobs = 10;
x = rand(ncandidates,njobs);
r1 = rand(1);
r2 = rand(1);
fx = zeros(ncandidates,1);

x_ = zeros(ncandidates,njobs);

for i = 1:ncandidates
    schedule = scheduler(x(i,:));
    [fx(i,1),CT,PM_Int] = MM_Cost();
end

S=xlsread('S_Original','sp'); %Input Parameters of Scheduling 
M=xlsread('M_integrated','mp'); %Input Parameters of Maintenance

gen_size = 4; 

for k=1:gen_size
    r1 = rand(1);
    r2 = rand(1);
    
    [fx,schedule,CT,PM_Int] = fx_builder(x);
    
    [f_best,x_best] = min(fx);
    [f_worst,x_worst] = max(fx);
    
    for i = 1:njobs
        for j = 1:ncandidates
            x_(j,i) = x(j,i)+r1*(x(j,x_best)-abs(x(j,i)))-r2*(x(j,x_worst)-abs(x(j,i)));
        end
    end
    
    [fx_,schedule_,CT_,PM_Int_] = fx_builder(x_);
    
    for t = 1:ncandidates
        if (fx_(t,1) < fx(t,1))
            x(t,:) = x_(t,:);
        end
    end
    
end