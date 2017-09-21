close all;
clear all;

%%%%%%%%parameters%%%%%%%%%%
delta=0.025;                 %depreciation rate
alpha=0.35;                  %capital share of production
beta=0.99;                   %discount rate
sigma=2;                     %utility function parameter
pi=[0.977 0.023;0.074 0.926];%transition matrix
pi_s=pi^1000;                %invariant distribution
ah=1e-10+1;                  %set a value for high ptoductivity shock
al=(1-ah*pi_s(1))/pi_s(3);   %calculate low prod shock

%%%%%%%%simulation%%%%%%%%%%
nStates = 1000;
initialProbabilityState=[1 0];
states = zeros(nStates,2);
states(1,:) = initialProbabilityState;
for ns = 2:nStates
    states(ns,:) =states(ns-1,:)*pi ;
end                          %using markov chain to generate a prob seq

for ns=1:nStates
    r=rand();
    if states(ns,1)>r;
        a(ns)=ah;
    else a(ns)=al;
    end
end                         %get a sequence of a, length of 1000


%%%%%%%set k grid%%%%%%%%%%
knum=1000;                   %set number of k equals 1000
kmin=0 ;                    
kmax=45;
k=linspace(kmin,kmax,knum);  %decide on k grid
kmat=repmat(k',[1 knum]); 

%%%%%%%std of output%%%%%%
for ns=1:1000
y(ns)=a(ns)*k(ns).^alpha;
y_std=std(y);
end