close all
clear all
%%%%%%%%parameters%%%%%%%%%%
delta=0.025;                 %depreciation rate
alpha=0.35;                  %capital share of production
beta=0.99;                   %discount rate
sigma=2;                     %utility function parameter
pi=[0.977 0.023;0.074 0.926];%transition matrix

%%%%%%%%k grid%%%%%%%%%%%%%%
ah = 1.1;                    %high productivity 
al=0.678;                    %low productivity
knum=1000;                   %set number of k equals 1000
kmin=0 ;                     %set k min equals 0
%?*syms kmax
%?eqn=kmax==ah*kmax^alpha+(1-delta)*kmax
%?kmax=solve(eqn,kmax)        %calculate the k max[tooooo big
%kmax = 1.1*(alpha*ah/(1/beta-1+delta))^(alpha/(1-alpha))...
    %+(1-delta)*(alpha*1.1/(1/beta-1+delta))^(1/(1-alpha));
kmax=45;
k=linspace(kmin,kmax,knum);  %decide on k grid
kmat=repmat(k',[1 knum]);    %this is a matrix which will be useful later

%%%%%%%%variables%%%%%%%%%%%
conh=ah*kmat.^alpha+(1-delta)*kmat-kmat';
conl=al*kmat.^alpha+(1-delta)*kmat-kmat';
reth=conh.^(1-sigma)/(1-sigma);
retl=conl.^(1-sigma)/(1-sigma);
reth(conh < 0)=-Inf;
retl(conl < 0)=-Inf;


%%%%%%%%Iteration%%%%%%%%%%%
dis = 1; tol = 1e-06;        % tolerance for stopping 
v_guess = zeros(2, knum);   %?????
while dis > tol
    
    vh_mat = reth + beta *(pi(1,1)* repmat(v_guess(1,:), [knum 1])...
        +pi(1,2)*repmat(v_guess(2,:), [knum 1]));
    vl_mat = retl + beta *pi(2,2)* repmat(v_guess(2,:), [knum 1])+...
        beta *pi(2,1)* repmat(v_guess(1,:), [knum 1]);
    
    [vfnh, ph_indxh] = max(vh_mat, [], 2);
     vfnh = vfnh';
    [vfnl, pl_indxl] = max(vl_mat, [], 2);
    vfnl = vfnl';
    
    
    dis =[max(abs(vfnl-v_guess(2,:)));max(abs(vfnh - v_guess(1,:)))]    %?????
    v_guess =[vfnh;vfnl];
    %?????
    
end

gh = k(ph_indxh);
gl = k(pl_indxl);

figure
plot(k,vfnh,'g');
xlabel('k');
ylabel('value function');
title('value functions for high state');

figure
plot (k,vfnl,'b');
xlabel('k');
ylabel('value function');
title('value functions for low state');


figure 
plot (k, gh,'g');
xlabel('k');
ylabel('policy function');
title('policy functions for high state and low state') ;
hold on;
plot(k , gl,'b');
hold off;