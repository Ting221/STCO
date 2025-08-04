function S = STCO1_run(data,lambda,gamma,rho,eta,tao)

%% Variables Initialization
[T,m]=size(data);
b0_hat=zeros(m,1);      %initialize the portfolio before transaction
b1=ones(m,1)/m;         %initialize the portfolio in t=1
S=zeros(T,1);           %record the cumulative wealth in each iteration
s0=1;                   %initialize the wealth before transaction
w0=w(b0_hat,b1,gamma);  %net wealth proportion, which from the formula 1=w_t-1+gamma*||b_t-1_hat-b_t*w_t-1||_1(refer paper[1])
b_t_hat=zeros(m,T);     %Record the proportion of assets after the end of the day's trading

%% main
%compute the cumulative wealth in t=1        
x1=data(1,:)';

run_ret =s0*w0; 
s1=run_ret*(b1'*x1);
S(1)=s1;   

b1_hat=(b1.*x1)/(b1'*x1);
b_t_hat(:,1)=b1_hat;  

%compute the cumulative wealth in t=2:T 
for t=2:T
    b_tcut1_hat=b_t_hat(:,t-1);
    s_tcut1=S(t-1);
    b_t=DENRPO1_fun(data,t-1,b_tcut1_hat,rho,lambda,eta,tao); 
    b_t=b_t/norm(b_t,1);
    w_tcut1=w(b_tcut1_hat,b_t,gamma);
    run_ret =s_tcut1*w_tcut1;
    x_t=data(t,:)';
    s_t=run_ret*(b_t'*x_t);
    S(t)=s_t;
    b_t_hat(:,t)=(b_t.*x_t)/(b_t'*x_t);
end
% if S(end)<10000
%     fprintf('\t %.2f \n',S(end));       
% else
%     fprintf('\t %.2e \n',S(end));         
% end
S=S(end);