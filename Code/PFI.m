clear, clc
l=1;
sigma_epsi=0.2;
rho=0.5;
m=3;
nz=5;
[z,zprob] =TAUCHEN(nz,rho,sigma_epsi,m);
PI=zprob;
P=ones(1,nz)/nz;
dis=1;
while abs(dis)>0;
    P1=P*zprob;
    dis=P1-P;
    P=P1;
end
Z=exp(z);
L=P*Z*l;
a_lo = 0;
a_hi = 80;
num_a = 500;
a = linspace(a_lo, a_hi, num_a);
alpha=1/3;
beta=0.99;
delta=0.025;
sigma=2;
r_max=1/beta;
r_min=1;
k_max=((r_min-1+delta)/alpha)^(1/(alpha-1))*L;
k_min=((r_max-1+delta)/alpha)^(1/(alpha-1))*L;
k_guess=1;
aggsav=0;
while abs(k_guess-aggsav)>=0.01;
    k_guess=1/2*(k_min+k_max);
    %r_guess=1/2*(r_min+r_max);
    %k_guess=((r_guess-1+delta)/alpha)^(1/(alpha-1))*L;
r=alpha*(k_guess/L)^(alpha-1)+1-delta;
w=(1-alpha)*(k_guess/L)^(alpha);
y_s=Z'*w*l;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, r*a', a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(nz, num_a);
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.000001;
     for i=1:nz;
        E(:,:,i)=repmat(PI(i,:)*v_guess, [num_a 1]);
     end
        v=ret+beta*E;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vmax, pol_indx]=max(v, [], 2);
        vmax=permute(vmax, [3 1 2]);
    
        V=reshape(vmax, [nz*num_a 1]);
        pol_indx1=permute(pol_indx, [3 1 2]);
        pol_fn1 = a(pol_indx1);
        cons1=r*repmat(a,[nz 1])-pol_fn1;
        cons1 = cons1+repmat(y_s', [1 num_a]);
        ret1 = (cons1 .^ (1-sigma)) ./ (1 - sigma);
        ret1(cons1<0)=-Inf;
        ret1=reshape(ret1, [nz*num_a,1]);
        Q = makeQmatrix(pol_indx1, PI);
        for j=1:30
        vNew=ret1+beta*Q*V;
        v=vNew;
        j=j+1;
        end
        vmax=reshape(v,[nz num_a]);
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:nz;
        E(:,:,i)=repmat(PI(i,:)*vmax, [num_a 1]);
        end
        v=ret+beta*E;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vmax, pol_indx]=max(v, [], 2);
        vmax=permute(vmax, [3 1 2]);
        v_tol = max(max(abs(vmax-v_guess)));
        v_guess=vmax;
    end;
     % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn = a(pol_indx);
    
   % SET UP INITITAL DISTRIBUTION
   Mu = zeros(nz,num_a);
  Mu(1, 4) = 1;
    %Mu=ones(nz, num_a)/(nz*num_a);
    % ITERATE OVER DISTRIBUTIONS
    dif=1;
    while dif >10^-7;
          [emp_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
          MuNew = zeros(size(Mu));
          for ii = 1:length(emp_ind)
              apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
              MuNew(:, apr_ind) = MuNew(:, apr_ind) + (PI(emp_ind(ii),:)*Mu(emp_ind(ii),a_ind(ii)))';
          end
          dif=max(max(abs(MuNew-Mu)));
          Mu=MuNew;
    end
    aggsav=sum(Mu*a');
    if aggsav>=k_guess;
       k_min=k_guess;
    else
       k_max=k_guess;
    end
    [k_guess,aggsav,k_max,k_min,r]
    if abs(k_max-k_min)<10^(-5);
       break;
   end
end
r_cm=1/beta
figure(1)
plot(a,pol_fn)
legend('State 1','State 2','State 3','State 4','State 5','location','northwest')
title(['Policy Function'])
Mu1=Mu';
pop=[Mu1(:,1);Mu1(:,2);Mu1(:,3);Mu1(:,4);Mu1(:,5)];
wealth=repmat(a',[nz 1]);
figure(2)
c_w=gini(pop, wealth,true);
title(['Wealth, gini=',num2str(c_w)])
M=sum(Mu)
figure(3)
plot(a,M)
title(['Distribution of wealth'])
Y=repmat(y_s',[1,num_a]);
A=repmat(a,[nz,1])
c=Y+r*A-pol_fn;
cf=c(:,pol_indx');
cf1=reshape(cf,[nz num_a nz]);
i=1;
while i < nz+1
c1(i,:)=PI(i,:)*cf1(:,:,i);
i=i+1;
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta*c1.^(-sigma)*r).*Mu))