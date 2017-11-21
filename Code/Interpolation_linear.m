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
    while v_tol >.00001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        for i=1:nz;
        E(:,:,i)=repmat(PI(i,:)*v_guess, [num_a 1]);
        end
        v=ret+beta*E;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vmax, pol_indx]=max(v, [], 2);
        vmax=permute(vmax, [3 1 2]);
        %pol_indx=permute(pol_indx, [3 1 2]);
        pol_indx_1=pol_indx-1;
        pol_indx_1(pol_indx_1==0)=1;
        pol_indx_2=pol_indx+1;
        pol_indx_2(pol_indx_2==501)=500;
        %pol_indx_1=reshape(pol_indx_1', [num_a, 1, nz]);
        %pol_indx_2=reshape(pol_indx_2', [num_a, 1, nz]);
        %pol_indx_m=reshape(pol_indx', [num_a, 1, nz]);
        a_new=[a(pol_indx_1),a(pol_indx),a(pol_indx_2)];
        i=1;
        while i<nz+1
            j=1;
            while j<num_a+1
              pol_indx_1(j,:,i)=pol_indx_1(j,:,i)+500*(j-1)+250000*(i-1);
              pol_indx_2(j,:,i)=pol_indx_2(j,:,i)+500*(j-1)+250000*(i-1);
              pol_indx_m(j,:,i)=pol_indx(j,:,i)+500*(j-1)+250000*(i-1);
              j=j+1;
            end
        ret_t(:,:,i)=ret(:,:,i)';
        E_t(:,:,i)=E(:,:,i)';
        i=i+1;
        end
        ret_1=ret_t(pol_indx_1);
        E_1=E_t(pol_indx_1);
        ret_2=ret_t(pol_indx_2);
        E_2=E_t(pol_indx_2);
        ret_m=ret_t(pol_indx_m);
        E_m=E_t(pol_indx_m);
        ret_new=[ret_1,ret_m,ret_2];
        E_new=[E_1,E_m,E_2];
        b=[1,2,3];
        b1=interpn(b,5);
        i=1;
        while i < nz+1
           ret_new1(:,:,i)=(interp1(b,ret_new(:,:,i)',b1))';
           E_new1(:,:,i)=(interp1(b,E_new(:,:,i)',b1))';
           a_new1(:,:,i)=(interp1(b,a_new(:,:,i)',b1))';
           i=i+1;
        end
        v_new=ret_new1+beta*E_new1;
        [vmax_new, pol_indx_new]=max(v_new, [], 2);
        vmax_new=permute(vmax_new, [3 1 2]);
        
        
        v_tol = max(max(abs(vmax_new-v_guess)));
        v_guess=vmax_new;
    end;
     % KEEP DECSISION RULE
    pol_indx_new=permute(pol_indx_new, [3 1 2]);
        i=1;
        while i<nz+1
            j=1;
            while j<num_a+1
              pol_indx_new1(i,j)=pol_indx_new(i,j)+65*(j-1)+500*65*(i-1);
              j=j+1;
            end
        a_new_t(:,:,i)=a_new1(:,:,i)';
        i=i+1;
        end
    
     pol_fn=a_new_t(pol_indx_new1);
     i=1; 
     while i <nz+1
     d(:,:,i) = abs(bsxfun(@minus, pol_fn(i,:),a'));
     [a_min, a_indx]=min(d(:,:,i) ,[],1);
     p(i,:)=a_indx;
     di(i,:)=pol_fn(i,:)-a(a_indx);
     i=i+1;
     end
   P= di./(a(2)-a(1))+p;
   
   
   % SET UP INITITAL DISTRIBUTION
    %Mu=ones(nz, num_a)/(nz*num_a);
    Mu = zeros(nz,num_a);
  Mu(1, 4) = 1; % initial guess: everyone employed, 0 assets
 % Mu = ones(nz, num_a); %alternative initial guess: same mass in all states
%Mu = Mu / sum(Mu(:)); % normalize total mass to 1

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1; 
while mu_tol > 1e-07
    [emp_ind, a_ind] = find(Mu > 0); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = P(emp_ind(ii), a_ind(ii)); 
        
         
        if apr_ind<500
            if apr_ind>1
               apr_ind1=fix(apr_ind);
               apr_ind2=apr_ind1+1;
               MuNew(:, apr_ind1) = MuNew(:, apr_ind1) + ...
                   (apr_ind2-apr_ind)*(PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
               MuNew(:, apr_ind2) = MuNew(:, apr_ind2) + ...
                   +(apr_ind-apr_ind1)*(PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
            else
            apr_ind1=fix(apr_ind)+1;
            MuNew(:, apr_ind1) = MuNew(:, apr_ind1) + ...
                (PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
            end
        else
            apr_ind1=fix(apr_ind);
            MuNew(:, apr_ind1) = MuNew(:, apr_ind1) + ...
                (PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
        end
    end

    mu_tol = max(max(abs(MuNew-Mu)));
    
    Mu = MuNew ;
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
cf=c(:,P');
cf1=reshape(cf,[nz num_a nz]);
i=1;
while i < nz+1
c1(i,:)=PI(i,:)*cf1(:,:,i);
i=i+1;
end
Eulererror=sum(sum(abs(c.^(-sigma)-beta*c1.^(-sigma)*r).*Mu))