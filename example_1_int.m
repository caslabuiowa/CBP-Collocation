function max_error = second_order_ode1
%% interval 
L = 3;
K = 50;
M = 3;
%% precompute
T0 = int_transform(0,K,M,L);
T1 = int_transform(1,K,M,L);
T2 = int_transform(2,K,M,L);
T0_one = int_transform(0,K,M,L,0)*poly1_matrix(M,K,M);
T1_one = int_transform(1,K,M,L,0)*poly1_matrix(M-1,K,M);
T2_one = int_transform(2,K,M,L,0)*poly1_matrix(M-2,K,M);
bx = BernsteinMatrix_a2b(M,linspace(0,1,3)); bx = bx(1:2,:);
bxd = BernsteinMatrix_a2b(M-1,linspace(0,1,3)); bxd = bxd(1:2,:);
bxdd = BernsteinMatrix_a2b(M-2,linspace(0,1,3)); bxdd = bxdd(1:2,:);
s_coll = linspace(0,1,3)*L/(K-1); s_coll = s_coll(1:2)';
s_coll(1)=1e-7;
%% decision var
theta0 = zeros(1,K+M-1);
options = optimoptions('fsolve', ...
                       'Algorithm', 'levenberg-marquardt', ... % Change algorithm
                       'MaxIterations', 1000000, ... % Increase max iterations
                       'MaxFunctionEvaluations', 5000000, ... % Increase max function evaluations
                       'StepTolerance', 1e-8, ... % Adjust step tolerance
                       'FunctionTolerance', 1e-8); % Adjust function tolerancefun = @(X)eqn(X);
fun = @(eta)eqn(eta);
theta0 = fsolve(fun,theta0,options);
thetaM = theta0*int_transform(0,K,M,L,0);
%% plotting
figure(1)
splot = linspace(0,L,100);
PBP = PiecewiseBernsteinPoly(thetaM(1:end-M),linspace(0,L,K),splot);
plot(linspace(0,L,100),1./sqrt(1+(linspace(0,L,100)).^2./3))
hold on
plot(splot,PBP)
[tnodes,w,Dm] = PiecewiseBeBOT(M,linspace(0,L,K));
plot(tnodes,thetaM(1:end-M),'o')
title('CBP Collocation')
legend('True Solution','CBP','Control Points')
error = abs(PBP-(1./sqrt(1+(splot).^2./3)));
max_error = max(error);
figure(2)
plot(linspace(0,L,100),error)
title('abs error')
%% Diff Eq
    function Z = eqn(theta0)
        x = theta0*T0;
        xd = theta0*T1;
        xdd = theta0*T2;
        %% BC
        Z_BC1 = x(1) - 1;
        Z_BC2 = xd(1);
        Z_BC = [Z_BC1(:);Z_BC2(:)];
        %% eqn
        s = linspace(0,L,K);
        s(1) = 1e-7;
        Z_EQ = xdd + (2./s).*xd + x.*x.*x.*x.*x;
        %% Extra
        Z_EX = 0;
        if M>3
            x_coll = bx*(theta0*T0_one)';
            xd_coll = bxd*(theta0*T1_one)';
            xdd_coll = bxdd*(theta0*T2_one)';
            Z_EX = xdd_coll + (2./s_coll).*xd_coll + x_coll.*x_coll.*x_coll.*x_coll.*x_coll;
        end
        Z = [Z_BC;Z_EQ(:);Z_EX(:)];
    end
%% additonal funcitons
    function CPone = poly1_matrix(n,K,M)
        CPone = zeros((K-1)*(n+1)+M,n+1);
        CPone(1:n+1,1:n+1)=eye(n+1);
    end
end