function max_error = second_order_ode1
%% clear paths and call BeBOT from folder
clear
clc
addpath('BeBOT')
%% set order
K = 50; % knots (including terminals)
m = 3; % poly order between knots
L = 3;
% when changing m, also change return variable Z to enforce appropriate
% continuity constraints
[tnodes,w,Dm] = PiecewiseBeBOT(m,linspace(0,L,K));
%% initil condition
x0 = 1; 
x_dot0 = 0;
%% initial guess control points
for k = 1:K-1
    X_guess(:,k)=linspace(0,0,m+1);
end
X_guess = X_guess(:)';
%% fsolve call
options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','FunctionTolerance',1e-10);
fun = @(X)eqn(X);
tic
X = fsolve(fun,X_guess,options);
toc
X_dot = X*Dm;
%% plot
figure(1)
PBP = PiecewiseBernsteinPoly(X,linspace(0,L,K),linspace(0,L,100));
plot(linspace(0,L,100),1./sqrt(1+(linspace(0,L,100)).^2./3))
hold on
plot(linspace(0,L,100),PBP);
plot(tnodes,X,'o')
title('CBP Collocation')
legend('True Solution','CBP','Control Points')
%% error plot
figure(2)
error = (PBP-(1./sqrt(1+(linspace(0,L,100)).^2./3)));
max_error = max(error);
plot(linspace(0,L,100),error)
title('abs error')
%% collocation function x'(s)=sin(s), x(0)=0
    function Z = eqn(X)
        %% derivatives
        X_deriv = X;
        derivs = X;
        for i = 1:m-1
            X_deriv = X_deriv*Dm;
            derivs(end+1:end+k*(m+1)) = X_deriv;
        end
        X_dot = X*Dm;
        X_ddot = X_dot*Dm;
        %% get knots
        for i = 1:m
            [X_knots,knots_left0,knots_right0] = getKnots(derivs(1+(i-1)*(k*(m+1)):i*(k*(m+1))));
            if i == 1
                knots_left = knots_left0;
                knots_right = knots_right0;
            else 
                knots_left(end+1:end+k-1) = knots_left0;
                knots_right(end+1:end+k-1) = knots_right0;
            end
        end
        [X_knots,knots_left0,knots_right0] = getKnots(X);
        [X_dot_knots,knots_left1,knots_right1] = getKnots(X_dot);
        [X_ddot_knots,knots_left1,knots_right1] = getKnots(X_ddot);
        a_knots = X_ddot_knots;
        v_knots = X_dot_knots;
        x_knots = X_knots;
        %% initial condition
        Z1 = X(1)-x0;
        Z2 = X_dot(1)-x_dot0;
        %% ODE at knots
        s = linspace(0,L,K);
        s(1) = 0.001;
        Z3 = a_knots + (2./s).*v_knots + x_knots.*x_knots.*x_knots.*x_knots.*x_knots;
        %% Natural Spline
        ZN1 = X_ddot(1);
        ZN2 = 0;X_ddot(end);
        %% Continuity
        ZC = knots_left-knots_right;
        %% return variable
        Z = [ZN1,ZN2,ZC,Z1,Z2,Z3];
    end
%% pull knot values
    function [knots,knots_left,knots_right] = getKnots(Cp)
        Cp_reshape = reshape(Cp,m+1,K-1);
        knots = [Cp(1), Cp_reshape(end,1:end)];
        knots_left = Cp_reshape(end,1:end-1);
        knots_right = Cp_reshape(1,2:end);
    end
end