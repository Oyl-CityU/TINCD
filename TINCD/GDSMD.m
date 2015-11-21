function [theta_star, score] = GDSMD(W, K, lambda, n_iter,repeat_times)
% GDSMD is the main algorithm of Graph regularized Doubly Stochastic Matrix Decomposition (GDSMD) model.  
%
% Inputs:
%   W: integrated similarity matrix. W is a N*N matrix.
%   K: maximum number of possible protein complexes.The default value is
%   2000.
%   lambda: rate parameter of exponential distribution. The default value
%   is 32.
%   n_iter: the number of iterations limited in GDSMD. The default value is
%   200.
%   repeat_times: the number of times that we repeat the entire calculation
%   to avoid local minimum. The default value is 20.
% Outputs:
%   theta_star: resultant protein-complex indication matrix.
%   score: the value of objective function in Equation (8).

    if nargin < 5
        repeat_times = 20;
    end

    if nargin < 4
        n_iter = 200;
    end

    if nargin < 3
        lambda = 32;
    end
        
    if nargin < 2
        error('Yor need input integrated similarity matrix W and maximum possible number K of protein complexes.');
    end
    
    n = size(W,1);
    W = W - diag(diag(W));
    for i = 1 : n
        if sum(W(i,:)) == 0
            W(i,i) = 1;
        end
    end
    W = diag(sum(W,2))*W;
    W = 0.5*W + 0.5*eye(n);
    W = (W + W')/2;
    
    lowest_score = inf;
    for i = 1: repeat_times
    fprintf(['This is the ',num2str(i), '-th repeat...'])
    fprintf('\n')
    H_old = rand(n,K);
    D = diag(sum(W));
	score0 = inf;
    for j  = 1: n_iter
        S = diag(1./(sum(H_old)+eps));
        U = H_old*S*H_old';
        Z = W./U;
        Dp = repmat(diag(H_old'*Z*H_old)',n,1)*(S.^2)  + 2*repmat(sum(H_old),n,1)*S  +  2*lambda*(D*H_old);
        Dn = 2*(Z*H_old)*S  + repmat(diag(H_old'*ones(n,n)*H_old)',n,1)*(S.^2) + 2*lambda*(W*H_old); 
        a = diag(sum(H_old./(Dp + eps),2));
        b = sum(H_old.*(Dn./(Dp + eps)),2);
        H = H_old.*((a*Dn + 1)./(a*Dp + repmat(b,1,K) + eps));
		
		% Calculate the  value of the objective function (9).
        score = - sum(sum(W.*log(U + eps))) + sum(sum(U)) + lambda*sum(sum(trace(H'*D*H) - trace(H'*W*H)));
		if  abs(score0 - score)/score< 1e-6
            break;
        else        
            H_old = H;
			score0 = score;
        end
    end
    

    if score < lowest_score
        final_H = H;
        lowest_score = score;
    end
    end
 
    fprintf('Detecting protein complex ...')
    fprintf('\n')

    theta = final_H;
    theta(theta < 10^(-3)) = 0;
    [theta_sort,~] = sort(theta,2,'descend');
    dif_theta = theta_sort(:,1:(size(theta_sort,2) - 1)) - theta_sort(:,2:size(theta_sort,2));
    [~,dif_ID] = max(dif_theta,[],2);
    act_val = zeros(n,1);
    for i = 1 : n
       act_val(i) = theta_sort(i,dif_ID(i) + 1);
    end
    act_matrix = repmat(act_val,1,size(theta,2));
    theta_star = zeros(n,size(theta,2));
    theta_star(theta > act_matrix) = 1;
    theta_star(:,sum(theta_star) <= 2) = [];

    
