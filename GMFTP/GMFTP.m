function theta_star = GMFTP(A, F, lambda, K, repeat_times, T, rho, tau)
% The core function of this paper. It detects protein complexes from a
% given PPI network with adjacent matrix A and a functional profile with
% association matrix F.

% Input:
%   A: adjacent matrix of PPI network with size N by N.
%   F: association matrix of functional profile with size N by C.
%   lambda: rate parameter of exponential distribution. The default value
%   is 4.
%   K: maximum number of possible protein complexes. The default value is
%   1000.
%   repeat_times: the number of times that we repeat the entire calculation
%   to avoid local minimum. The default value is 100.
%   T: the number of iterations limited in GMFTP. The default value is 400.
%   rho: the tolerance threshold of the stop criterion. The default
%   value is 1e-6.
%	tau: the threshold used to obtain protein complexes from estimator of theta. The default
%   value is 0.2.


% Outputs:
%   theta_star: the protein-complex membership indication matrix.


if nargin < 8
    tau = 0.2;
end

if nargin <7
    rho = 1e-6;
end

if nargin < 6
    T = 400;
end

if nargin < 5
    repeat_times = 100;
end

if nargin < 4
    K = 1000;
end



if nargin < 3
    lambda = 4;
end

if nargin < 2
    error('You need input A and F');
end

% parameter estimation
[theta, psi, score] = GMFTP_ParameterEstimation(A, F, lambda, K, repeat_times, T, rho);

% protein complex detection
theta_star= GMFTP_ComplexDetection(A, theta, tau);








function [theta, psi, score] = GMFTP_ParameterEstimation(A, F, lambda, K, repeat_times, T, rho)
% Estimate model parameters theta and psi by updating them
% according to Equations (11) and (12) alternately. It repeats the entire
% updating procedure multiple times and chooses the result that gives
% the lowest value of objective function in Equation (10) as the final estimator.


% Input:
%   A: adjacent matrix of PPI network with size N by N.
%   F: association matrix of functional profile with size N by C.
%   lambda: rate parameter of exponential distribution. The default value
%   is 4.
%   K: maximum number of possible protein complexes.  The default value is
%   1000.
%   repeat_times: the number of times that we repeat the entire calculation
%   for parameter estimation to avoid local minimum. The default value is 100.
%   T: the number of iterations limited in GMFTP. The default value is 400.
%   rho: the tolerance threshold of the stop criterion.  The default
%   value is 1e-6.


% Outputs:
%   theta: the estimator of protein-complex affinity matrix.
%   psi: the estimator of complex-function preference matrix.
%   score: the value of objective function in Equation (10).

if nargin < 7
    rho = 1e-6;
end

if nargin < 6
    T = 400;
end

if nargin < 5
    repeat_times = 100;
end

if nargin < 4
    K = 1000;
end


if nargin < 3
    lambda = 4;
end

if nargin < 2
    error('You need input A and F');
end

fprintf('Estimating parameters...')
fprintf('\n')

N = size(A,1);
C = size(F,2);

% Introduce S to represent whether the corresponding proteins are functional
% characterized, where S(i) =1 means protein i is annotated, and S(i)=0
% otherwise.
S = sum(F,2) >0;

lowest_score = inf;

% Repeat the entire updating procedure multiple times.
for i = 1: repeat_times
    fprintf(['This is the ',num2str(i), '-th repeat...'])
    fprintf('\n')
    theta = rand(N,K); %Initialize  theta randomly.
     psi = rand(K,C); %Initialize  psi randomly.
    score_old = inf;
    
    for j = 1 : T
        % Update theta according to Equation (11).
        theta = theta .* ( (  diag(S)*((F./( 1-exp(-theta* psi)+eps ) )* psi') + A./( 1-exp(-theta*theta') +eps)*theta ) ./ ( diag(S)*repmat(sum( psi,2)',N,1) + repmat(sum(theta),N,1) + lambda  + eps ) );
        % Update psi according to Equation (12).
        psi =  psi .* ( (theta'*diag(S)*(F./( 1-exp(-theta* psi)+eps))) ./ ( repmat(sum(diag(S)*theta)',1,C) + lambda + eps ) );
        
        
        % Calculate the value of the objective function (10).
        score =  - sum(sum(diag(S)*(F.*log(1-exp(-theta*psi)+eps)))) + sum(sum(diag(S)*((1-F).*(theta*psi)))) + 0.5*( -sum(sum(A.*log(1-exp(-theta*theta')+eps))) + sum(sum((1-A).*(theta*theta')) ))  + lambda*sum(sum(theta)) ...
            + lambda*sum(sum(psi)) ;
        
        
        if abs(score-score_old) / abs(score) < rho
            break;
        else
            score_old = score;
        end
        theta(theta<eps) = 0;
        psi(psi<eps) = 0;
    end
    
    
    
    % Choose the results that give the lowest value of objective
    % function.
    if score < lowest_score
        final_theta = theta;
        final_psi = psi;
        lowest_score = score;
    end
    
end

theta = final_theta;
psi = final_psi;
score = lowest_score;


function theta_star = GMFTP_ComplexDetection(A, theta, tau)
% Detect protein complexes from the estimators of model parameters.

% Obtain theta_star according to Equation (13)
fprintf('Detecting promplex detection...')
fprintf('\n')

theta_star = theta;

% Normalize the rows of theta
theta_star = theta_star ./ repmat(sum(theta_star,2)+eps,1,size(theta_star,2));

theta_star(theta_star >= tau) = 1;
theta_star(theta_star < tau) = 0;

% If the detected complex candidate is composed of several isolated subnetworks, process it such that
% each connected subnetwork is a complex.
[N, K] = size(theta_star);
t = 1;
for k = 1:K
    k_index = find(theta_star(:,k));
    [S, C] = graphconncomp(sparse(A(k_index,k_index)));
    for i = 1:S
        temp_theta_star(:,t) = zeros(N,1);
        temp_theta_star(k_index(C==i),t) = theta_star(k_index(C==i),k);
        t = t+1;
    end
end
theta_star = temp_theta_star;

% Filter out the detected complexes which include less than three proteins.
small_size_indices = sum(theta_star)<=2;
theta_star = theta_star(:,~small_size_indices);

