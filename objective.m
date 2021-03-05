
function [Y]=objective(X)

disp(X);

scale=[7.00  8.83886  10.27  16.95 0.1   0.1   0.1   0.1   0.1   0.1     0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1   0.1     0.1   0.1   0.1   0.1   0.1   0.1  0.1   0.1   0.1   0.1   0.1   0.1     0.1   0.1   0.1   0.1   0.1   0.1 10 10 10 16 89000 24600];

X=X.*scale;

R = 5e6;            % range
C_t = 1.8639e-4;    % spec. fuel consumption
V = 0.76*295.2;     % cruise velocity

CL_CD = X(44);
W_TO_max = X(45);

Wfuel = (1-0.938*(1/(exp(R*C_t/(V*CL_CD)))))*W_TO_max;

Y = Wfuel/24600;

toc
end
