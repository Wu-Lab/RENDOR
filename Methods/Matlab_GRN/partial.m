function [net_new]=partial(W_input)
n = length(W_input);
net_new = -diag(1./diag(W_input))*W_input*diag(1./diag(W_input));
net_new = net_new - diag(diag(net_new)) + eye(n);
