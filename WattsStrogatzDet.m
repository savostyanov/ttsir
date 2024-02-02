function h = WattsStrogatzDet(N,K,source,target)
% H = WattsStrogatzDet(N,K,source,target) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and deterministic rewiring of one
% edge (source->source+1) to (source->target)
%

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;
t(source,1) = target;

h = graph(s,t);
end
