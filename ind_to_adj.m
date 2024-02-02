% Convert a network configuration index to adjacency matrix
function [W] = ind_to_adj(d, ind)
W = sparse(d,d);
cnt = 1;
for i=1:d
    for j=i+1:d
        W(i,j) = double(ind(cnt)==2);
        cnt = cnt+1;
    end
end
W = W+W';
end
