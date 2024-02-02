% Convert an adjacency matrix to the network configuration index
function [ind] = adj_to_ind(W)
d = size(W,1);
cnt = 1;
ind = ones(1, d*(d-1)/2);
for i=1:d
    for j=i+1:d
        ind(cnt) = double(W(i,j)>0)+1;
        cnt = cnt+1;
    end
end
end
