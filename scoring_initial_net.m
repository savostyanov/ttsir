function [ind] = scoring_initial_net(Xobs)
% Initial guess of the network based on probability of 10 -> 11 transitions
[M,d,nt] = size(Xobs);

scores = zeros(d,d);
for j=1:M
    state = Xobs(j,:,1);
    for i=2:nt
        if any(Xobs(j,:,i)~=state)
            iprev1 = find(state==1);
            inext1 = find((state==0) & (Xobs(j,:,i)==1));
            scores(iprev1, inext1) = scores(iprev1, inext1) + 1/numel(iprev1);
            state = Xobs(j,:,i);
        end
    end
end
scores = scores + scores';

cnt = 1;
ind = zeros(1, d*(d-1)/2);
for i=1:d
    for j=i+1:d
        ind(cnt) = scores(i,j);
        cnt = cnt+1;
    end
end
scores = ind;
scores = scores/sum(scores);
ind = double(scores > 2/(d*(d-1))) + 1;
end
