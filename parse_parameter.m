% If batch is empty, asks for keyboard input, replacing void by default.
% Otherwise returns batch
function [P] = parse_parameter(prompt, default, batch)
if (nargin<3)||(isempty(batch))
    P = input([prompt, sprintf(' (default %g): ', default)]);
    if (isempty(P))
        P = default;
    end
else
    P = batch;
end
end
