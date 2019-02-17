function ASE = ASE(clear,recovered)
% This function calculates the Average Squared Error of two input signals
% Input: 
% clear - clear signal, that has no noise
% recovered - signal that has gone through a noise reduction process
% mu - µ the step-size
% Output:
% ASE - Average Squared Error
Ns = length(clear);
cc = clear(Ns/2:Ns);
ce = cc-recovered(Ns/2:Ns);
ASE = (ce'*ce)/(cc'*cc);

% upper = clear(end/2:end) - recovered(end/2:end);
% lower = clear(end/2:end)'*clear(end/2:end);
% ASE = (upper'*upper)/lower;
% upper = mean(mean(clear(end/2+1:end)-recovered(end/2+1:end))).^2;
% lower = mean(mean(clear(end/2+1:end))).^2;
% ASE = upper/lower;

end

