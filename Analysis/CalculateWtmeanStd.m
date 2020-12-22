function [wtInfo] = CalculateWtmeanStd(invec,wtvec)

n=length(wtvec);
wtmean = sum(invec.*wtvec) ./ sum(wtvec);

% wtstd = sqrt( ( sum(wtvec.*((invec-wtmean).^2)) ) / ....
%      ( ((length(wtvec)-1)/length(wtvec)) * sum(wtvec)) );

wtstd = std(invec,wtvec);

wtInfo = [wtmean,wtstd,n];
end
