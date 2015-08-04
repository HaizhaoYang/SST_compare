function Dist = EMDMat(showMat,InstFreq,R_high)
[sm sn] = size(showMat);
Dist = 0;
for cnt = 1:sn
    distA = showMat(:,cnt);
    distA = distA/sum(distA);
    distB = zeros(sm,1);
    distB(round(sm*InstFreq(cnt)/R_high)) = 1;
    Dist = Dist + distVec(distA,distB)/sn;
end