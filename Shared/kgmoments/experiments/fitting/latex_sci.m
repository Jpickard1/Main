function strs=latex_sci(nums,k)
% Make the number nicely formatted for latex

strs = {};
if ~exist('k','var') || isempty(k), k=2; 
else k=abs(k); 
end
k = int2str(k);
for n=nums(:)'
    tmp=num2str(n,['%.' k 'e']);
    idx=find(tmp=='e');
    exnum=str2double(tmp(idx+1:end)); % exponent with no insignificant zeroes
    if exnum<0
        exnumstr = ['\mbox{-}' num2str(-exnum)];
    else
        exnumstr =  num2str(exnum);
    end
    numtxt=['$',tmp(1:idx-1),'\tdot 10^{',exnumstr,'}$'];
    strs{end+1} = numtxt;
end