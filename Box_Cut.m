function [PAR_Rg] = Box_Cut(PAR,N,XMAX,r)
Thold = XMAX - 1;
for i = 1:N
    if sum(sum(PAR{i}<=Thold,2)==3)==size(PAR{i},1)  % only get out the clusters that are completely inside the box
        PAR_Rg{i} = PAR{i};
    elseif isempty(PAR{i})
        PAR_Rg{i} = [];
    else
        PAR_Rg{i} = [];
    end
end

