function [PAR,POS] = CenterMaker(XMIN,XMAX,N,DIM,r,PAR)
POS = zeros(N,DIM);
POS(1,:) = randa(XMIN,XMAX,1,DIM); % initial position set randomly 
InCheck = 0;
r2 = (2*r)^2;
k=1;
while k <= N && InCheck <= 50000
    NewPosition = randa(XMIN,XMAX,1,DIM);
    d2 = sum((POS(1:k-1,:) - NewPosition(1,:)).^2,2);
    if sum(d2(1:end,1) > r2) == length(d2)*1
        POS(k,:) = NewPosition;
        k = k + 1;  
    end
    InCheck = InCheck + 1;
end
if InCheck >= 50000
             error('change your dimensions or your spheres it is taking to much time')
    end
for n=1:N
    PAR{1,n} = POS(n,:); % Turning the initial particle positions to a cell ( each cell is a particle initially )
end
end

