function [PARout,PAR_Rgout,Movement,VELout] = Move_Rg(VEL,PAR,PAR_Rg,N,StepSize,MAX,MIN);
Movement = 0; % while stuff
VELout = VEL; % velocity ( new postion )
PARout = PAR; % Particle position
PAR_Rgout = PAR_Rg; % particle postiton with no boundary condition
for pos = 1:N % particle movement part 
    if isempty(PAR{pos})
       continue % if the cell is empty we move to the next one 
    elseif norm(VEL(pos,:)) >= StepSize 
        PARout{pos} = mod(PAR{pos}+(VEL(pos,:).*(StepSize/norm(VEL(pos,:)))),MAX);  % the moving of the particles, mod is used for the boundary conditions       
        if sum(sum(PARout{pos} == PAR{pos}+(VEL(pos,:).*(StepSize/norm(VEL(pos,:)))),2) == 3) == 0 % if all of the particle in a cluster are out of the boundary 
           PAR_Rgout{pos} = PARout{pos};
        else
           PAR_Rgout{pos} = PAR_Rg{pos} + (VEL(pos,:).*(StepSize/norm(VEL(pos,:))));
        end
        temp = (PAR{pos}(1,:)+VEL(pos,:))-PARout{pos}(1,:); % this part is used to determine how much of the distance between the old and new position is moved
        VELout(pos,:) = mod(abs(temp),MAX); % making sure that the remaining path is not out of the boundary
        VELout(pos,:)=VELout(pos,:) .* sign(temp); % new next velocity when we take the sign into account since we took the mod
        Movement = 1; % keep going with the while loop
    elseif norm(VEL(pos,:)) < StepSize && norm(VEL(pos,:)) ~= 0 % if the remaining distance is less than the stepsize
       PARout{pos} =  mod(PAR{pos}+VEL(pos,:),MAX); % new position, taking the BC into account 
       if sum(sum(PARout{pos} == PAR{pos}+VEL(pos,:),2) == 3) == 0 % if all of the particle in a cluster are out of the boundary 
           PAR_Rgout{pos} = PARout{pos};
       else
           PAR_Rgout{pos} = PAR_Rg{pos} + VEL(pos,:);
       end
       VELout(pos,:) = 0; 
       Movement = 1; % keep going with the while loop
    end
    
end
end
