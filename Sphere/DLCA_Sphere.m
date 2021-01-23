close all
clear
clc
set(0,'DefaultLineLineWidth',2)
% DLCA simulation made by SAA
% randa is a random number generator function that I made  @ SAA
% move is a function that I made @ SAA
% Center maker is a function I made @ SAA
% calcRg is a function made by Dr. Ganesan @ SAA
% We use a periodic boundary condition so to calculate Rg we are mapping everything
% out of the box and calculate it relative to every other cluster since it
% only depends on the center of mass and orientation
%% making N particles
N=300; % choose the number of particles for your system ( we can make this as an input)
PAR=cell(1,N); % Making a cell with N particles in it each cell being a sinlge particle initially
PAR_C=ones(1,N); % number of partiles in each cell to track how many particles are in each cells and we also use it in the algorithm
ITR=75000; %input('Number of iterations: '); % Monte carlo, the number of iterations that the code will run
r=0.5; %input('Radius of Particle: '); % The radius of the particles
StepSize = 1/6; % number of steps that a particle takes before it reaches the next position
%% Making the lattice
L = 25; % size of lattice is L*L*L
MIN = [0 0 0];   %lower bound   % x y z min, don't change this
MAX = [L L L];   %upper bound   % x y z max you can set to any
DIM=length(MIN);   % dimensions of the solution region
[PAR, POS] = CenterMaker(0,L,N,DIM,r,PAR); % Making N number of nonoverlapping particles between 0 and L
PAR_Rg = PAR; % PAR_Rg is used to calculate the radius of gyration read the SOP to understand why we are doing this
filename = 'choose a name.gif'; % gif stuff
h = figure;
PAR3 = cell(1,1); % This will be used to save data at each iteration
Rg3 = cell(1,1);
mobility3 = cell(1,1);
PAR_Rg3 = cell(1,1);
Update_PAR=0;
[x,y,z]=sphere(50);  % MATLAB sphere maker with (0,0,0) center
tic
itr = 0;
NLE = 1; % not large enough, to check if we ended up with one cluster or if our Rg is high enough to stop the simulation
while  NLE == 1 && itr < ITR % &&  testing the termination condition Rg/4
    itr = itr + 1  % you can put a ; if you want but this shows you how many iterations have passed and also used later for the plot
    %% Particle velocity maker
    phi = randa(-pi,pi,N,1); % each particle moves randomly a distance equal to 1/3 of the radius in any direction
    theta = randa(-pi,pi,N,1);
    [VEL(:,1),VEL(:,2),VEL(:,3)] = sph2cart(phi,theta,randa(0,r/3,N,1));
    %% Checking for clusters and how many particles in each cluster
    Stuck = 1; % for the second while loop
    Movement = 1; % for the first while loop
    while Movement == 1 % to make sure that we go through all of the movement steps that we set
        while Stuck == 1 % if a cluster has formed we need to make sure that we are checking all of the particles each time and not skipping any
            Stuck = 0; % while stuff
            for p1 = 1:N % going through particles 1 by one and the one next to them in the line below
                for p2 = p1+1:N
                    if PAR_C(p1) ~= 0 && PAR_C(p2) ~= 0  % if PAR_C is at least 1 there is a cluster to work with if not we don't
                        L1 = size(PAR{p1},1);
                        for ss = 1:L1
                            L2=size(PAR{p2},1);
                            for ss2 = 1:L2
                                if vecnorm(PAR{p1}(ss,:) - PAR{p2}(ss2,:)) <= 2*r % can be changed with pdist2, regardless we see if two particles collide
                                    diff_p1 = PAR_Rg{p1}(ss,:)-PAR{p1}(ss,:); % find the difference between mapped and unmapped
                                    diff_p2 = PAR_Rg{p2}(ss2,:)-PAR{p2}(ss2,:);
                                    if isequal(PAR{p1},PAR_Rg{p1}) && isequal(PAR{p2},PAR_Rg{p2}) % all of the particles are inside the boundary (periodic BC did not apply)
                                        PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                        PAR_Rg{p2} = [];
                                        Update_PAR = 1;
                                    elseif ~isequal(PAR{p1},PAR_Rg{p1}) && isequal(PAR{p2},PAR_Rg{p2}) % some part of p1 is out of the boundary and all of p2 is inside
                                        if isequal(PAR{p1}(ss,:),PAR_Rg{p1}(ss,:)) % collision happened with the part that is inside
                                            PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}]; % normal way of combining clusters
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        else % collision happened with the part out of the boundary
                                            PAR_Rg{p2} = PAR_Rg{p2}+(diff_p1.*ones(size(PAR_Rg{p2},1),3)); % shift the desired part outside
                                            PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        end
                                    elseif ~isequal(PAR{p2},PAR_Rg{p2}) && isequal(PAR{p1},PAR_Rg{p1})  % some part of p2 is out of the boundary and all of p1 is inside
                                        if isequal(PAR{p2}(ss2,:),PAR_Rg{p2}(ss2,:)) % collision happened with the part that is inside
                                             PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}]; % normal way of combining clusters
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        else % collision happened with the part out of the boundary
                                            PAR_Rg{p1} = PAR_Rg{p1} + (diff_p2.*ones(size(PAR_Rg{p1},1),3)); % Shift the desired part outside
                                            PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        end
                                    elseif ~isequal(PAR{p1},PAR_Rg{p1}) && ~isequal(PAR{p2},PAR_Rg{p2}) % both p1 and p2 have parts that are out of the boundary
                                        if isequal(PAR{p1}(ss,:),PAR_Rg{p1}(ss,:)) && isequal(PAR{p2}(ss2,:),PAR_Rg{p2}(ss2,:)) % if the parts that collided were all inside
                                            PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}]; % normal way
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        elseif ~isequal(PAR{p1}(ss,:),PAR_Rg{p1}(ss,:)) && isequal(PAR{p2}(ss2,:),PAR_Rg{p2}(ss2,:)) % the part of p1 is outside and p2 is inside
                                            PAR_Rg{p2} = PAR_Rg{p2}+(diff_p1.*ones(size(PAR_Rg{p2},1),3)); % shift the desired part outside
                                            PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        elseif isequal(PAR{p1}(ss,:),PAR_Rg{p1}(ss,:)) && ~isequal(PAR{p2}(ss2,:),PAR_Rg{p2}(ss2,:)) % the part of p1 is inside and p2 is outside
                                            PAR_Rg{p1} = PAR_Rg{p1} + (diff_p2.*ones(size(PAR_Rg{p1},1),3)); % Shift the desired part
                                            PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                            PAR_Rg{p2} = [];
                                            Update_PAR = 1;
                                        elseif ~isequal(PAR{p1}(ss,:),PAR_Rg{p1}(ss,:)) && ~isequal(PAR{p2}(ss2,:),PAR_Rg{p2}(ss2,:)) % both parts are outside
                                            if dot(diff_p1,diff_p2) == 0 % the parts are in 2 different faces
                                                PAR_Rg{p2} = PAR_Rg{p2}-(diff_p2.*ones(size(PAR_Rg{p2},1),3))+(diff_p1.*ones(size(PAR_Rg{p2},1),3)); % shift the desired part outside
                                                PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                                PAR_Rg{p2} = [];
                                                Update_PAR = 1;
                                            else % the parts are in the same face
                                                PAR_Rg{p1} = [PAR_Rg{p1};PAR_Rg{p2}];
                                                PAR_Rg{p2} = [];
                                                Update_PAR = 1;
                                            end
                                        end
                                    end
                                    if Update_PAR == 1
                                        PAR{p1} = [PAR{p1};PAR{p2}]; % sticking the particles to each other if they touch (unmapped, actual)
                                        PAR{p2} = []; % if a particles sticks to another particles we set the cell that was holding it to empty
                                        PAR_C(p1) = PAR_C(p1) + PAR_C(p2); % keeping track of how many particles we have in each cell (cluster)
                                        PAR_C(p2) = 0; % if a particle sticks to another cluster we set it's place holder to zero
                                        Update_PAR = 0;
                                        Stuck=1; % at least 1 particle is attached to a cluster so we need to redo the whole process to make sure we are not skipping what was adjacent to that particle
                                        break
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        [PAR_Rg_Box] = Box_Cut(PAR,N,L,r); % make a box and see what clusters are inside it
        PAR_CMin = min(nonzeros(PAR_C'),[],'all');
        PAR_CMax = max(nonzeros(PAR_C'),[],'all');
        mobility = PAR_CMin./PAR_C'; % probability of movement is minimum number of particels over number of particles in the selected cluster
        VEL = mobility .* VEL; % probability of displacement
        [PAR,Movement,VEL] = Move_Rg(VEL,PAR,N,StepSize,MAX,MIN); % function made to move the particle around
        PAR3{itr} = PAR; % keeping data from each iteration
        PAR_Rg3{itr} = PAR_Rg;
        mobility3{itr} = mobility;
        if  nnz(PAR_C) == 1 %|| max(Rg) >= L/4
            NLE = 0   % if prticles end up in one cluster or get too big stop the simulation
            disp('Change ITR you ended up with one cluster')
        end
    end
end
toc
tic
for j=0:100:itr % change the steps here to the get desired video frames, you can play around with it to get the frames that you want
    if j == 0
        ii = 1; % making sure loop matches if you go by increments other than 1
    else
        ii = j;
    end
    PAR2 = [];   % making an empty matrix to use for plotting later
    for i= 1:size(PAR3{ii},2)
        if size(PAR3{ii}(i),1) ~= 0
            PAR2 = [PAR2;PAR3{ii}{i}]; % putting the particles in a matrix to plot them
        end
    end
    for i=1:size(PAR2,1)
        X=r*x+PAR2(i,1);
        Y=r*y+PAR2(i,2);
        Z=r*z+PAR2(i,3);
        b=surf(X,Y,Z,'FaceColor',[0 0 1], ...
            'FaceAlpha',0.5,'FaceLighting','none','EdgeColor','none') ;
        hold on
        %set(b,'FaceColor',[0 0 1],'FaceAlpha',0.5,'FaceLighting','gouraud','EdgeColor','none')
    end
    hold off
    grid on % plot stuff
    box on % plot stuff
    axis square
    title(['iteration=',num2str(ii)]); % gif stuff
    axis([0 L 0 L 0 L]);
    %     xticks(XMIN:5:XMAX);
    %     yticks(YMIN:5:YMAX);
    %     zticks(ZMIN:5:ZMAX);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    drawnow;
    frame = getframe(h); % gif stuff
    im = frame2im(frame); % gif stuff
    [imind,cm] = rgb2ind(im,256); % gif stuff
    % Write to the GIF File
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
toc
%save('NN(put number)_r(put radius)_10k.mat','mobility','PAR','PAR3','PAR_C','POS','mobility3','PAR_Rg','PAR_Rg_Box','PAR_Rg3') use this
%line to save the specs that are important
