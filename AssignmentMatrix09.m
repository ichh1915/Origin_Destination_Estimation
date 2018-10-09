function [ Ao,Aod,Aod_path,EList,OList,ODList,V,Frac ] = AssignmentMatrix09( A,tau_max )
%% Revised Elastic with traffic travel at 0.5 link/unit time
% Fixed a link, a proportion of the vehicles has speed of 2_links/unit_time while
% the other proportion has speed of 1_link/unit_time

% Generate random traffic assignment matrices 
% Input: 
%   A: Adjacent matrix
%   tau_max: maximum number of steps
% Output: 
%   Ao: O-Flow traffic assignment matrix
%       In Array form
%   Aod: OD-Flow traffic assignemnt matrix
%       In Array form

% Initialization
[~,PathM,ODList] = PathGenerate01( A,tau_max );

[PathN,~] = size(PathM); %number of all possible paths that take maximum 4 steps=252
[ODN,~] = size(ODList); %number of OD pairs=72
OList = union(ODList(:,1),[]); %list of all origins 1-9
OListLen = length(OList); %=9

PathNList = ODList(:,4) - ODList(:,3) + ones(ODN,1);
path_max = max(PathNList); %7

EList = AdjacentM2Edges01( A ); %list of all links
[EdgeN,~] = size(EList); %=24


%% Constructing variables for Elastic Speed Model
vN_max = 10;  %maximum number of different speeds
vN_min = 1;   %minimum number of different speeds
vN = round((vN_max-vN_min)*rand+vN_min);  %Number of different speeds 
                                          %/Single Discrete uniformly distributed random number range from 1 to vN_max
                                          
v_max = 4;  %maximum possible speed
v_min = 0.2;    %minimum possible speed
V = zeros(vN,1);
for i=1:vN
    V(i,1) = (v_max-v_min)*rand+v_min;  %Values of speeds /Continuous uniformly distributed random number range from 0 to 4
end

Frac = randfixedsum(vN,1,1,0,1);  %Number of diff. fractions = Number of diff. speeds
                                  %Sum of the fractions  = 1
t_max = ceil(tau_max/min(V));     %maximum amount of time taken to reach any destination                         
                                  

Ao = zeros(EdgeN,OListLen,t_max); %[24*9*t_max]
Aod = zeros(EdgeN,ODN,t_max); %[24*72*t_max]
Aod_path = zeros(EdgeN,ODN,t_max,path_max);  %[24*72*t_max*7] path_max = maximum number of possible paths for OD pairs

PathProb = 1 + rand(PathN,1)*10; %[252 by 1] matrix of uniform random vairables from 1 to 11
PathProb = PathProb/sum(PathProb); %sum(PahtProb) = 1

ODProb = zeros(ODN,1); %[72 by 1]
for odn = 1:ODN
    ODProb(odn) = sum( PathProb(ODList(odn,3):ODList(odn,4)) ); % [72 by 1] matrix of ODProb
end
%sum(ODProb) = 1


%% Construct OD traffic Assignment Matrix                              
for odn = 1:ODN %OD pair index
    
    for pn = ODList(odn,3):ODList(odn,4) %Pathes index for OD pair odn
        offset = ODList(odn,3)-1;
        p_index = pn-offset;   %Eg: pn=3:6  =>  p_index=1:4  
        stepN = sum(PathM(pn,:)~=0)-1; %Number of steps for path pn
        
        for tn = 1:t_max %Time index 
            %Round each distance to the nearest integer
            %less or equal to that distance
            %Vehicle sensor located at the end of each link:
            %(Eg:Traffic at distance 1.5 links is registered as traffic on link one)
            for vn = 1:vN %Speed index
                if(V(vn,1)*tn<=stepN)
                    if( V(vn,1)*tn>=1 )
                    Edgen = find( EList(:,1)==PathM(pn,floor(V(vn,1)*tn)) & ...     
                        EList(:,2)==PathM(pn,floor(V(vn,1)*tn)+1) ); 
                    Aod_path(Edgen,odn,tn,p_index) = Aod_path(Edgen,odn,tn,p_index)+ ...
                        Frac(vn,1)*PathProb(pn)/ODProb(odn); 
                    end
                end          
            end 
            
        end
        
    end
    
end

for p_i=1:path_max
    Aod(:,:,:) = Aod(:,:,:) + Aod_path(:,:,:,p_i);
end

%% Construct O-Flow traffic Assignment Matrix
OProb = zeros(OListLen,1);
for oNode = 1:OListLen
    o_ODList = find( ODList(:,1)==OList(oNode) ); %List of OD pairs that originate from 'o'
    OProb(oNode) = sum(ODProb(o_ODList)); 
    for odn = o_ODList
        for tn = 1:t_max
            Ao(:,oNode,tn) = Ao(:,oNode,tn) + ...
                Aod(:,odn,tn)*ODProb(odn)/OProb(oNode); 
        end
    end
end

end
