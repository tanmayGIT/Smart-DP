%This code is for ESC algo
% Compute the self similarity join of time series


% For details of the algorithm, see:
% Mondal, T., Ragot, N., Ramel, J. Y., & Pal, U. (2016). Flexible Sequence 
% Matching technique: An effective learning-free approach for word spotting. 
% Pattern Recognition, 60, 596–612.


% Usage:
% [pathCost,pathTarget,indxcol,indxrow,distSum,jumpcost] = FSM_Algo(refSample,testSample,weight)

% Input:
%     refSample: The reference signal (or vector) to be matched
%     testSample: The target signal (or vector) with whom the reference signal would be matched
%     weight: The additional weight to be applied with the jumpcost

% Output:
%     pathCost: The path cost matrix
%     indxcol: The column or the indexes of matched target elements
%     indxrow: The rows or the indexes of matched query elements
%     distSum: The total distance value
%     jumpcost: The cost value accumulated due to number of jumps either in
%     query or in target
%%


function [pathCost,indxcol,indxrow,distSum,jumpcost] = FSM_Algo(refSample,testSample,weight)

[noOfSamplesInRefSample,N] = size(refSample);
[noOfSamplesInTestSample,M] = size(testSample);

if(noOfSamplesInRefSample == 0)
    disp('This is unwanted/unknown error');
end

if(N == M)
    Dist = zeros(noOfSamplesInRefSample,noOfSamplesInTestSample);
    pathCost = inf(noOfSamplesInRefSample,noOfSamplesInTestSample);
    pathTargetRw = zeros(noOfSamplesInRefSample,noOfSamplesInTestSample);
    pathTargetCol = zeros(noOfSamplesInRefSample,noOfSamplesInTestSample);
    
    for i=1:noOfSamplesInRefSample
        for j=1:noOfSamplesInTestSample
            total = zeros(N,1);
            for goFeature = 1:N
                total(goFeature,1) = (double((refSample(i,goFeature)-testSample(j,goFeature))^2));
            end
            Dist(i,j) = sqrt(sum(total));
        end
    end
    elasticity = (noOfSamplesInTestSample - noOfSamplesInRefSample);
    if(noOfSamplesInTestSample == noOfSamplesInRefSample)
        elasticity = noOfSamplesInRefSample;
    end
    
    if (~exist('jumpcost','var'))
        statmatx = calcJumpCost(Dist, elasticity);
        %         statmatx = min(Dist,[],2);
        [~,~,statmatx] = find(statmatx);
        if(isempty(statmatx))
            jumpcost = 0;
        else
            jumpcost = ( (mean(statmatx)+ 3*std(statmatx)) );
            smalJC = mean(statmatx);
        end
    end
    
    if(noOfSamplesInTestSample >= noOfSamplesInRefSample)
        pathCost(1,1) = Dist(1,1);
        for ji = 2:1:(noOfSamplesInTestSample)
            pathCost(1,ji) = Dist(1,ji) ;
            
            if( (pathCost(1,ji) >= (pathCost(1,ji-1) + Dist(1,ji)))  )
                pathCost(1,ji) =   pathCost(1,ji-1) + Dist(1,ji) ;
                pathTargetRw(1,ji) = 1;
                pathTargetCol(1,ji) = ji-1;
            end
        end
        for i = 2:1:noOfSamplesInRefSample
            
            stopMotherRight = min((i-1+(elasticity)),noOfSamplesInTestSample);
            stopMotherLeft = max(((i-1)-(elasticity)),1);
            
            for k = stopMotherLeft:1:stopMotherRight
                stopj = min(((k+1+elasticity)-(max(0,(k-(i-1))))),noOfSamplesInTestSample); % change the abs 
                stopLeft  = max(k,1);
                for j = (stopLeft):1:stopj
                    if (((j-(k+1))) <= 0)
                        if ((j-(k+1)) == -1) % for vertical link
                            costJump = smalJC;
                        else
                            costJump = 0;
                        end
                    else
                        costJump = jumpcost*(weight*(abs(j-(k+1))));
                    end
                    if ((pathCost(i,j)) >  (( (pathCost(i-1,k) + ((Dist(i,j)))) + costJump ) ))
                        pathCost(i,j) = (( (pathCost(i-1,k) + ((Dist(i,j)))) + costJump )) ;
                        
                        pathTargetRw(i,j) = i-1;
                        pathTargetCol(i,j) = k;
                    end
                    takeMax = max(1,j-1);
                    if( (pathCost(i,j) > (pathCost(i,takeMax) + smalJC + Dist(i,j)))  )
                        pathCost(i,j) =   pathCost(i,takeMax) + smalJC + Dist(i,j) ;
                        pathTargetRw(i,j) = i;
                        pathTargetCol(i,j) = j-1;
                    end
                end
            end
        end
        
        lastElasticity = max((noOfSamplesInRefSample),1);
        [minVal,~] = min(pathCost(noOfSamplesInRefSample,lastElasticity:noOfSamplesInTestSample));
        tempArr = pathCost(noOfSamplesInRefSample,lastElasticity:noOfSamplesInTestSample);
        ind = find(tempArr == minVal);
        indices = max(ind);
        
        mincol = (lastElasticity-1)+indices;
        minrow = noOfSamplesInRefSample;
        
        Wrapped =[];
        
        while (minrow>=1 && mincol>=1)
            Wrapped = cat(1,[minrow,mincol],Wrapped);
            if(minrow == 1)&&(mincol == 1)
%                 break;
            end
            mincolTemp = pathTargetCol(minrow,mincol);
            minrow = pathTargetRw(minrow,mincol);
            mincol = mincolTemp;
        end
        indxrow = Wrapped(:,1);
        indxcol = Wrapped(:,2);
        
        distSum = pathCost(indxrow(end,1),indxcol(end,1));
        distSum = normalizationTech_0(indxrow,indxcol,distSum,noOfSamplesInRefSample,noOfSamplesInTestSample,jumpcost);      
        return
    end
end
end


function myMinArr = calcJumpCost(distMat, elasticity)
myMinArr = zeros(1,1);
elasticity = 2;%round(elasticity/2);
for iiii = 1:1:size(distMat,1)
    tempMinArr =  sort(distMat(iiii,:));
    takeSome = tempMinArr(1,1:elasticity);
    myMinArr(iiii,1) = mean(takeSome);
end
end

function distSum = normalizationTech_0(~,indxcol,distSum,~,~,~)
distSum = ((distSum))/(size(indxcol,1));
end
