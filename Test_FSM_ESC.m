% This algorithm is about checking FSM and ESC algorithm with a set of numbers
clear ;
% close all;
clc;

query1 = [1,2,9,16,9,25,37]';
target1 = [1,2,9,16,9,25,37]';

query2 = [1,2,8,4,6,8,5,6,7]';
target2 = [1,2,3,8,6,8,9,12,8]';

query3 = [1,2,8,8,8]';
target3 = [1,2,9,95,79,26,39,31]';

query4 = [1,2,8,8]';
target4 = [1,2,95,95,95,8,8]';


[pathCostFSM,indxcolFSM,indxrowFSM,distSumFSM,jumpcostFSM] = FSM_Algo(query2, target2,1);

[pathCostESC,indxcolESC,indxrowESC,distSumESC,jumpcostESC] = ESC_Algo(query4, target4,1);
 

