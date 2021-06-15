% Generate percentages for P1:9 for training
sd = 1.5; %adjust this as needed
Ps = 2*normpdf([0.5:8.5],0,sd) %use this output
sumPs = sum(2*normpdf([0.5:8.5],0,sd)) %check it adds up

%Estimate trials for 30 minutes
oneside = normpdf([0.5:8.5],0,sd);
%Both sides evenly
(oneside+fliplr(oneside)).*300
%Subset in even for oneside
oneside.*300