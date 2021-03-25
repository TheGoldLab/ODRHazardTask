% Generate percentages for P1:9 for training
sd = 1.5; %adjust this as needed
Ps = 2*normpdf([0.5:8.5],0,sd) %use this output
sumPs = sum(2*normpdf([0.5:8.5],0,sd)) %check it adds up