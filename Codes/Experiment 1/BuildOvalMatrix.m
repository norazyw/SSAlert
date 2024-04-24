%% ovalMatrix.index:  number of the trials in main experiment
%% ovalMatrix.oval:  shape change in every second, in each block first lasts 30s,others 5s
%% ovalMatrix.oval 1:  oval orientation : 2horizontal, 1vertical
%% ovalMatrix.oval 2:  oval color : 1red, 2black
%% ovalMatrix.oval 3:  orientation*color, if =2, the target

function ovalMatrix=BuildOvalMatrix
trials=length(BuildMatrix);
trialPerBlock=trials/4;%2 conditions, big/small
ovalOriet=[1 2];  
ovalColor=[1 2];
for i=1:trials
    if mod(i,trialPerBlock)==1
        ovalnums=20;
    else
        ovalnums=4;
    end
    for n=1:ovalnums
        ovalMatrix(i).index=i;
        ovalMatrix(i).oval(n,1)=randsample(ovalOriet,1);
        ovalMatrix(i).oval(n,2)=randsample(ovalColor,1);
        if n>1
            while (ovalMatrix(i).oval(n,1)==ovalMatrix(i).oval(n-1,1))&&(ovalMatrix(i).oval(n,2)==ovalMatrix(i).oval(n-1,2))
                ovalMatrix(i).oval(n,1)=randsample(ovalOriet,1);
                ovalMatrix(i).oval(n,2)=randsample(ovalColor,1);
            end
        end
        ovalMatrix(i).oval(n,3)=ovalMatrix(i).oval(n,1).*ovalMatrix(i).oval(n,2);
        
    end
end
end



















