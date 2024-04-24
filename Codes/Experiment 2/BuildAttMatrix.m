%% attMatrix.index:  number of the trials in main experiment
%% attMatrix.att:  shape change in every second, in each block first lasts 20s,others 4s
%% attMatrix.att 1:  whether change in this second：0 don‘t，1 change
%% attMatrix.att 2:  change direction：1 clockwise/up/bigger，-1 anticlockwise/down/smaller 
%% attMatrix.att 3:  change grating：1 left grating, 2 right grating
%% attMatrix.att 4:  quest staircase of attended stimuli: 1 up staircase, 2 down staircase
%% attMatrix.att 5:  change amplitude
%% attMatrix.att 6:  whether response right
%% attMatrix.att 7:  unrelated feature1: whether change in this second：0 don‘t，1 change
%% attMatrix.att 8:  unrelated feature2: whether change in this second：0 don‘t，1 change

function attMatrix=BuildAttMatrix
ppp=BuildMatrix;
trials=length(ppp);
conditions=length(unique(ppp(:,2)));
trialPerBlock=trials/(2*conditions);
% changeProb=[0.25 0.75];
whetherChange=[0 1];  
changeDire=[-1 1];
changeGrat=[1 2];
totalChangeNum1=0;
totalChangeNum2=0;
stairCaseNum=[1;2];
for i=1:trials
    if mod(i,trialPerBlock)==1
        attnums=20;
    else
        attnums=4;
    end
    for n=1:attnums
        if n>1%first second don't change
            attMatrix(i).index=i;
            attMatrix(i).att(n,1)=randsample(whetherChange,1);
            attMatrix(i).att(n,7)=randsample(whetherChange,1);
            attMatrix(i).att(n,8)=randsample(whetherChange,1);
            if n>1
                if attMatrix(i).att(n-1,1)==1
                    attMatrix(i).att(n,1)=0;
                end
                if attMatrix(i).att(n-1,7)==1
                    attMatrix(i).att(n,7)=0;
                end
                if attMatrix(i).att(n-1,8)==1
                    attMatrix(i).att(n,8)=0;
                end
            end
            attMatrix(i).att(n,2)=randsample(changeDire,1);
            attMatrix(i).att(n,3)=randsample(changeGrat,1);

            if attMatrix(i).att(n,1)==1
                if i<=2*trialPerBlock
                    totalChangeNum1=totalChangeNum1+1;
                else 
                    totalChangeNum2=totalChangeNum2+1;
                end
            end
        end
    end
end

%Find when change happens in different conditions
%condition1
temp1=[];
tempNum1=round(totalChangeNum1/(length(stairCaseNum)));
for i=1:tempNum1
    temp1=[temp1;stairCaseNum];
end
r=randperm(size(temp1,1) );   
temp1B=temp1(r, :);
%condition2
temp2=[];
tempNum2=round(totalChangeNum2/(length(stairCaseNum)));
for i=1:tempNum2
    temp2=[temp2;stairCaseNum];
end
r=randperm(size(temp2,1) );   
temp2B=temp2(r, :);
stairCaseMatrix=[temp1B;temp2B];

%put stairCaseMatrix into attMatrix
s=0;
firstAttNum=[];
for  i=1:trials
    if mod(i,trialPerBlock)==1
        attnums=20;
    else
        attnums=4;
    end
    for n=1:attnums
        if attMatrix(i).att(n,1)==1
            s=s+1;
            attMatrix(i).att(n,4)=stairCaseMatrix(s);
        end
    end
    %find the first att event in every condition
    if mod(i,trialPerBlock*2)==1
        for n=1:20
            if attMatrix(i).att(n,1)==1
                firstAttNum=[firstAttNum,n];
                break
            end
        end
    end

end
            

totalChangeNum1
totalChangeNum2
end



















