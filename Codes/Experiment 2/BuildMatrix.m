%%%%%% BuildMatrix.m
%%%%%% 建立对比度适应SS参数矩阵
%%%%%% column 1：trial number
%%%%%% column 2：适应光栅大小，10为大光栅直径，4为小光栅直径
%%%%%% column 3：在2IFC哪一个出现，1在第一个间隔出现，2在第二个间隔出现（完全随机）
%%%%%% column 4：staircase编号，1为下降序列，2为上升序列
%%%%%% column 5：stimuli 测试光栅contrast值
%%%%%% column 6：被试判断结果，1在第一个间隔出现，2在第二个间隔出现
%%%%%% column 7：判断是否正确，正确1，错误0
%%%%%% column 8：是否为拐点，1代表是拐点，0代表不是拐点


function Matrix=BuildMatrix

TrialPerCondition=120;  %相同条件试次的数量
reverseTrials=0;
quartertrials=0.25*TrialPerCondition;

GratingSize=[2 5.5];  %适应光栅大小
StaircaseNum=[1 2];  %staircase编号
Question2IFC=[1 2]; %which interval has target

blockcondition=[GratingSize(1);GratingSize(2)];

%%%%生成实验不同条件组合情况
[x0,x1]=ndgrid(GratingSize,StaircaseNum);%生成参数组合矩阵，带上每一列的参数范围组合成矩阵
CombinePara=[x0(:),x1(:)];%实验trial所有的条件组合，每一种条件为1行，共8行

col=8;%参数矩阵一共7列。temp是8×7矩阵，所有刺激排列组合在一起
temp=zeros(length(CombinePara(:,1)),col);%length获取向量CombinePara的元素个数，即行数
%temp(:,2)=CombinePara(:,1);%grating size  
temp(:,4)=CombinePara(:,2);  %阶梯序列

paramatrix_temp1=[];
for j=1:quartertrials
    paramatrix_temp1=[paramatrix_temp1;temp];%   ;是竖向组合   ,是横向组合
end    

paramatrix_temp2=zeros(reverseTrials,col);
length1=length(paramatrix_temp1(:,1));
length2=length(paramatrix_temp2(:,1));
length0=length1+length2;

paramatrix=[];
blockindex=randperm(length(blockcondition(:,1)));%把四种block条件的顺序随机排列
for kk=1:2
      index=randperm(length(paramatrix_temp1(:,1)));%%%打乱刺激的次序
      for i=1+(kk-1)*length0:(kk-1)*length0+length1 %for i=1+(kk-1)*16:kk*16
            a=i-(kk-1)*length0;
            paramatrix(i,:)=paramatrix_temp1(index(a),:);%使得阶梯序列随机（在每64个trial内随机）
            paramatrix(i,2)=blockcondition(blockindex(kk),1);
            paramatrix(i,3)=Question2IFC(:,randperm(2,1));
      end    
      for i=kk*length0-length2+1 :kk*length0 %for i=1+(kk-1)*16:kk*16
%             a=i-(kk-1)*length(paramatrix_temp2(:,1));
            paramatrix(i,2)=blockcondition(blockindex(kk),1);
            paramatrix(i,3)=Question2IFC(:,randperm(2,1));
            paramatrix(i,5)=0.1;
      end  
      
      
 end
        
paramatrix(:,1)=[1:length(paramatrix(:,1))];

clear TrialPerCondition
clear GratingSize
clear Orientation
clear StaircaseNum
clear condition1
clear condition2
clear condition3
clear condition4
clear CombinePara
clear col j
clear paramatrix_temp1
clear temp
clear blockcondition0
clear StaircaseNum
clear blockcondition
clear index
clear kk
clear a
clear i
clear x0;clear x1;clear x2;

Matrix=paramatrix;

save('paramatrix');
end

