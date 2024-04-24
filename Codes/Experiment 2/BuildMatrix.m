%%%%%% BuildMatrix.m
%%%%%% �����Աȶ���ӦSS��������
%%%%%% column 1��trial number
%%%%%% column 2����Ӧ��դ��С��10Ϊ���դֱ����4ΪС��դֱ��
%%%%%% column 3����2IFC��һ�����֣�1�ڵ�һ��������֣�2�ڵڶ���������֣���ȫ�����
%%%%%% column 4��staircase��ţ�1Ϊ�½����У�2Ϊ��������
%%%%%% column 5��stimuli ���Թ�դcontrastֵ
%%%%%% column 6�������жϽ����1�ڵ�һ��������֣�2�ڵڶ����������
%%%%%% column 7���ж��Ƿ���ȷ����ȷ1������0
%%%%%% column 8���Ƿ�Ϊ�յ㣬1�����ǹյ㣬0�����ǹյ�


function Matrix=BuildMatrix

TrialPerCondition=120;  %��ͬ�����Դε�����
reverseTrials=0;
quartertrials=0.25*TrialPerCondition;

GratingSize=[2 5.5];  %��Ӧ��դ��С
StaircaseNum=[1 2];  %staircase���
Question2IFC=[1 2]; %which interval has target

blockcondition=[GratingSize(1);GratingSize(2)];

%%%%����ʵ�鲻ͬ����������
[x0,x1]=ndgrid(GratingSize,StaircaseNum);%���ɲ�����Ͼ��󣬴���ÿһ�еĲ�����Χ��ϳɾ���
CombinePara=[x0(:),x1(:)];%ʵ��trial���е�������ϣ�ÿһ������Ϊ1�У���8��

col=8;%��������һ��7�С�temp��8��7�������д̼����������һ��
temp=zeros(length(CombinePara(:,1)),col);%length��ȡ����CombinePara��Ԫ�ظ�����������
%temp(:,2)=CombinePara(:,1);%grating size  
temp(:,4)=CombinePara(:,2);  %��������

paramatrix_temp1=[];
for j=1:quartertrials
    paramatrix_temp1=[paramatrix_temp1;temp];%   ;���������   ,�Ǻ������
end    

paramatrix_temp2=zeros(reverseTrials,col);
length1=length(paramatrix_temp1(:,1));
length2=length(paramatrix_temp2(:,1));
length0=length1+length2;

paramatrix=[];
blockindex=randperm(length(blockcondition(:,1)));%������block������˳���������
for kk=1:2
      index=randperm(length(paramatrix_temp1(:,1)));%%%���Ҵ̼��Ĵ���
      for i=1+(kk-1)*length0:(kk-1)*length0+length1 %for i=1+(kk-1)*16:kk*16
            a=i-(kk-1)*length0;
            paramatrix(i,:)=paramatrix_temp1(index(a),:);%ʹ�ý��������������ÿ64��trial�������
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

