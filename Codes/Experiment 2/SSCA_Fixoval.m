% by zhoayuwan
% SScontrast with attention on grating orientation during adaptation, subject has to response to any orientation change as soon as possible
% Using Quest staircase to measure contrast threshold after adaptation
% Using Eyelink to record subjects' eyemovement and pupil size

function SSCA_Fixoval(date,subjectname,start,over)
% Screen('Preference', 'SkipSyncTests', 1);     %正式实验的时候去掉
if nargin<1                                                         %如果没有输入参数的话
    date='111101';
    subjectname='zhaoyuwan';
    start=1;
    over=1;

end

subjectid=strcat(date,subjectname);
startnumber=num2str(start);
% Screen('Preference', 'SkipSyncTests', 1);
load newclut;
load oldclut;
Screen('LoadNormalizedGammaTable',0,newclut);%write CLUT, screen normalization  clut，颜色查找表
WaitSecs(2);

if start==1
    paramatrix=BuildMatrix;  %建立参数矩阵
    attMatrix=BuildAttMatrix;
else
    load([subjectid  '_f_paramatrix']);  
end
HideCursor;%隐藏鼠标光标
blocks=length(unique(paramatrix(:,2)));
TrialperBlock=length(paramatrix)/blocks;
 %%%%%%%%%%%%%%%%%%%%%%打开屏幕%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 background=128;%background,光栅背景亮度
 [window,windowRect]=Screen('OpenWindow', 0, background,[0 0 1600 1200]);
 frame_rate=Screen('FrameRate',window);  %帧率，单位赫兹，每秒有多少帧
 framedur=1000/frame_rate;%单位ms，每一帧有多长时间
 
dummymode=1;
%% -------------------------------------------------
EyelinkInit(0); % Initialize EyeLink connection
status = Eyelink('IsConnected');
if status == 1 % If EyeLink is connected
    dummymode = 0; 
end
%%-------------------------------------------------
%-------------------------------------------------
%%eyelink settings
if ~dummymode
    edfFile = [date '.edf'];
    Eyelink('OpenFile', edfFile);%创建记录眼动数据的文件open file to record data to
    el = EyelinkInitDefaults(window);    
    el.calibrationtargetsize = 3;% Outer target size as percentage of the screen
    el.calibrationtargetwidth = 0.7;% Inner target size as percentage of the screen
    el.backgroundcolour = [1 1 1]*background;% RGB black
    el.calibrationtargetcolour = [255 255 255];% RGB black
    el.msgfontcolour = [255 255 255];% RGB black
    EyelinkUpdateDefaults(el);
    % This command is crucial to map the gaze positions from the tracker to screen pixel positions to determine fixation
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, windowRect(3)-1,  windowRect(4)-1);
    Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, windowRect(3)-1, windowRect(4)-1);    
    % set calibration type.
    Eyelink('Command', 'calibration_type = HV9'); % horizontal-vertical 9-points
    EyelinkDoTrackerSetup(el);%眼动仪校正
%     preambleText = sprintf('RECORDED BY Psychtoolbox demo %s session name: %s', mfilename, edfFile);
%     Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);%记录文件修改历史
end
%%-------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%基本刺激参数%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%周期
 DistanceToScreen=66;                                  %屏幕距离 厘米
 WidthOfScreen=40;                                        %屏幕宽度 厘米
 XResolution=windowRect(3);                        %屏幕水平分辨率。全屏幕的最右的像素点坐标，即水平方向有多少个像素点，即屏幕的水平分辨率。类似地，windowRect(4)屏幕的垂直分辨率
 PixelSize=atan((WidthOfScreen/XResolution)/DistanceToScreen)*180/pi;%每个像素的视角度数
 PixelsPerDeg=1/PixelSize;                              %每个像素的视角度数的倒数，即在一度视角中有多少个像素
 fx=windowRect(3)/2; %屏幕中心横坐标   
 fy=windowRect(4)/2; %屏幕中心纵坐标
 fixRange=2;
 Eccentricity(1,:)=round([5 0]/PixelSize);%偏离幕中心的距离（x y）方向
 %%%%注视点 
 FixLum=0;  %注视点的灰度
 FixationRect=CenterRect([0 0 2 2],windowRect);%注视点居中
 FixationRect1=OffsetRect(FixationRect,Eccentricity(1,1),Eccentricity(1,2));%%适应光栅
 FixationRect2=OffsetRect(FixationRect,-Eccentricity(1,1),-Eccentricity(1,2));
FixDuration=300;
%  if eye==1
%       FixDuration=0;
%  else
%      FixDuration=500;
%  end
 %%duration of blanks
 AdaptTestBlank=400;%测试光栅和适应光栅之间0.3秒黑屏间隔 blank screen between adaptation and tests 0.3s
 InterTestBlank=400;%2IFC中间的间隔时间 blank screen between 2IFC
 durationOfTest=200;% 测试光栅呈现时间200ms,其中75ms出现，75ms消失 duration of test 200ms,with 75ms on and 75ms off

 %%按键设置
 keycode=zeros(1,256);
 KbName('UnifyKeyNames') 
     qkey=KbName('q');
     leftkey=KbName('leftarrow');
     rightkey=KbName('rightarrow'); 
     spacekey=KbName('space');

%% -----------------------2IFC测试任务阶梯法参数（QUEST）----------------
%%Quest初始参数设置
StrengthRange(1,:)=[0.01,0.15];%测试光栅对比度
tGuess(1:2)=[StrengthRange(1,1) StrengthRange(1,2)];% 先验阈限估计值
tGuessSd=2;% 先验阈限估计值的标准差(注：快稳定的时候，标准差是多少)
pThreshold=0.75;% 阈限临界值处反应正确的概率
beta=1.5;% 控制心理测量函数的斜率
delta=0.01;% 被试盲目反应比例
gamma=0.5;% 随机概率

% %%orientation任务：光栅朝向变化刺激强度（偏转的角度）
% %%Quest初始参数设置
% maxstrength2(2)=20;
% minstrength2(2)=0.5;
% tGuess2(1:2)=[minstrength2(2) maxstrength2(2)];% 先验阈限估计值
% tGuessSd2=2;% 先验阈限估计值的标准差(注：快稳定的时候，标准差是多少)
% pThreshold2=0.83;% 阈限临界值处反应正确的概率
% beta2=1.5;% 控制心理测量函数的斜率
% delta2=0.01;% 被试盲目反应比例
% gamma2=0.5;% 随机概率

%%fixoval任务：椭圆的胖瘦变化（椭圆length变化的大小）
%%Quest初始参数设置
maxstrength2(2)=0.2;
minstrength2(2)=0.001;
tGuess2(1:2)=[minstrength2(2) maxstrength2(2)];% 先验阈限估计值
tGuessSd2=2;% 先验阈限估计值的标准差(注：快稳定的时候，标准差是多少)
pThreshold2=0.83;% 阈限临界值处反应正确的概率
beta2=1.5;% 控制心理测量函数的斜率
delta2=0.01;% 被试盲目反应比例
gamma2=0.5;% 随机概率

%find the first att event in every condition
firstAttNum=[];
for i=1:length(paramatrix)
    if mod(i,TrialperBlock)==1
        for n=1:20
            if attMatrix(i).att(n,1)==1
                firstAttNum=[firstAttNum,n];
                break
            end
        end
    end
end


 %% 提示音参数 Error sound parameter
ff=700;%%声音频率
ff1=1000;
duration=100;%%声音持续的时间，单位毫秒
samplefrq=10000;%%声音信号的采样频率
errorsound_data=MakeBeep(ff,duration/1000,samplefrq);
errorsound_data1=MakeBeep(ff1,200/1000,samplefrq);
FramesofSound=round(duration/framedur);
 %% 椭圆参数 Oval parameter
oval_length=round(0.6*PixelsPerDeg);%椭圆长度
oval_width=round(0.3*PixelsPerDeg);%椭圆宽度
oval_red=[200 50 50];%红色椭圆
oval_black=[50 50 50];%灰色椭圆

ovalrect1=[0 0 oval_width oval_length];%椭圆1大小
ovalrect1=CenterRect(ovalrect1,windowRect);%%位置

ovalrect2=[0 0 oval_length oval_width];%椭圆2大小
ovalrect2=CenterRect(ovalrect2,windowRect);%%位置

oval_dure=1000;%椭圆呈现时长
img_dure=100;%每幅画面呈现时长
oval_frames=oval_dure/1000*frame_rate;%椭圆帧数
img_frames=img_dure/1000*frame_rate;%每幅画面帧数

ovalColor=oval_black;
ovalRect0=ovalrect1;
oval_width_change=1;
upordown=[1 -1];

global probe kd
     
%% 光栅相关参数 Grating parameter
 SpatialPeriod=round(PixelsPerDeg/4);            %PixelsPerDeg取整，除以几就是几度视角
 SpatialPeriodDithering=round(SpatialPeriod/2);
 CounterPhaseFrequency=2;%反向闪烁2hz
 FramesofOneFlickerPeriod=round(frame_rate/CounterPhaseFrequency);% 一个周期500ms，即50帧（一个周期是指光栅波峰/波谷交替后又恢复的整段时间） 
 FramesofResPerPeriod=FramesofOneFlickerPeriod-FramesofSound;
 %% %%%%%%%%%%%%%%%%%%%按任意键开始实验%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   while KbCheck; end
   Screen('DrawText',window,'Fixation oval task',FixationRect(1)-160,FixationRect(2)-50,[0 0 0]);%提示开始  
   Screen('DrawText',window,'Press any key to begin',FixationRect(1)-160,FixationRect(2),[0 0 0]);%提示开始  
   Screen('Flip', window);   
   KbWait;
 
    for mm=start:over
        
        isbreak=0;%先赋值为0，如果isbreak值为1，则跳出此循环
        
        q_up=QuestCreate(log10(tGuess(1)),tGuessSd,pThreshold,beta,delta,gamma);% 上升序列
        q_down=QuestCreate(log10(tGuess(2)),tGuessSd,pThreshold,beta,delta,gamma);% 下降序列

        qAtt_up=QuestCreate(log10(tGuess2(1)),tGuessSd2,pThreshold2,beta2,delta2,gamma2);% 上升序列
        qAtt_down=QuestCreate(log10(tGuess2(2)),tGuessSd2,pThreshold2,beta2,delta2,gamma2);% 下降序列

        
        for i=1+(mm-1)*TrialperBlock:TrialperBlock*mm
%         for i=63:80
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

             %%设置Quest初始值%% 
            StaircaseNum=paramatrix(i,4);
            if StaircaseNum==1
              contrast=max([min([10^(QuestMean(q_up)) 0.99]) 0.0001]);
            elseif StaircaseNum==2
              contrast=max([min([10^(QuestMean(q_down)) 0.99]) 0.0001]);
            end
            paramatrix(i,5)=contrast;
            condition1=round(paramatrix(i,2));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%生成光栅%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%光栅大小，直径，像素
            adaptGratingSize=round(paramatrix(i,2)*PixelsPerDeg);%%适应光栅大小,这里需要加载paramatrix
            testGratingSize=round(1.75*PixelsPerDeg);%%测试光栅大小
            testGratingSizeDithering=round(0.5*testGratingSize);  
            adaptGratingRect=CenterRect([0 0 adaptGratingSize adaptGratingSize],windowRect);
            testGratingRect=CenterRect([0 0 testGratingSize testGratingSize],windowRect);

            %光栅呈现位置        
            adaptGratingRect1(1,:)=OffsetRect(adaptGratingRect,Eccentricity(1,1),Eccentricity(1,2));%%适应光栅
            adaptGratingRect1(2,:)=OffsetRect(adaptGratingRect,-Eccentricity(1,1),-Eccentricity(1,2));

            testGratingRect1(1,:)=OffsetRect(testGratingRect,Eccentricity(1,1),Eccentricity(1,2));%%测试光栅
            testGratingRect1(2,:)=OffsetRect(testGratingRect,-Eccentricity(1,1),-Eccentricity(1,2));

            %%%%适应光栅         
            adaptGratingPhase1=randsample(0:pi/6:2*pi,1);%每个trial的相位随机     start:step:end         
            %中央
            i_size1=adaptGratingSize; %i_size,中央光栅大小，像素
            i_Speriod1=SpatialPeriod;%i_Speriod,中央光栅的空间周期，像素
            i_flicker1=CounterPhaseFrequency;%i_tf,中央光栅闪动频率，赫兹
            i_phase1=adaptGratingPhase1;%i_phase,中央光栅的相位
            i_contrast1=0.75;%i_contrast,中央光栅的对比度
            i_orientation1=0;%i_orientation,中央光栅的朝向            
            i_blur1=round(0.15/PixelSize);%i_blur,中央光栅边缘的正弦模糊像素数
            i_TRU1=0;%i_TRU,中央光栅以余弦方式平缓出现的帧数
            i_TRD1=0;%i_TRD,中央光栅以余弦方式平缓消失的帧数
            
            %外周
            o_size1=adaptGratingSize; %o_size,外周光栅大小，像素
            o_Speriod1=SpatialPeriod;%o_Speriod,外周光栅的空间周期，像素
            o_flicker1=0;%o_tf,外周光栅闪动频率，赫兹
            o_phase1=0;%o_phase,外周光栅的相位
            o_contrast1=0;%o_contrast,外周光栅的对比度
            o_orientation1=0;%o_orientation,外周光栅的朝向
            o_blur1=10;%o_blur,外周光栅边缘的正弦模糊像素数
            o_TRU1=0;%o_TRU,外周光栅以余弦方式平缓出现的帧数
            o_TRD1=0;%o_TRD,外周光栅以余弦方式平缓消失的帧数

            %%%%测试光栅
%             DurationofTest=durationOfTest;
            testGratingFrames=round((durationOfTest/2)/framedur);
            gratingPhase2=randsample(0:pi/6:2*pi,1);%每个trial的相位随机
            %中央
            i_size2=testGratingSizeDithering;%i_size,中央光栅大小，直径，像素
            i_Speriod2=SpatialPeriodDithering;%i_Speriod,中央光栅的空间周期，像素
            i_flicker2=0;%i_tf,中央光栅闪动频率，赫兹
            i_phase2=gratingPhase2;%i_phase,中央光栅的相位
            i_contrast2=paramatrix(i,5);%要调用paramatrix 
            i_orientation2=0;%i_orientation,中央光栅的朝向
            i_blur2=round(0.33/PixelsPerDeg);%i_blur,中央光栅边缘的正弦模糊像素数
            i_TRU2=testGratingFrames/2;%i_TRU,中央光栅以余弦方式平缓出现的帧数,仅生成一半帧数，后面倒序呈现
            i_TRD2=0;%i_TRD,中央光栅以余弦方式平缓消失的帧数
            gap=0;%中央外周之间的间隔
            %外周
            o_size2=testGratingSizeDithering;%o_size,光栅大小，直径，像素
            o_Speriod2=SpatialPeriodDithering;%o_Speriod,光栅的空间周期，像素
            o_flicker2=0;%o_tf,外周光栅闪动频率，赫兹
            o_phase2=gratingPhase2;%o_phase,外周光栅的相位
            o_contrast2=0;%o_contrast,外周光栅的对比度
            o_orientation2=0;%o_orientation,外周光栅的朝向
            o_blur2=round(0.15/PixelsPerDeg);%o_blur,外周光栅边缘的正弦模糊像素数
            o_TRU2=testGratingFrames/2;%o_TRU,外周光栅以余弦方式平缓出现的帧数，仅生成一半帧数，后面倒序呈现
            o_TRD2=0;%o_TRD,外周光栅以余弦方式平缓消失的帧数
            %%%%%%%%生成光栅矩阵

            %适应光栅
            adaptGratingmatrix=flickergratingmatrix(frame_rate,background,FramesofOneFlickerPeriod,o_size1,o_Speriod1,o_flicker1,o_phase1,o_contrast1,o_orientation1,o_blur1,o_TRU1,o_TRD1,...
                                          gap,i_size1,i_Speriod1,i_flicker1,i_phase1,i_contrast1,i_orientation1,i_blur1,i_TRU1,i_TRD1);  
            adaptGratingmatrix=round(adaptGratingmatrix);       
            
            %测试光栅
            testGratingmatrix=flickergratingmatrix(frame_rate,background,testGratingFrames,o_size2,o_Speriod2,o_flicker2,o_phase2,o_contrast2,o_orientation2,o_blur2,o_TRU2,o_TRD2,...
                                          gap,i_size2,i_Speriod2,i_flicker2,i_phase2,i_contrast2,i_orientation2,i_blur2,i_TRU2,i_TRD2);  
            testGratingmatrix=dithering(testGratingmatrix); 
            testGratingmatrix=round(testGratingmatrix);   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%呈现刺激%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% -------------------------------------------------
            if ~dummymode
                Eyelink('Command', 'set_idle_mode');
                WaitSecs(0.05);
                Eyelink('startrecording');  
                Eyelink('Message', 'TrialNumber %d; Size %d;StaircaseNumber %d',i,condition1,StaircaseNum);
                error=Eyelink('CheckRecording');
                if(error~=0)
                    break;
                end
                eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                if eye_used == el.BINOCULAR % if both eyes are tracked
                    eye_used = el.LEFT_EYE; % use left eye
                end
            end
          %%-------------------------------------------------
          %%-------------------------------------------------
            if ~dummymode
                Eyelink('Message', 'Adapt %d',condition1);
            end
          %% -------------------------------------------------
            %呈现适应光栅
            if mod(i,TrialperBlock/2)==1%每种条件下第一个trial呈现时间为20s，之后呈现时间为4s
                %注视点
                for ii=1:round(FixDuration/framedur)  
                    Screen('FillOval',window,FixLum,FixationRect);%fixation 
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);
                    Screen('Flip', window);
                end
                for aa=1:40
                    if mod(aa,4)==1
                        Eyemoved=0;
                    end
                    n=floor((aa-1)/2)+1;%att target second
                    if mod(aa,2)==1%start of a att target
                        probe=1;
                        kd=0;%是否按键
                    end
                    leftRotation=0;
                    rightRotation=0;
                    
                    if attMatrix(i).att(n,1)==1
                        attStairCaseNum=attMatrix(i).att(n,4);
                        if attStairCaseNum==1
                          attChangeAmp=max([min([10^(QuestMean(qAtt_up)) 0.2]) 0.0001]);
                        elseif attStairCaseNum==2
                          attChangeAmp=max([min([10^(QuestMean(qAtt_down)) 0.2]) 0.0001]);
                        end
                        attMatrix(i).att(n,5)=attChangeAmp;
                        attStrength=round(attChangeAmp*attMatrix(i).att(n,2)*PixelsPerDeg);
                        oval_width1=oval_width+attStrength;
                        ovalRect1=[0 0 oval_width1 oval_length];
                        ovalRect=CenterRect(ovalRect1,windowRect);   
                    else 
                        ovalRect=CenterRect(ovalRect0,windowRect);
                    end
                    if attMatrix(i).att(n,7)==1
                        rotation=9*(upordown(:,randperm(2,1)));
                        if attMatrix(i).att(n,3)==1
                            leftRotation=rotation;
                            rightRotation=0;
                        elseif attMatrix(i).att(n,3)==2
                            leftRotation=0;
                            rightRotation=rotation;
                        end
                    end
                     if attMatrix(i).att(n,8)==1
                         changeCon=0.05*(upordown(:,randperm(2,1)));
                        if attMatrix(i).att(n,3)==1
                            leftContrast=0.75+changeCon;
                        elseif attMatrix(i).att(n,3)==2
                            rightContrast=0.75+changeCon;
                        end
                    end
                    for pp=1:round(FramesofOneFlickerPeriod)
                        if mod(aa,2)==0 && pp==FramesofResPerPeriod%check response, 100ms before every second
                            if (attMatrix(i).att(n,1) == 1 && kd==1) || (attMatrix(i).att(n,1) == 0 && kd==0)    
                                attMatrix(i).att(n,6)=1;%反应正确
                            else
                                attMatrix(i).att(n,6)=0;%反应错误
                                if kd==0 && probe
                                    attMatrix(i).att(n,6)=0;
                                    errorsound_player=audioplayer(errorsound_data,samplefrq);%%声音的handle
                                    play(errorsound_player);
                                end
                            end
                            if attMatrix(i).att(n,1) == 1          
                               if (mod(i,TrialperBlock)~=1)||((mod(i,TrialperBlock)==1)&&(n>firstAttNum(mm)))
                                    attJudge=attMatrix(i).att(n,6);%%正确为1，错误为0
                                    if attStairCaseNum==1 % 上升序列
                                      qAtt_up=QuestUpdate(qAtt_up,log10(attChangeAmp),attJudge);
                                    elseif attStairCaseNum==2 % 下降序列
                                      qAtt_down=QuestUpdate(qAtt_down,log10(attChangeAmp),attJudge);
                                    end
                               end
                            end
                            attMatrix(i).att(n,9)=kd;
                        end
                        AdaptGratingTexture(pp)=Screen('MakeTexture',window,adaptGratingmatrix(:,:,pp)); 
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(1,:),leftRotation);
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(2,:),rightRotation);
                        Screen('FillOval',window,ovalColor,ovalRect); 
                        Screen('FillOval',window,FixLum,FixationRect1);
                        Screen('FillOval',window,FixLum,FixationRect2);
                        Screen('Flip', window);  
                        Screen('Close',AdaptGratingTexture(pp));
                        %%收集按键反应
                        if keycode(qkey)~=1 && probe %&& mod(aa,2)>0
                            [keyisdown,~,keycode]=KbCheck;
                            if keycode(spacekey)==1
                                kd=1;
                                if attMatrix(i).att(n,1) == 0
                                    attMatrix(i).att(n,6)=-1;
%                                     Beeper(1000);
                                errorsound_player=audioplayer(errorsound_data,samplefrq);%%声音的handle
                                play(errorsound_player);
                                end
                                probe=0;  
                            end
                        elseif keycode(qkey)==1
                            %%-------------------------------------------------
                            if ~dummymode
                                Eyelink('Command', 'set_idle_mode');
                                WaitSecs(0.5);
                                Eyelink('CloseFile');
                            
                                try
                                    fprintf('Receiving data file ''%s''\n', edfFile);
                                    status = Eyelink('ReceiveFile',edfFile,[subjectid startnumber '_f_Eyelink.edf']);
                                    if status > 0
                                        fprintf('ReceiveFile status %d\n', status);
                                    end
                            %             if 2==exist(edfFile, 'file')
                            %                 fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
                            %             end
                                catch %#ok<*CTCH>
                                    fprintf('Problem receiving data file ''%s''\n', edfFile );
                                end
                            end
                            sca;
                            return
                        end
                        %% -------------------------------------------------
                        if ~dummymode %
                            if pp==1
                                conuntEye=0;
                                xgaze=[];
                                ygaze=[];
                            end
                %             error=Eyelink('CheckRecording');
                %             if(error~=0)
                %                 break;
                %             end

                            if Eyelink( 'NewFloatSampleAvailable') > 0
                                % get the sample in the form of an event structure
                                evt = Eyelink( 'NewestFloatSample');
                                if eye_used ~= -1 % do we know which eye to use yet?
                                    % if we do, get current gaze position from sample

                                    eyex = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                    eyey = evt.gy(eye_used+1);
                                    eyeppd = evt.rx;


                                    % do we have valid data and is the pupil visible?
                                    if eyex~=el.MISSING_DATA && eyey~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                                        conuntEye=conuntEye+1;
                                        xgaze(conuntEye)=eyex;
                                        ygaze(conuntEye)=eyey;
                                    end
                                end
                            end
                            if pp==round(FramesofOneFlickerPeriod)
                                mx=mean(xgaze);
                                my=mean(ygaze);
                                if (mx-fx)^2>=(fixRange*eyeppd)^2
                                    Eyemoved=Eyemoved+1;
                                    attMatrix(i).att(n,10)=mx;
                                    attMatrix(i).att(n,11)=my;
                                end
                            end
%                         evtype = Eyelink('GetNextDataType');            
%                         if evtype ~= el.STARTBLINK && evtype ~= el.ENDBLINK
%                             evt = Eyelink('GetFloatData', evtype); % access the sample structure             
%                             % Save current gaze x y sample fields in variables. See EyeLink Programmers Guide manual > Data Structures > FEVENT
%                             if eye_used ~= -1
%                                 x_gaze = evt.gx(eye_used+1); % +1 as we are accessing an array
%                                 y_gaze = evt.gy(eye_used+1);                      
%                                 if (x_gaze-fx)^2>=fixRange^2
%                                     Eyemoved=Eyemoved+1;
%                                     attMatrix(i).att(n,10)=x_gaze;
%                                     attMatrix(i).att(n,11)=y_gaze;
%                                 end
%                             end
%                         end   
                        end                    
                    %% -------------------------------------------------
                        
                    end
                   if mod(aa,4)==0
                        if Eyemoved>=1
                            [x,fs] = audioread('keepfixed.wav');
                            sound(x,fs);
                            paramatrix(i,8)=Eyemoved;
                        end
                    end 
                end
            else
                for aa=1:8
                    n=floor((aa-1)/2)+1;%att target second
                    if mod(aa,2)==1%start of a att target
                        probe=1;
                        kd=0;%是否按键
                    end
                    if mod(aa,4)==1
                        Eyemoved=0;
                    end
                    leftRotation=0;
                    rightRotation=0;
                    if attMatrix(i).att(n,1)==1
                        attStairCaseNum=attMatrix(i).att(n,4);
                        if attStairCaseNum==1
                          attChangeAmp=max([min([10^(QuestMean(qAtt_up)) 0.2]) 0.0001]);
                        elseif attStairCaseNum==2
                          attChangeAmp=max([min([10^(QuestMean(qAtt_down)) 0.2]) 0.0001]);
                        end
                        attMatrix(i).att(n,5)=attChangeAmp;
                        attStrength=round(attChangeAmp*attMatrix(i).att(n,2)*PixelsPerDeg);
                        oval_width1=oval_width+attStrength;
                        ovalRect1=[0 0 oval_width1 oval_length];
                        ovalRect=CenterRect(ovalRect1,windowRect);   
                    else 
                        ovalRect=CenterRect(ovalRect0,windowRect);
                    end
                    if attMatrix(i).att(n,7)==1
                        rotation=9*(upordown(:,randperm(2,1)));
                        if attMatrix(i).att(n,3)==1
                            leftRotation=rotation;
                            rightRotation=0;
                        elseif attMatrix(i).att(n,3)==2
                            leftRotation=0;
                            rightRotation=rotation;
                        end
                    end
                    for pp=1:round(FramesofOneFlickerPeriod)
                        if mod(aa,2)==0 && pp==FramesofResPerPeriod%check response, 100ms before every second
                            if (attMatrix(i).att(n,1) == 1 && kd==1) || (attMatrix(i).att(n,1) == 0 && kd==0)    
                                attMatrix(i).att(n,6)=1;%反应正确
                            else
                                attMatrix(i).att(n,6)=0;%反应错误
                                if kd==0 && probe
                                    attMatrix(i).att(n,6)=0;
                                    errorsound_player=audioplayer(errorsound_data,samplefrq);%%声音的handle
                                    play(errorsound_player);
                                end
                            end
                            attMatrix(i).att(n,9)=kd;
                            if attMatrix(i).att(n,1) == 1
                                if (mod(i,TrialperBlock)~=1)||((mod(i,TrialperBlock)==1)&&(n>firstAttNum(mm)))
                                   attJudge=attMatrix(i).att(n,6);%%正确为1，错误为0
                                    if attStairCaseNum==1 % 上升序列
                                      qAtt_up=QuestUpdate(qAtt_up,log10(attChangeAmp),attJudge);
                                    elseif attStairCaseNum==2 % 下降序列
                                      qAtt_down=QuestUpdate(qAtt_down,log10(attChangeAmp),attJudge);
                                    end
                                end
                            end
                        end
                        AdaptGratingTexture(pp)=Screen('MakeTexture',window,adaptGratingmatrix(:,:,pp)); 
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(1,:),leftRotation);
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(2,:),rightRotation);
                        Screen('FillOval',window,ovalColor,ovalRect); 
                        Screen('FillOval',window,FixLum,FixationRect1);
                        Screen('FillOval',window,FixLum,FixationRect2);    
                        Screen('Flip', window);  
                        Screen('Close',AdaptGratingTexture(pp));
                        %%收集按键反应
                        if keycode(qkey)~=1 && probe %&& mod(pp,10)==0%check keyboard every 100 second
                            [keyisdown,~,keycode]=KbCheck;
                            if keycode(spacekey)==1
                                kd=1;
                                if attMatrix(i).att(n,1) == 0
                                    attMatrix(i).att(n,6)=-1;
                                errorsound_player=audioplayer(errorsound_data,samplefrq);%%声音的handle
                                play(errorsound_player);
                                end
                                probe=0;  
                            end
                        elseif keycode(qkey)==1
                            %%-------------------------------------------------
                            if ~dummymode
                                    Eyelink('Command', 'set_idle_mode');
                                    WaitSecs(0.5);
                                    Eyelink('CloseFile');
                                
                                    try
                                        fprintf('Receiving data file ''%s''\n', edfFile);
                                        status = Eyelink('ReceiveFile',edfFile,[subjectid startnumber '_f_Eyelink.edf']);
                                        if status > 0
                                            fprintf('ReceiveFile status %d\n', status);
                                        end
                                %             if 2==exist(edfFile, 'file')
                                %                 fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
                                %             end
                                    catch %#ok<*CTCH>
                                        fprintf('Problem receiving data file ''%s''\n', edfFile );
                                    end
                            end
                            sca;
                            return
                        end
                         %% -------------------------------------------------
                         if ~dummymode   
                            if pp==1
                                conuntEye=0;
                                xgaze=[];
                                ygaze=[];
                            end
                %             error=Eyelink('CheckRecording');
                %             if(error~=0)
                %                 break;
                %             end

                            if Eyelink( 'NewFloatSampleAvailable') > 0
                                % get the sample in the form of an event structure
                                evt = Eyelink( 'NewestFloatSample');
                                if eye_used ~= -1 % do we know which eye to use yet?
                                    % if we do, get current gaze position from sample

                                    eyex = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                    eyey = evt.gy(eye_used+1);
                                    eyeppd = evt.rx;


                                    % do we have valid data and is the pupil visible?
                                    if eyex~=el.MISSING_DATA && eyey~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                                        conuntEye=conuntEye+1;
                                        xgaze(conuntEye)=eyex;
                                        ygaze(conuntEye)=eyey;
                                    end
                                end
                            end
                            if pp==round(FramesofOneFlickerPeriod)
                                mx=mean(xgaze);
                                my=mean(ygaze);
                                if (mx-fx)^2>=(fixRange*eyeppd)^2
                                    Eyemoved=Eyemoved+1;
                                    attMatrix(i).att(n,10)=mx;
                                    attMatrix(i).att(n,11)=my;
                                end
                            end  
                         end

                    %% -------------------------------------------------
                    end 
                    if mod(aa,4)==0
                        if Eyemoved>=1
                            [x,fs] = audioread('keepfixed.wav');
                            sound(x,fs);
                            paramatrix(i,8)=Eyemoved;
                        end
                    end
                end
            end 

            
            %呈现0.3秒的黑屏、注视点
            for ii=1:round(AdaptTestBlank/framedur)  
                 Screen('FillOval',window,FixLum,FixationRect);%fixation 
                 Screen('FillOval',window,FixLum,FixationRect1);
                 Screen('FillOval',window,FixLum,FixationRect2);    
                 Screen('Flip', window);
            end
          %% -------------------------------------------------
            if ~dummymode
                Eyelink('Message', 'Test %d',i);
            end
          %%  -------------------------------------------------
            %呈现测试光栅
            if paramatrix(i,3)==1
                %第一个间隔呈现光栅
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%声音的handle
                play(errorsound_player);
                for pp=1:round(testGratingFrames)
                    testGratingTexture(pp)=Screen('MakeTexture',window,testGratingmatrix(:,:,pp)); 
                    Screen('DrawTexture',window,testGratingTexture(pp),[],testGratingRect1(1,:));
                    Screen('DrawTexture',window,testGratingTexture(pp),[],testGratingRect1(2,:));
                    Screen('FillOval',window,FixLum,FixationRect);      
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);  
                    %% -------------------------------------------------
                     if ~dummymode   
                        if pp==1
                            testEyemoved=0;
                            conuntTestEye=0;
                            xtestgaze=[];
                            ytestgaze=[];
                        end
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            evt = Eyelink( 'NewestFloatSample');
                            if eye_used ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                eyex = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                eyey = evt.gy(eye_used+1);
                                eyeppd = evt.rx;
                                % do we have valid data and is the pupil visible?
                                if eyex~=el.MISSING_DATA && eyey~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                                    conuntTestEye=conuntTestEye+1;
                                    xtestgaze(conuntTestEye)=eyex;
                                    ytestgaze(conuntTestEye)=eyey;
                                end
                            end
                        end
                        if pp==round(testGratingFrames)
                            mx=mean(xtestgaze);
                            my=mean(ytestgaze);
                            if (mx-fx)^2>=(fixRange*eyeppd)^2
                                testEyemoved=testEyemoved+1;
                                attMatrix(i).att(n,10)=mx;
                                attMatrix(i).att(n,11)=my;
                            end
                        end  
                     end
                    %% -------------------------------------------------
                    Screen('Close',testGratingTexture(pp));
                end
                %中间间隔400ms
                for ii=1:round(InterTestBlank/framedur) 
                    Screen('FillOval',window,FixLum,FixationRect);     
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);
                end
                %第二个间隔不呈现光栅
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%声音的handle
                play(errorsound_player);
                for pp=1:round(testGratingFrames)
                    Screen('FillOval',window,FixLum,FixationRect);     
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);
                    %% -------------------------------------------------
                     if ~dummymode   
                        if pp==1
                            conuntTestEye=0;
                            xtestgaze=[];
                            ytestgaze=[];
                        end
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            evt = Eyelink( 'NewestFloatSample');
                            if eye_used ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                eyex = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                eyey = evt.gy(eye_used+1);
                                eyeppd = evt.rx;
                                % do we have valid data and is the pupil visible?
                                if eyex~=el.MISSING_DATA && eyey~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                                    conuntTestEye=conuntTestEye+1;
                                    xtestgaze(conuntTestEye)=eyex;
                                    ytestgaze(conuntTestEye)=eyey;
                                end
                            end
                        end
                        if pp==round(testGratingFrames)
                            mx=mean(xtestgaze);
                            my=mean(ytestgaze);
                            if (mx-fx)^2>=(fixRange*eyeppd)^2
                                testEyemoved=testEyemoved+1;
                                attMatrix(i).att(n,10)=mx;
                                attMatrix(i).att(n,11)=my;
                            end
                        end  
                     end
                    %% -------------------------------------------------
                end
            elseif paramatrix(i,3)==2
                %第一个间隔不呈现光栅
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%声音的handle
                play(errorsound_player);
                for pp=1:round(testGratingFrames)
                    Screen('FillOval',window,FixLum,FixationRect);
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);
                    %% -------------------------------------------------
                     if ~dummymode   
                        if pp==1
                            testEyemoved=0;
                            conuntTestEye=0;
                            xtestgaze=[];
                            ytestgaze=[];
                        end
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            evt = Eyelink( 'NewestFloatSample');
                            if eye_used ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                eyex = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                eyey = evt.gy(eye_used+1);
                                eyeppd = evt.rx;
                                % do we have valid data and is the pupil visible?
                                if eyex~=el.MISSING_DATA && eyey~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                                    conuntTestEye=conuntTestEye+1;
                                    xtestgaze(conuntTestEye)=eyex;
                                    ytestgaze(conuntTestEye)=eyey;
                                end
                            end
                        end
                        if pp==round(testGratingFrames)
                            mx=mean(xtestgaze);
                            my=mean(ytestgaze);
                            if (mx-fx)^2>=(fixRange*eyeppd)^2
                                testEyemoved=testEyemoved+1;
                                attMatrix(i).att(n,10)=mx;
                                attMatrix(i).att(n,11)=my;
                            end
                        end  
                     end
                    %% -------------------------------------------------
                end
                %中间间隔400ms
                for ii=1:round(InterTestBlank/framedur) 
                    Screen('FillOval',window,FixLum,FixationRect);
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);
                end
                %第二个间隔呈现光栅
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%声音的handle
                play(errorsound_player);
                for pp=1:round(testGratingFrames)
                    testGratingTexture(pp)=Screen('MakeTexture',window,testGratingmatrix(:,:,pp)); 
                    Screen('DrawTexture',window,testGratingTexture(pp),[],testGratingRect1(1,:));
                    Screen('DrawTexture',window,testGratingTexture(pp),[],testGratingRect1(2,:));
                    Screen('FillOval',window,FixLum,FixationRect);
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);  
                    %% -------------------------------------------------
                     if ~dummymode   
                        if pp==1
                            conuntTestEye=0;
                            xtestgaze=[];
                            ytestgaze=[];
                        end
                        if Eyelink( 'NewFloatSampleAvailable') > 0
                            % get the sample in the form of an event structure
                            evt = Eyelink( 'NewestFloatSample');
                            if eye_used ~= -1 % do we know which eye to use yet?
                                % if we do, get current gaze position from sample
                                eyex = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                                eyey = evt.gy(eye_used+1);
                                eyeppd = evt.rx;
                                % do we have valid data and is the pupil visible?
                                if eyex~=el.MISSING_DATA && eyey~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                                    conuntTestEye=conuntTestEye+1;
                                    xtestgaze(conuntTestEye)=eyex;
                                    ytestgaze(conuntTestEye)=eyey;
                                end
                            end
                        end
                        if pp==round(testGratingFrames)
                            mx=mean(xtestgaze);
                            my=mean(ytestgaze);
                            if (mx-fx)^2>=(fixRange*eyeppd)^2
                                testEyemoved=testEyemoved+1;
                                attMatrix(i).att(n,10)=mx;
                                attMatrix(i).att(n,11)=my;
                            end
                        end  
                     end
                    %% -------------------------------------------------
                    Screen('Close',testGratingTexture(pp));
                end
            end

            %% -------------------------------------------------
            if ~dummymode   
                if testEyemoved>=1
                    [x,fs] = audioread('keepfixed.wav');
                    sound(x,fs);
                    paramatrix(i,8)=Eyemoved+testEyemoved;
                end
            end
            %% -------------------------------------------------

                     
           %呈现问题，判断刺激方向
            while true
                Screen('DrawText',window,'?',FixationRect(1),FixationRect(2),[0 0 0]);
                Screen('Flip', window);      
                [keyDown,keyTime,KeyCode]=KbCheck;%读取键盘
                if KeyCode(rightkey) || KeyCode(leftkey) || KeyCode(qkey)
                    break;
                end
            end
                        
            if KeyCode(qkey)==1%跳出小循环
                isbreak=1;%后面用于跳出大循环
                break                       
            elseif KeyCode(rightkey)==1
                paramatrix(i,6)=2;
            elseif KeyCode(leftkey)==1
                paramatrix(i,6)=1;
            end
            
            %%判断被试反应是否正确
            if paramatrix(i,3)==paramatrix(i,6)
                paramatrix(i,7)=1;
            else
                paramatrix(i,7)=0;
                if paramatrix(i,4)~=0
                    Beeper(500,0.8,0.1);   %反馈音效
                end
            end     
            
            %%%%更新Quest程序
            if mod(i,TrialperBlock)~=1 %%排除第一个试次          
               judgement=paramatrix(i,7);%%正确为1，错误为0
                if StaircaseNum==1 % 上升序列
                  q_up=QuestUpdate(q_up,log10(contrast),judgement);
                elseif StaircaseNum==2 % 下降序列
                  q_down=QuestUpdate(q_down,log10(contrast),judgement);
                end
            end
                        
             %%next trial
             Screen('FillOval',window,FixLum,FixationRect);%fixation 
             Screen('FillOval',window,FixLum,FixationRect1);
             Screen('FillOval',window,FixLum,FixationRect2);    
             Screen('Flip', window);
             WaitSecs(0.01);
             
        %% -------------------------------------------------
        if ~dummymode
            Eyelink('Message', 'Trial %d End',i);
            Eyelink('Stoprecording');
        end
        %% -------------------------------------------------
        save([subjectid '_f_paramatrix'],'paramatrix','attMatrix');  
        
        if TrialperBlock>=80
            if i==TrialperBlock/2+(mm-1)*TrialperBlock
                Screen('DrawText',window,'Please have a rest...',FixationRect(1)-160,FixationRect(2),[0 0 0]);
                Screen('Flip', window); 
                WaitSecs(20);
                while KbCheck; end 
                Screen('DrawText',window,'Press any key to continue',FixationRect(1)-160,FixationRect(2),[0 0 0]);
                Screen('Flip', window);   
                KbWait;
            end
        end 

        
        end  %到此为止一个condition结束
                
        clear responsematrix strengthmatrix  responsematrix1 strengthmatrix1  responsematrix2 strengthmatrix2 %清除数据，每一个block重新开始阶梯法测α值
                           
        %%%%%%存储数据%%%%%%%
        
    
        if isbreak==1%跳出大循环
            break
        end
        
        if mm==over
            Screen('DrawText',window,'You have finished the experiment',FixationRect(1)-160,FixationRect(2),[0 0 0]);
            Screen('Flip', window); 
            WaitSecs(1);
            break; 
        end
        
        %%每一个condition结束之后休息2分钟 
        Screen('DrawText',window,'Please have a rest...',FixationRect(1)-160,FixationRect(2),[0 0 0]);
        Screen('Flip', window); 
        WaitSecs(300);

        %%按键继续实验
        while KbCheck; end 
        Screen('DrawText',window,'Press any key to continue',FixationRect(1)-160,FixationRect(2),[0 0 0]);
        Screen('Flip', window);   
        KbWait;

    
    end%四个condition的实验全部结束
    
    
    Screen('CloseAll');
    Screen('LoadNormalizedGammaTable',0,oldclut);%write CLUT, screen normalization
%%-------------------------------------------------
if ~dummymode
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');

    try
        fprintf('Receiving data file ''%s''\n', edfFile);
        status = Eyelink('ReceiveFile',edfFile,[subjectid startnumber '_f_Eyelink.edf']);
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
%             if 2==exist(edfFile, 'file')
%                 fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
%             end
    catch %#ok<*CTCH>
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
end
