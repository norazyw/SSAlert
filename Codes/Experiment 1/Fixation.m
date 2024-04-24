%%%%%% High aletness condition
%%%%%% Center Go no-Go RT task
%%%%%% 1up2down to measure contrast threshold


function Fixation(date,subjectname,start,over)
% Screen('Preference', 'SkipSyncTests', 1);     %��ʽʵ���ʱ��ȥ��
if nargin<1                                                         %���û����������Ļ�
    date='050800';   
    subjectname='test';
    start=1;
    over=1;
end

subjectid=strcat(date,subjectname);
startnumber=num2str(start);

load newclut;
load oldclut;
Screen('LoadNormalizedGammaTable',0,newclut);%write CLUT, screen normalization  clut����ɫ���ұ�
WaitSecs(2);

if start==1
    paramatrix=BuildMatrix;  %������������
    ovalMatrix=BuildOvalMatrix;
else
    load(['Data/' subjectid  '_f_paramatrix']);  
end
HideCursor;%���������
TrakingEyes=eye;
TrialperBlock=0.5*length(paramatrix);%ʵ���Ϊ4��block
 %%%%%%%%%%%%%%%%%%%%%%����Ļ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 background=128;%background,��դ��������
 [window,windowRect]=Screen('OpenWindow', 0, background,[0 0 1600 1200]);
 frame_rate=Screen('FrameRate',window);  %֡�ʣ���λ���ȣ�ÿ���ж���֡
 framedur=1000/frame_rate;%��λms��ÿһ֡�ж೤ʱ��
 

%% -------------------------------------------------
dummymode=1;
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
    Eyelink('OpenFile', edfFile);%������¼�۶����ݵ��ļ�open file to record data to
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
    EyelinkDoTrackerSetup(el);%�۶���У��
%     preambleText = sprintf('RECORDED BY Psychtoolbox demo %s session name: %s', mfilename, edfFile);
%     Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);%��¼�ļ��޸���ʷ
end
%%-------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%�����̼�����%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%����
 DistanceToScreen=78;                                  %��Ļ���� ����
 WidthOfScreen=36.5;                                        %��Ļ��� ����
 XResolution=windowRect(3);                        %��Ļˮƽ�ֱ��ʡ�ȫ��Ļ�����ҵ����ص����꣬��ˮƽ�����ж��ٸ����ص㣬����Ļ��ˮƽ�ֱ��ʡ����Ƶأ�windowRect(4)��Ļ�Ĵ�ֱ�ֱ���
 PixelSize=atan((WidthOfScreen/XResolution)/DistanceToScreen)*180/pi;%ÿ�����ص��ӽǶ���
 PixelsPerDeg=1/PixelSize;                              %ÿ�����ص��ӽǶ����ĵ���������һ���ӽ����ж��ٸ�����
 fx=windowRect(3)/2; %��Ļ���ĺ�����   
 fy=windowRect(4)/2; %��Ļ����������
 fixRange=2;
 xCenter=windowRect(3)/2; %��Ļ���ĺ�����   
 yCenter=windowRect(4)/2; %��Ļ����������
 Eccentricity(1,:)=round([5 0]/PixelSize);%ƫ��Ļ���ĵľ��루x y������
 %%%%ע�ӵ� 
 FixLum=0;  %ע�ӵ�ĻҶ�
 FixationRect=CenterRect([0 0 2 2],windowRect);%ע�ӵ����
 FixationRect1=OffsetRect(FixationRect,Eccentricity(1,1),Eccentricity(1,2));%%��Ӧ��դ
 FixationRect2=OffsetRect(FixationRect,-Eccentricity(1,1),-Eccentricity(1,2));
 if eye==1
      FixDuration=300;
 else
     FixDuration=500;
 end
 %%duration of blanks
 AdaptTestBlank=400;%���Թ�դ����Ӧ��դ֮��0.3�������� blank screen between adaptation and tests 0.3s
 InterTestBlank=300;%2IFC�м�ļ��ʱ�� blank screen between 2IFC
 durationOfTest=200;% ���Թ�դ����ʱ��200ms,����75ms���֣�75ms��ʧ duration of test 200ms,with 75ms on and 75ms off

 %%��������
 keycode=zeros(1,256);
 KbName('UnifyKeyNames') 
     qkey=KbName('q');
     leftkey=KbName('leftarrow');
     rightkey=KbName('rightarrow'); 
     spacekey=KbName('space');

%% -----------------------2IFC����������ݷ�������QUEST��----------------
%%Quest��ʼ��������
StrengthRange(1,:)=[0.01,0.15];%���Թ�դ�Աȶ�
tGuess(1:2)=[StrengthRange(1,1) StrengthRange(1,2)];% �������޹���ֵ
tGuessSd=2;% �������޹���ֵ�ı�׼��(ע�����ȶ���ʱ�򣬱�׼���Ƕ���)
pThreshold=0.75;% �����ٽ�ֵ����Ӧ��ȷ�ĸ���
beta=1.5;% �����������������б��
delta=0.01;% ����äĿ��Ӧ����
gamma=0.5;% �������



 %% ��ʾ������ Error sound parameter
ff=700;%%����Ƶ��
ff1=1000;
ff0=500;
duration=100;%%����������ʱ�䣬��λ����
samplefrq=10000;%%�����źŵĲ���Ƶ��
errorsound_data0=MakeBeep(ff0,duration/1000,samplefrq);
errorsound_data=MakeBeep(ff,duration/1000,samplefrq);
errorsound_data1=MakeBeep(ff1,200/1000,samplefrq);
FramesofSound=round(duration/framedur);
 %% ��Բ���� Oval parameter
oval_length=round(0.4*PixelsPerDeg);%��Բ����
oval_width=round(0.3*PixelsPerDeg);%��Բ���
oval_red=[200 50 50];%��ɫ��Բ
oval_black=[50 50 50];%��ɫ��Բ

ovalrect0_1=[0 0 oval_width oval_length];%��Բ1��С
ovalrect1=CenterRect(ovalrect0_1,windowRect);%%λ��

ovalrect0_2=[0 0 oval_length oval_width];%��Բ2��С
ovalrect2=CenterRect(ovalrect0_2,windowRect);%%λ��

oval_dure=1000;%��Բ����ʱ��
img_dure=100;%ÿ���������ʱ��
oval_frames=oval_dure/1000*frame_rate;%��Բ֡��
img_frames=img_dure/1000*frame_rate;%ÿ������֡��

global probe kd
     
%% ��դ��ز��� Grating parameter
 SpatialPeriod=round(PixelsPerDeg/4);            %PixelsPerDegȡ�������Լ����Ǽ����ӽ�
 SpatialPeriodDithering=round(SpatialPeriod/2);
 CounterPhaseFrequency=2;%������˸2hz
 FramesofOneFlickerPeriod=round(frame_rate/CounterPhaseFrequency);% һ������500ms����50֡��һ��������ָ��դ����/���Ƚ�����ָֻ�������ʱ�䣩 
 FramesofResPerPeriod=FramesofOneFlickerPeriod-FramesofSound;
 %% %%%%%%%%%%%%%%%%%%%���������ʼʵ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   while KbCheck; end
   Screen('DrawText',window,'Only fixation',FixationRect(1)-100,FixationRect(2)-50,[175 0 0]);%��ʾ��ʼ  
   Screen('DrawText',window,'Press any key to begin',FixationRect(1)-160,FixationRect(2),[0 0 0]);%��ʾ��ʼ  
   Screen('Flip', window);   
   KbWait;
 
    for mm=start:over
        
        isbreak=0;%�ȸ�ֵΪ0�����isbreakֵΪ1����������ѭ��
        
        %%�������ݵĳ�ʼֵ
        q_up1=QuestCreate(log10(tGuess(1)),tGuessSd,pThreshold,beta,delta,gamma);% ��������
        q_down1=QuestCreate(log10(tGuess(2)),tGuessSd,pThreshold,beta,delta,gamma);% �½�����

        for i=1+(mm-1)*TrialperBlock:TrialperBlock*mm
%         for i=63:80
            Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

             %%����Quest��ʼֵ%% 
            StaircaseNum=paramatrix(i,4);
            if StaircaseNum==1
              contrast=max([min([10^(QuestMean(q_up1)) 0.99]) 0.0001]);
            elseif StaircaseNum==2
              contrast=max([min([10^(QuestMean(q_down1)) 0.99]) 0.0001]);
            end
            paramatrix(i,5)=contrast;
            condition1=round(paramatrix(i,2));
%              %��ʼֵ
%              maxstrength(1)=0.1;%�̼����ֶԱȶȱ��ʵ����ֵ
%              minstrength(1)=0.02;%�̼����ֶԱȶȱ��ʵ���Сֵ
%              initialstrength(1,:)=[minstrength(1) maxstrength(1)];
%             %������ݷ�����ֵ
%             if paramatrix(i,4)~=0
%                 StaircaseNum=paramatrix(i,4);  %2,�������У�1,�½�����
%                 if StaircaseNum==1 %�������½������ǲ�ͬ�Ľ�������
%                     strengthmatrix=strengthmatrix1;
%                     responsematrix=responsematrix1;
%                 elseif StaircaseNum==2
%                     strengthmatrix=strengthmatrix2;
%                     responsematrix=responsematrix2;
%                 end
%                 contrast=updownstaircase(staircasetype,strengthmatrix,responsematrix,minstrength,maxstrength,initialstrength(StaircaseNum));
%                 paramatrix(i,5)=contrast;
%                 if StaircaseNum==1
%                     strengthmatrix1=[strengthmatrix;paramatrix(i,5)];
%                 elseif StaircaseNum==2
%                     strengthmatrix2=[strengthmatrix;paramatrix(i,5)];
%                 end
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%���ɹ�դ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%��դ��С��ֱ��������
            adaptGratingSize=round(paramatrix(i,2)*PixelsPerDeg);%%��Ӧ��դ��С,������Ҫ����paramatrix
            testGratingSize=round(1.75*PixelsPerDeg);%%���Թ�դ��С
            testGratingSizeDithering=round(0.5*testGratingSize);  
            adaptGratingRect=CenterRect([0 0 adaptGratingSize adaptGratingSize],windowRect);
            testGratingRect=CenterRect([0 0 testGratingSize testGratingSize],windowRect);

            %��դ����λ��        
            adaptGratingRect1(1,:)=OffsetRect(adaptGratingRect,Eccentricity(1,1),Eccentricity(1,2));%%��Ӧ��դ
            adaptGratingRect1(2,:)=OffsetRect(adaptGratingRect,-Eccentricity(1,1),-Eccentricity(1,2));

            testGratingRect1(1,:)=OffsetRect(testGratingRect,Eccentricity(1,1),Eccentricity(1,2));%%���Թ�դ
            testGratingRect1(2,:)=OffsetRect(testGratingRect,-Eccentricity(1,1),-Eccentricity(1,2));

            %%%%��Ӧ��դ
            if mod(i,TrialperBlock)==1
                DurationofAdapt=20*1000;
            else 
                DurationofAdapt=4000;
            end
            
            adaptGratingPhase1=randsample(0:pi/6:2*pi,1);%ÿ��trial����λ���     start:step:end         
            %����
            i_size1=adaptGratingSize; %i_size,�����դ��С������
            i_Speriod1=SpatialPeriod;%i_Speriod,�����դ�Ŀռ����ڣ�����
            i_flicker1=CounterPhaseFrequency;%i_tf,�����դ����Ƶ�ʣ�����
            i_phase1=adaptGratingPhase1;%i_phase,�����դ����λ
            i_contrast1=0.75;%i_contrast,�����դ�ĶԱȶ�
            if paramatrix(i,4)==0
                i_contrast1=0.1;
            end
            i_orientation1=0;%i_orientation,�����դ�ĳ���            
            i_blur1=round(0.15/PixelSize);%i_blur,�����դ��Ե������ģ��������
            i_TRU1=0;%i_TRU,�����դ�����ҷ�ʽƽ�����ֵ�֡��
            i_TRD1=0;%i_TRD,�����դ�����ҷ�ʽƽ����ʧ��֡��
            %����
            o_size1=adaptGratingSize; %o_size,���ܹ�դ��С������
            o_Speriod1=SpatialPeriod;%o_Speriod,���ܹ�դ�Ŀռ����ڣ�����
            o_flicker1=0;%o_tf,���ܹ�դ����Ƶ�ʣ�����
            o_phase1=0;%o_phase,���ܹ�դ����λ
            o_contrast1=0;%o_contrast,���ܹ�դ�ĶԱȶ�
            o_orientation1=0;%o_orientation,���ܹ�դ�ĳ���
            o_blur1=10;%o_blur,���ܹ�դ��Ե������ģ��������
            o_TRU1=0;%o_TRU,���ܹ�դ�����ҷ�ʽƽ�����ֵ�֡��
            o_TRD1=0;%o_TRD,���ܹ�դ�����ҷ�ʽƽ����ʧ��֡��

            %%%%���Թ�դ
%             DurationofTest=durationOfTest;
            testGratingFrames=round((durationOfTest/2)/framedur);
            gratingPhase2=randsample(0:pi/6:2*pi,1);%ÿ��trial����λ���
            %����
            i_size2=testGratingSizeDithering;%i_size,�����դ��С��ֱ��������
            i_Speriod2=SpatialPeriodDithering;%i_Speriod,�����դ�Ŀռ����ڣ�����
            i_flicker2=0;%i_tf,�����դ����Ƶ�ʣ�����
            i_phase2=gratingPhase2;%i_phase,�����դ����λ
            i_contrast2=paramatrix(i,5);%Ҫ����paramatrix 
            i_orientation2=0;%i_orientation,�����դ�ĳ���
            i_blur2=round(0.33/PixelsPerDeg);%i_blur,�����դ��Ե������ģ��������
            i_TRU2=testGratingFrames/2;%i_TRU,�����դ�����ҷ�ʽƽ�����ֵ�֡��,������һ��֡�������浹�����
            i_TRD2=0;%i_TRD,�����դ�����ҷ�ʽƽ����ʧ��֡��
            gap=0;%��������֮��ļ��
            %����
            o_size2=testGratingSizeDithering;%o_size,��դ��С��ֱ��������
            o_Speriod2=SpatialPeriodDithering;%o_Speriod,��դ�Ŀռ����ڣ�����
            o_flicker2=0;%o_tf,���ܹ�դ����Ƶ�ʣ�����
            o_phase2=gratingPhase2;%o_phase,���ܹ�դ����λ
            o_contrast2=0;%o_contrast,���ܹ�դ�ĶԱȶ�
            o_orientation2=0;%o_orientation,���ܹ�դ�ĳ���
            o_blur2=round(0.15/PixelsPerDeg);%o_blur,���ܹ�դ��Ե������ģ��������
            o_TRU2=testGratingFrames/2;%o_TRU,���ܹ�դ�����ҷ�ʽƽ�����ֵ�֡����������һ��֡�������浹�����
            o_TRD2=0;%o_TRD,���ܹ�դ�����ҷ�ʽƽ����ʧ��֡��
            
            %%%%%%%%���ɹ�դ����

            %��Ӧ��դ
            adaptGratingmatrix=flickergratingmatrix(frame_rate,background,FramesofOneFlickerPeriod,o_size1,o_Speriod1,o_flicker1,o_phase1,o_contrast1,o_orientation1,o_blur1,o_TRU1,o_TRD1,...
                                          gap,i_size1,i_Speriod1,i_flicker1,i_phase1,i_contrast1,i_orientation1,i_blur1,i_TRU1,i_TRD1);  
            adaptGratingmatrix=round(adaptGratingmatrix);       
            
            %���Թ�դ
            testGratingmatrix=flickergratingmatrix(frame_rate,background,testGratingFrames,o_size2,o_Speriod2,o_flicker2,o_phase2,o_contrast2,o_orientation2,o_blur2,o_TRU2,o_TRD2,...
                                          gap,i_size2,i_Speriod2,i_flicker2,i_phase2,i_contrast2,i_orientation2,i_blur2,i_TRU2,i_TRD2);  
            testGratingmatrix=dithering(testGratingmatrix); 
            testGratingmatrix=round(testGratingmatrix);   
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %%%%%%%%%%%%%%%%%%%���ִ̼�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %������Ӧ��դ
            if mod(i,TrialperBlock/2)==1%ÿ�������µ�һ��trial����ʱ��Ϊ20s��֮�����ʱ��Ϊ4s
                %ע�ӵ�
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
%                 for aa=1:25
                    n=floor((aa-1)/2)+1;%center target number
                    if mod(aa,2)==1%start of a center target
                        probe=1;
                        kd=0;%�Ƿ񰴼�
                    end
                    if ovalMatrix(i).oval(n,1)==1
                        ovalRect=ovalrect1;
                    elseif ovalMatrix(i).oval(n,1)==2
                        ovalRect=ovalrect2;
                    end
                    if ovalMatrix(i).oval(n,2)==1
                        ovalColor=oval_red;
                    elseif ovalMatrix(i).oval(n,2)==2
                        ovalColor=oval_black;
                    end
                    
                    for pp=1:round(FramesofOneFlickerPeriod)
                        if mod(aa,2)==0 && pp==FramesofResPerPeriod% end and start of one oval
%                             if (ovalMatrix(i).oval(n,1) == 1 && kd==1) || (ovalMatrix(i).oval(n,1) ~= 1 && kd==0)    
%                                 ovalMatrix(i).oval(n,4)=1;%��Ӧ��ȷ
%                             else
%                                 ovalMatrix(i).oval(n,4)=-1;%��Ӧ����
%                                 if kd==0 && probe
%                                     ovalMatrix(i).oval(n,4)=0;
%                                     errorsound_player=audioplayer(errorsound_data,samplefrq);%%������handle
%                                     play(errorsound_player);
%                                 end
%                             end
                            ovalMatrix(i).oval(n,5)=kd;
                        end
                        %%start of one oval
                        if mod(aa,2)==1 && pp==1
                            tic;
                        end
                        AdaptGratingTexture(pp)=Screen('MakeTexture',window,adaptGratingmatrix(:,:,pp)); 
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(1,:));
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(2,:));
                        Screen('FillOval',window,ovalColor,ovalRect); 
                        Screen('FillOval',window,FixLum,FixationRect1);
                        Screen('FillOval',window,FixLum,FixationRect2);
                        Screen('Flip', window);  
                        Screen('Close',AdaptGratingTexture(pp));
                        %%�ռ�������Ӧ
                        if keycode(qkey)~=1 && probe %&& mod(aa,2)>0
                            [keyisdown,~,keycode]=KbCheck;
                            if keycode(spacekey)==1
                                kd=1;
                                kdtime=toc;
                                ovalMatrix(i).oval(n,6)=kdtime; 
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
                                    status = Eyelink('ReceiveFile');
                                    filename=['Data/' subjectid startnumber '_f_Eyelink.edf'];
                                    copyfile(edfFile,filename);
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
                            %%-------------------------------------------------
%                             QuitCode=0;
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
%                 for aa=1:25
                for aa=1:8
                    if mod(aa,4)==1
                        Eyemoved=0;
                    end
                    n=floor((aa-1)/2)+1;
                    if mod(aa,2)==1
                        probe=1;
                        kd=0;
                    end
                    if ovalMatrix(i).oval(n,1)==1
                        ovalRect=ovalrect1;
                    elseif ovalMatrix(i).oval(n,1)==2
                        ovalRect=ovalrect2;
                    end
                    if ovalMatrix(i).oval(n,2)==1
                        ovalColor=oval_red;
                    elseif ovalMatrix(i).oval(n,2)==2
                        ovalColor=oval_black;
                    end
                    for pp=1:round(FramesofOneFlickerPeriod)
                        if mod(aa,2)==0 && pp==FramesofResPerPeriod
%                             if (ovalMatrix(i).oval(n,1) == 1 && kd==1) || (ovalMatrix(i).oval(n,1) ~= 1 && kd==0)    
%                                 ovalMatrix(i).oval(n,4)=1;%��Ӧ��ȷ
%                             else
%                                 ovalMatrix(i).oval(n,4)=-1;%��Ӧ����
%                                 if kd==0 && probe
%                                     ovalMatrix(i).oval(n,4)=0;
%                                     errorsound_player=audioplayer(errorsound_data,samplefrq);%%������handle
%                                     play(errorsound_player);
%                                 end
%                             end
                            ovalMatrix(i).oval(n,5)=kd;
                        end
                        %%start of one oval
                        if mod(aa,2)==1 && pp==1
                            tic;
                        end
                        AdaptGratingTexture(pp)=Screen('MakeTexture',window,adaptGratingmatrix(:,:,pp)); 
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(1,:));
                        Screen('DrawTexture',window,AdaptGratingTexture(pp),[],adaptGratingRect1(2,:));
                        Screen('FillOval',window,ovalColor,ovalRect); 
                        Screen('FillOval',window,FixLum,FixationRect1);
                        Screen('FillOval',window,FixLum,FixationRect2);    
                        Screen('Flip', window);  
                        Screen('Close',AdaptGratingTexture(pp));
                        %%�ռ�������Ӧ
                        if keycode(qkey)~=1 && probe %&& mod(pp,10)==0%check keyboard every 100 second
                            [keyisdown,~,keycode]=KbCheck;
                            if keycode(spacekey)==1
                                kd=1;
                                kdtime=toc;
                                ovalMatrix(i).oval(n,6)=kdtime; 
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
                                        status = Eyelink('ReceiveFile');
                                        filename=['Data/' subjectid startnumber '_f_Eyelink.edf'];
                                        copyfile(edfFile,filename);
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
                            %%-------------------------------------------------
%                             QuitCode=0;
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

            
            %����0.3��ĺ�����ע�ӵ�
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

            %���ֲ��Թ�դ
            if paramatrix(i,3)==1
                %��һ��������ֹ�դ
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%������handle
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
                            testEyemoved=0;
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
                %�м���400ms
                for ii=1:round(InterTestBlank/framedur) 
                    Screen('FillOval',window,FixLum,FixationRect);     
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window); 
                end
                %�ڶ�����������ֹ�դ
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%������handle
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
                %��һ����������ֹ�դ
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%������handle
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
                            testEyemoved=0;
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
                %�м���400ms
                for ii=1:round(InterTestBlank/framedur) 
                    Screen('FillOval',window,FixLum,FixationRect);
                    Screen('FillOval',window,FixLum,FixationRect1);
                    Screen('FillOval',window,FixLum,FixationRect2);    
                    Screen('Flip', window);
                end
                %�ڶ���������ֹ�դ
%                 Beeper(1000,0.7,0.2);
                errorsound_player=audioplayer(errorsound_data1,samplefrq);%%������handle
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
            
           %�������⣬�жϴ̼�����
            while true
                Screen('DrawText',window,'?',FixationRect(1),FixationRect(2),[0 0 0]);
                Screen('Flip', window);      
                [keyDown,keyTime,KeyCode]=KbCheck;%��ȡ����
                if KeyCode(rightkey) || KeyCode(leftkey) || KeyCode(qkey)
                    break;
                end
            end
                       
            
            if KeyCode(qkey)==1%����Сѭ��
                isbreak=1;%��������������ѭ��
                break                       
            elseif KeyCode(rightkey)==1
                paramatrix(i,6)=2;
            elseif KeyCode(leftkey)==1
                paramatrix(i,6)=1;
            end
            
            %%�жϱ��Է�Ӧ�Ƿ���ȷ
            if paramatrix(i,3)==paramatrix(i,6)
                paramatrix(i,7)=1;
            else
                paramatrix(i,7)=0;
                if paramatrix(i,4)~=0
                    errorsound_player=audioplayer(errorsound_data0,samplefrq);%%������handle
                    play(errorsound_player);
                end
            end     
            
            %%%%����Quest����
            if mod(i,TrialperBlock)~=1 %%�ų���һ���Դ�          
               judgement=paramatrix(i,7);%%��ȷΪ1������Ϊ0
                if StaircaseNum==1 % ��������
                  q_up1=QuestUpdate(q_up1,log10(contrast),judgement);
                elseif StaircaseNum==2 % �½�����
                  q_down1=QuestUpdate(q_down1,log10(contrast),judgement);
                end
            end
%             responsematrix=[responsematrix;paramatrix(i,7)];
%             if StaircaseNum==1
%                 responsematrix1=responsematrix;%�Ѵ˴�trial���Եķ�Ӧ����ŵ����ݷ��½����У�1���ķ�Ӧ������ȥ
%             elseif StaircaseNum==2
%                 responsematrix2=responsematrix;
%             end
                        
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

        save(['Data/' subjectid '_f_paramatrix'],'paramatrix','ovalMatrix');  
        
        if TrialperBlock>=80
            if i==TrialperBlock/2+(mm-1)*TrialperBlock
                Screen('DrawText',window,'Only fixation',FixationRect(1)-100,FixationRect(2)-50,[175 0 0]);%��ʾ��ʼ  
                Screen('DrawText',window,'Please have a rest...',FixationRect(1)-160,FixationRect(2),[0 0 0]);
                Screen('Flip', window); 
                WaitSecs(20);
                while KbCheck; end 
                Screen('DrawText',window,'Only fixation',FixationRect(1)-100,FixationRect(2)-50,[175 0 0]);%��ʾ��ʼ  
                Screen('DrawText',window,'Press any key to continue',FixationRect(1)-160,FixationRect(2),[0 0 0]);
                Screen('Flip', window);   
                KbWait;
            end
        end 

        
        end  %����Ϊֹһ��condition����
                
        clear responsematrix strengthmatrix  responsematrix1 strengthmatrix1  responsematrix2 strengthmatrix2 %������ݣ�ÿһ��block���¿�ʼ���ݷ����ֵ
                           
        %%%%%%�洢����%%%%%%%
        
    
        if isbreak==1%������ѭ��
            break
        end
        
        if mm==over
            Screen('DrawText',window,'You have finished the experiment',FixationRect(1)-160,350,[0 0 0]);
            Screen('Flip', window); 
            WaitSecs(1);
            break; 
        end
        
        %%ÿһ��condition����֮����Ϣ2���� 
        Screen('DrawText',window,'Only fixation',FixationRect(1)-100,FixationRect(2)-50,[175 0 0]);%��ʾ��ʼ  
        Screen('DrawText',window,'Please have a rest...',FixationRect(1)-160,FixationRect(2),[0 0 0]);
        Screen('Flip', window); 
        WaitSecs(300);

        %%��������ʵ��
        while KbCheck; end 
        Screen('DrawText',window,'Only fixation',FixationRect(1)-100,FixationRect(2)-50,[175 0 0]);%��ʾ��ʼ  
        Screen('DrawText',window,'Press any key to continue',FixationRect(1)-160,FixationRect(2),[0 0 0]);
        Screen('Flip', window);   
        KbWait;

    
    end%�ĸ�condition��ʵ��ȫ������
    
    
    Screen('CloseAll');
    Screen('LoadNormalizedGammaTable',0,oldclut);%write CLUT, screen normalization
    
%%-------------------------------------------------
if ~dummymode
    Eyelink('Command', 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');

    try
        fprintf('Receiving data file ''%s''\n', edfFile);
        status = Eyelink('ReceiveFile');
        filename=['Data/' subjectid startnumber '_f_Eyelink.edf'];
        copyfile(edfFile,filename);
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