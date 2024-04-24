

%  flickergratingmatrix(frame_rate,background,frames,o_size,o_Speriod,o_tf,o_phase,o_contrast,o_orientation,o_blur,o_TRU,o_TRD,
%                                          gap,i_size,i_Speriod,i_tf,i_phase,i_contrast,i_orientation,i_blur,i_TRU,i_TRD)
%  frame_rate=59;%��Ļˢ����
%  background=128;%background,��դ��������
%  frames=300; %frames,�̼�������ʱ��
%  o_size=200; %o_size,���ܹ�դ��С������
%  o_Speriod=20;%o_Speriod,���ܹ�դ�Ŀռ����ڣ�����
%  o_tf=0;%o_tf,���ܹ�դ����Ƶ�ʣ�����
%  o_phase=0;%o_phase,���ܹ�դ����λ
%  o_contrast=0.5;%o_contrast,���ܹ�դ�ĶԱȶ�
%  o_orientation=0;%o_orientation,���ܹ�դ�ĳ���
%  o_blur=10;%o_blur,���ܹ�դ��Ե������ģ��������
%  o_TRU=100;%o_TRU,���ܹ�դ�����ҷ�ʽƽ�����ֵ�֡��
%  o_TRD=50;%o_TRD,���ܹ�դ�����ҷ�ʽƽ����ʧ��֡��
%  
%  gap=10;%��������֮��ļ��
%  
%  i_size=80; %i_size,�����դ��С������
%  i_Speriod=20;%i_Speriod,�����դ�Ŀռ����ڣ�����
%  i_tf=0;%i_tf,�����դ����Ƶ�ʣ�����
%  i_phase=0;%i_phase,�����դ����λ
%  i_contrast=0.5;%i_contrast,�����դ�ĶԱȶ�
%  i_orientation=0;%i_orientation,�����դ�ĳ���
%  i_blur=10;%i_blur,�����դ��Ե������ģ��������
%  i_TRU=120;%i_TRU,�����դ�����ҷ�ʽƽ�����ֵ�֡��
%  i_TRD=40;%i_TRD,�����դ�����ҷ�ʽƽ����ʧ��֡��

%%2015 10 19
%�ú������������Ĺ�դ����դ�����ܹ�դ�����Ĺ�դ��ɣ�����֮�����һ��gap�������դ�ı�Ե����ģ����
%���ĺ����ܹ�դ�Ĵ�С���ռ�Ƶ�ʡ�����ʱ��Ƶ�ʡ���λ���Աȶȡ����򡢱�Եģ��������������ramp���½�ramp�Ȳ�������ֱ�ӵ��ڡ�
 
function output=flickergratingmatrix(frame_rate,background,frames,o_size,o_Speriod,o_flicker,o_phase,o_contrast,o_orientation,o_blur,o_TRU,o_TRD,...
                                  gap,i_size,i_Speriod,i_flicker,i_phase,i_contrast,i_orientation,i_blur,i_TRU,i_TRD)
   
 
%  o_blur=o_blur+0.0001;
%  i_blur=i_blur+0.0001;

frames=round(frames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y]=meshgrid((1-o_size/2):(o_size/2),(1-o_size/2):(o_size/2));%�̼��Ĵ�С
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if o_blur==0 && i_blur==0
   %��û��ģ����ʱ��������ģ������صļ��㣬���ٳ�������ʱ��
   %%boundary2 ���ܻ�������%%%%%
   r_s=round(o_size/2);%���ܹ�դ�뾶
   r_c=round(i_size/2);%�����դ�뾶
   bound2_1=sign(sign(r_s.^2-x.^2-y.^2)+1);%���ܹ�դ��Բ�α�
   bound2_2=sign(sign(-(r_c+gap).^2+x.^2+y.^2)+1);%�����դ��Բ�α�
   boundary2=bound2_1.*bound2_2;%���ܻ�������
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%boundary3 ����Բ������%%%%
    boundary3=1-sign(sign(-(r_c).^2+x.^2+y.^2)+1);%�����դ��Բ�α�
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
else
    
  %�����դ���ҵ�����%%%%%%%%
   r_c=round(i_size/2);%�����դ�뾶
   rampperiod=2*i_blur;%�����դ����ramp������
   bound1_i_1=(sin(2*pi*(1/rampperiod)*(sqrt(x.^2+y.^2)-r_c)-pi/2)+1)/2;%����Բ�̵����ҵݼ���
   bound1_i_2=sign(sign((r_c).^2-x.^2-y.^2)+1);
   bound1_i_3=sign(sign(-(r_c-rampperiod/2).^2+x.^2+y.^2)+1);
   bound1_i_4=sign(sign((r_c-rampperiod/2).^2-x.^2-y.^2)+1);
   bound1_i=bound1_i_1.*bound1_i_2.*bound1_i_3+bound1_i_4;
   bound1_i=min(bound1_i,1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%

  %���ܹ�դ�ڲ����ҵݼ��ߣ�gap
  rampwidth_s=2*o_blur;%���ܹ�դ�ڱ�����ramp������
  bound1_o_1=(sin(2*pi*(1/rampwidth_s)*(sqrt(x.^2+y.^2)-r_c-gap)-pi/2)+1)/2;
  bound1_o_2=sign(sign(-(r_c+gap).^2+x.^2+y.^2)+1);
  bound1_o_3=sign(sign((r_c+rampwidth_s/2+gap).^2-x.^2-y.^2)+1);
  bound1_o_4=sign(sign(-(r_c+rampwidth_s/2+gap).^2+x.^2+y.^2)+1);
  bound1_o=bound1_o_1.*bound1_o_2.*bound1_o_3+bound1_o_4;
  bound1_o=min(bound1_o,1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %boundary1=bound1_i+bound1_o;

  %%boundary2 ���ܻ�������%%%%%
  r_s=round(o_size/2);%���ܹ�դ�뾶
  bound2_1=sign(sign(r_s.^2-x.^2-y.^2)+1);%���ܹ�դ��Բ�α�
  bound2_2=sign(sign(-r_c.^2+x.^2+y.^2)+1);%�����դ��Բ�α�
  boundary2=bound2_1.*bound2_2.*bound1_o;%���ܻ��������ڻ�����gap+rampģ��
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%boundary3 ����Բ������%%%%
  bound3=sign(sign((r_c).^2-x.^2-y.^2)+1);%�����դ��Բ�α�
  boundary3=bound3.*bound1_i;%����Բ�����������rampģ��
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sin ramp_up window
o_TRU_timepoint=0:o_TRU/(o_TRU-1):o_TRU; %%o_TRU��ʱ���
o_TRU_window=(sin(2*pi*(1/(2*o_TRU))*o_TRU_timepoint-pi/2)+1)/2;

i_TRU_timepoint=0:i_TRU/(i_TRU-1):i_TRU; %%o_TRU��ʱ���
i_TRU_window=(sin(2*pi*(1/(2*i_TRU))*i_TRU_timepoint-pi/2)+1)/2;

%%sin ramp_down window
o_TRD_timepoint=0:o_TRD/(o_TRD-1):o_TRD; %%o_TRU��ʱ���
o_TRD_window=(sin(2*pi*(1/(2*o_TRD))*o_TRD_timepoint+pi/2)+1)/2;

i_TRD_timepoint=0:i_TRD/(i_TRD-1):i_TRD; %%o_TRU��ʱ���
i_TRD_window=(sin(2*pi*(1/(2*i_TRD))*i_TRD_timepoint+pi/2)+1)/2;

%%%%%%��������ܵ�ʱ��ƽ����
o_timewindow=ones(1,frames);
i_timewindow=ones(1,frames);

if o_TRU~=0
    o_timewindow(1:o_TRU)=o_TRU_window;
end

if o_TRD~=0
    o_timewindow((frames-o_TRD+1):frames)=o_TRD_window;
end

if i_TRU~=0
    i_timewindow(1:i_TRU)=i_TRU_window;
end

if i_TRD~=0
    i_timewindow((frames-i_TRD+1):frames)=i_TRD_window;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%���� �ռ�Ƶ�ʵ�%%%%%%%%%%%%%%%
    angle_s=o_orientation;%���ܹ�դ���򣬵�λ ��
    angle_c=i_orientation;%�����դ����, ��λ ��
    ang_s=angle_s*pi/180;%ת�� �� Ϊ ����
    ang_c=angle_c*pi/180;%ת�� �� Ϊ ����
    fre_s=1/o_Speriod;
    fre_c=1/i_Speriod;
    
    aa_s=2*pi*fre_s*cos(ang_s);
    bb_s=2*pi*fre_s*sin(ang_s);
    aa_c=2*pi*fre_c*cos(ang_c);
    bb_c=2*pi*fre_c*sin(ang_c);
    
    con_s=o_contrast;
    con_c=i_contrast;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
for jj=1:frames
    gratingmatrix(:,:,jj)=background*(1+o_timewindow(jj)*cos(2*pi*o_flicker*jj/frame_rate)*con_s*cos(aa_s*x+bb_s*y+o_phase).*boundary2+ ... %���ܹ�դ
                                        i_timewindow(jj)*cos(2*pi*i_flicker*jj/frame_rate)*con_c*cos(aa_c*x+bb_c*y+i_phase).*boundary3);    %�����դ
end

output=gratingmatrix;

% secs2=GetSecs;
% secs2-secs1



% [window,windowRect]=Screen('OpenWindow', 0, background);
% 
% 
% 
% for ii=1:frames
% gratings=Screen('MakeTexture',window,gratingmatrix(:,:,ii));
% Screen('DrawTexture', window,gratings,[],[100 100 300 300]);%��դ
% tt(ii)=Screen('Flip', window);
% end
% 
% 
% Screen('CloseAll')
% 
% 

