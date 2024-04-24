% %%%%%%%%%%%%%%%%%% 子業娼業制割 %%%%%%%%%%%%
function outputt=dithering(a)
 b(:,:,1)=[0 0;0 0];
 b(:,:,2)=[0 1;0 0];
 b(:,:,3)=[0 1;1 0];
 b(:,:,4)=[1 0;1 1];
 a=round(a*1000)/1000;
 temp1=floor(a);
 temp2=(a-temp1)+0.0001;
 temp3=ceil(temp2*4);
 temp4=zeros(size(a)*2);
  for i=1:length(a(:,1))
      for j=1:length(a(1,:))
        temp4((i*2-1):(i*2),(j*2-1):(j*2))=temp1(i,j)+b(:,:,temp3(i,j));
      end 
  end
 outputt=temp4;
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

