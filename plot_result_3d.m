function plot_result_3d(dirname,time,pre)
%close all

% 
% 
dirname='../pointnet.pytorch/half_real/results/20170806_231442/';
time='20170806_235646';
pre='ft_signfirst0xyz31000000000boundary20_margin2000002000000000000boundary20cos100000';

% comparing hte following two, finetuned is slightly better (flipped area is less)
% #1:time='20170803_183844';
% dirname='../pointnet.pytorch/half_real/results/20170803_170951/';
% pre='xyz21000000000.0boundary20barrier20002000.0boundary0cos100000';
% ft of #1
% dirname='../pointnet.pytorch/half_real/results/20170806_220515/';
% pre='ft_signfirst0xyz21000000000boundary20_margin2000002000000000boundary20cos100000';
% time='20170806_224344';

% area of faces in the center become larger, but why?
% dirname='../pointnet.pytorch/half_real/results/20170806_231442/';
% time='20170806_235646';
% pre='ft_signfirst0xyz31000000000boundary20_margin2000002000000000000boundary20cos100000';

%% to check: l_area.basepara=2
% dirname='../pointnet.pytorch/half_real/results/20170807_002104/';
% time='20170807_011729';
% pre='ft_signfirst0xyz21000000000boundary20_margin2000002000000000000boundary20cos100000';

%% to check: l_area.basepara=0.2(right bottom)
% dirname='../pointnet.pytorch/half_real/results/20170807_002128/';
% time='20170807_011821';
% pre='ft_signfirst0xyz01000000000boundary20_margin2000002000000000000boundary20cos100000';
% for the above two , don't know why all scaled down? ->10^(-3)
dirname='../pointnet.pytorch/harmonic/results/20170807_233409/';
time='20170808_000749';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
%% harmonic

dirname='../pointnet.pytorch/harmonic/results/20170810_102632/';
time='20170810_110114';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
% +26f (three and human are slightly better, but others are not)
dirname='../pointnet.pytorch/harmonic/results/20170810_234902/';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_024911';
% boundary prediction good, only some cases fail to predict those on head/
% with negative x
dirname='../pointnet.pytorch/harmonic/results/20170811_015144/';
pre='stn1har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_023225';
% stn2 not better either for loss or image
dirname='../pointnet.pytorch/harmonic/results/20170811_100651/';
pre='stn2har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_104139';
%visualize trained samples
dirname='../pointnet.pytorch/harmonic/results/20170811_161624/';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_161625';

dirname='../pointnet.pytorch/harmonic/results/20170811_162026/';
pre='har_signfirst0xyz21000000000boundary20_margin200000.020boundary20cos100000';
time='20170811_162027';

% combine harmonic
%overwelming margin(left up)
dirname='../pointnet.pytorch/harmonic/results/20170812_000854/';
pre='stn1har_signfirst0xyz210000boundary20_margin200000.020000000000boundary20cos100000';
time='20170812_043814';
%overwelming margin+har(right up)
dirname='../pointnet.pytorch/harmonic/results/20170812_001332/';
pre='stn1har_signfirst0xyz210000boundary20_margin200000.020000000000boundary20cos100000';
time='20170812_044106';

%overwelming cos(right down)
dirname='../pointnet.pytorch/harmonic/results/20170811_235955/';
pre='stn1har_signfirst0xyz21000boundary20_margin200000.0200000boundary20cos100000';
time='20170812_040540';
%overoverwelming margin+not weighted average for margin (left down)
dirname='../pointnet.pytorch/harmonic/results/20170812_003706/';
pre='stn1har100000_signfirst0xyz210000boundary20_margin200000.020000000000000N20cos100000';
time='20170812_045312';


%har+margin   better than only har (get mapping with more non-flipping area, but worse for regular 3d mesh like sphere)
dirname='../pointnet.pytorch/harmonic/results/20170812_162434/';
pre='res0stn1har100000_signfirst0N00N20_margin200000.02000N20cos0';   %2e5 is too large
time='20170812_170329';
%har+area    +a bit margin or (less area loss?)   shape is changed,area
%distortion is less
%according to original shape   
dirname='../pointnet.pytorch/harmonic/results/20170812_141649/';
pre='res0stn1har100000_signfirst0N010000N20_margin200000.00N20cos0';
time='20170812_164900';
%har+cos a bit small?  ==>should implement + total area the same
dirname='../pointnet.pytorch/harmonic/results/20170812_142831/';
pre='res0stn1har100000_signfirst0N00N20_margin200000.00N20cos100';
time='20170812_165620';
%first time try to combine area +area_sum: not obvious influence
dirname='../pointnet.pytorch/harmonic/results/20170812_210334/';
time='20170812_211649';
pre='res0stn1har100000_signfirst0N0[10000.0, 1000.0]N20_margin200000.00N20cos0';
%weird: overlapping but not flipped  
dirname='../pointnet.pytorch/harmonic/results/20170812_213607/';
time='20170812_214705';
pre='res0stn1har100000_signfirst0N0[10000.0, 1000.0]N20_margin200000.02000N20cos100';


dirname='../pointnet.pytorch/harmonic/results/20170813_202654/';
pre='res0stn1har100000_signfirst0N0[1, 0]har0.001_margin200000.00N20cos0';
time='20170813_203704';
pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';



%use har margin ft on area+har_24 :
filename='20170813_204753/res0stn1har100000_signfirst0N0[0, 0]N0.001_margin20000.02000har0.001cos0_face_20170813_205801.txt';
%use har margin ft on the above:  it learns to autonomously align the 3d
%half to the same direction, squeeze the highest points to a line
filename='20170813_211526/res0stn1har100000_signfirst0N0[0, 0]N0.001_margin20000.020000har0.01cos0_face_20170813_212540.txt';
filename='20170813_223735/res0stn1har100000_signfirst0N0[10000.0, 0]N0.001_margin200000.02000har0.01cos100_face_20170813_225045.txt';
% double map
%(cos)2nd map performs better on those unseen, but worse on those seen
filename='20170814_013037/res0stn1har100000_signfirst0N0[0.0, 0]N0.001_margin200000.00har0.01cos100_face_20170814_042228.txt';
%(margin) for some squeezed 1st map, the 2nd spread the map again (best
%human#251 till now). (not good for cases backward,margin loss weightperhaps is too large)
filename='20170814_012925/res0stn1har100000_signfirst0N0[0.0, 0]N0.001_margin200000.02000har0.01cos0_face_20170814_041334.txt';
%area:slightly better than single
filename='20170814_013013/res0stn1har100000_signfirst0N0[10000.0, 0]N0.001_margin200000.00har0.01cos0_face_20170814_042146.txt';
%mostly second better than first. a totally flipped: case #171
filename='20170814_130703/res0stn1har100000_signfirst0N0[10000.0, 0]N0.001_margin200000.00har0.01cos0_face_20170814_151506.txt';
%har+area(per-face+sum)+margin+cos  second is better
filename='20170814_142425/res0stn1har100000N0[10000.0, 1000.0]N0.001_margin200000.0200har0.01cos10b0_face_20170814_164323.txt';

filename='20170814_144815/2mapres0stn1har100000N0[10000.0, 1000.0]N0.001_margin200000.00har0.01cos0b10000_face_20170814_170048.txt';
%fairly good
filename='20170815_143556/2mapres0stn1har100000N0[10000.0, 1000.0]N0.001_margin200000.00har0.01cos100b10000_face_20170815_155910.txt';
%% half3400
% best till now
filename='20170816_132424/2mapres0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0_face_20170816_163437.txt';
%a combined loss: some totally flipped
% train vs test's comparing difference : 
% 	training and testing data are not that identical. for similar result, training is slightly better tha ntesting
filename='20170816_132936/2mapres0stn1har100000000N0[10000000.0, 0.0]N0.001_margin200000.0200000har0.01cos100000b10000000_face_20170816_184953.txt';
% left down: harmonic L1+area on real_half (if t is better, apply it) yes ,
% perhaps it is better
filename='20170817_020000/2mapl1res0stn1har100000N0[10000.0, 0.0]N0.001_margin200000.00har0.01cos0b0_face_20170817_034713.txt';
%combo loss, rectified /1+2/ dataset:worse than only area loss! ??? loss
%weight improper?
filename='20170817_011958/2mapres0stn1har100000000N0[10000000.0, 0.0]N0.001_margin200000.0200har0.01cos100000b10000000_face_20170817_074544.txt';
%left down but one: only area loss, rectified /1+2/ dataset
%relatively good
filename='20170817_012707/2mapres0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0_face_20170817_072633.txt';
% left first: har area, only /1/  good
filename='20170817_175111/2mapl1res0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.0200har0.01cos0b0_face_20170817_204940.txt';
%left 2nd: har area+boundary not better
%filename='20170817_175554/2mapl1res0stn1har0N0[10000000.0, 0.0]N0.001_margin200000.0200har0.01cos0b100000_face_20170817_204218.txt';
%left 1st target_pt='har',combo
filename='20170817_230824/2mapl1res0stn1har10000000N0[10000000.0, 0.0]N0.001_margin200000.0200har0.01cos100000b0_face_20170818_010704.txt';
%?
filename='20170818_204637/2mapl1res0stn1har10000000N0[10000000.0, 1000.0]N0.001_margin200000.0200har0.01cos100000b0_face_20170818_233754.txt';
%to check stacked mappers
%observation: loss decrease is postponed
filename='20170818_175341/smooth0.01har0N0[100000000.0, 0.0]N0.001_margin200000.00har0.01cos0b0_face_20170819_145213.txt'; 
% no har is better that has har
filename='20170819_215157/2mapl1res0stn1har0N0[10000000.0, 10000000.0]N0.001_margin200000.02e+11N20000cos100000b0_face_20170820_025426.txt';
% larger area_sum is better:best till now
filename='20170819_214303/2mapl1res0stn1har0N0[10000000.0, 10000.0]N0.001_margin200000.020000000000.0N20000cos100000b0_face_20170820_024744.txt';
%right down new dataset  : triangle like
filename='20170820_214454/2mapl1res0stn1har0N0[10000000.0, 10000000.0]N0.001_margin200000.020000000000.0N20000cos100000b0_face_20170820_233156.txt';

%right down: ft with sign_first(best till now) c #176
filename='20170820_115512/2mapl1res0stn1har0N0[10000000.0, 10000000.0]N0.001_margin200000.02e+11N20000cos100000b0_face_20170820_162710.txt';
%filename='20170820_214454/2mapl1res0stn1har0N0[10000000.0, 10000000.0]N0.001_margin200000.020000000000.0N20000cos100000b0_face_20170820_233156.txt';

filename='20170822_100212/2mapl1res0stn1har0N0[10000000.0, 10000.0]N0.001_margin200000.020000000000.0N20000cos100000b0_face_20170823_063338.txt';
%best till now for hard cases
filename='20170822_141050/2mapl1res0stn1har0N0[10000000.0, 10000.0]N0.001_margin200000.020000000000.0N20000cos100000b0_face_20170823_120556.txt';
%3d
% filename='20170822_105404/smooth13d1000000000000har0N0[100000000.0, 1000.0]N0.001_margin200000.00har0.01cos0b0_face_20170824_035308.txt';
% filename='20170821_131910/smooth13d100000000000har0N0[100000000.0, 1000.0]N0.001_margin200000.00har0.01cos0b0_face_20170823_192752.txt';
% filename='20170825_224953_/g_2dgres0har0N0[10000000000.0, 0.0]N0.001_margin200000.00.0N20000cos0b0_face_20170828_123415.txt';
map2=0;


%to make sure: bn before regression? relu between first conv and first
%gconv
dirname='../pointnet.pytorch/harmonic/results/';
pre=filename(17:end-25);
fidin=fopen(strrep([dirname filename],'face','input'));
fidface=fopen([dirname filename]);
fidout=fopen(strrep([dirname filename],'face','output'));
fidout2=fopen(strrep([dirname filename],'face','output2'));
fidhar=fopen(strrep([dirname filename],'face','har'));
has_har=0;
smooth=0;

% if has_har
%     fidhar=fopen([dirname pre '_har_' time '.txt']);
% end
% fidin=fopen([dirname pre '_input_' time '.txt']);
% fidout=fopen([dirname pre '_output_' time '.txt']);
% fidface=fopen([dirname pre '_face_' time '.txt']);
if has_har
har=textscan(fidhar,'%f%f', 'TreatAsEmpty','*','EmptyValue',10);
har=[har{1,1},har{1,2}];
end
in=textscan(fidin,'%f64 %f64 %f64');
in=[in{1,1},in{1,2},in{1,3}];
face=textscan(fidface,'%f %f %f');
face=[face{1,1},face{1,2},face{1,3}];
if smooth==1
    a=textscan(fidout,'%f%f%f', 'TreatAsEmpty','*','EmptyValue',10);
    a=[a{1,1},a{1,2},a{1,3}];
else % map2==1 ||
    a=textscan(fidout,'%f%f', 'TreatAsEmpty','*','EmptyValue',10);
    a=[a{1,1},a{1,2}];
end    
[r,c]=find(a==10);
[rr,~]=find(a(:,2)~=10);
last=1;
out=a;
if map2==1
    a2=textscan(fidout2,'%f%f', 'TreatAsEmpty','*','EmptyValue',10);
    out2=[a2{1,1},a2{1,2}];
end
% for i = 1:size(r,1)
%     out(r(i)+1:r(i)+250,:)=a(r(i)+1:r(i)+250,:)/100;
% end
% out=out(rr,:);
winsize=[600 50 2000 1300];
winsize=[1550 650 1330 1530];
fig = figure;u = fig.Units;fig.OuterPosition=(winsize); 
num_points=1868;%1840;%704; %700; %500 ;%1030;
num_face=3625;%3400;%1270; %1270; %1000; %2052;

labels = cellstr( num2str([1:num_points]') );
for i=1:num_points
    if mod(i,200)~=0
        labels{i}=' ';
    end
end

plotr=4;plotc=3+map2*2;
pos=-1;
for k =4:10:1000 %(size(in,1)/num_points)-37:(-3):1 %78:1:1000 % 
    %  356,237F,302face,409human,318 flower,355bear,%353T,284k,384Y,235F,250H,5 one
% 5three,10 plane4two,24chess,172T,125C
    
    inc=1;
    pos=pos+1;
    a_ori=in((k-1)*num_points+1:k*num_points,:);
    b_ori=out((k-1)*num_points+1:k*num_points,:);
    tri=face((k-1)*num_face+1:k*num_face,:);  
    if has_har
    b_har=har((k-1)*num_points+1:k*num_points,:);
    end
    if map2==1
        b_ori2=out2((k-1)*num_points+1:k*num_points,:);
    end
    %% process face
    l=0;
    for j=100:num_face
        if tri(j,:)==zeros(1,3)
           l=j-1;
           break;
        end
    end
    if l==0
        l=num_face;
    end
    tri=tri(1:l,:);
    tri=tri+1;
    a=a_ori(1:max(max(tri)),:);
    b2=b_ori(1:max(max(tri)),:);
    if map2==1
            b=b_ori2(1:max(max(tri)),:);
    else
        b=b2;
    end
    if map2==1
        tmp=b;
        b=b2;b2=tmp;
    end
    if has_har
        b_har=b_har(1:max(max(tri)),:);    
    end
    %% calculate area

    
%     for j=1:size(Y,1)
%         colormap(I(j),:)=[j/size(Y,1),0.5*j/size(Y,1),1-j/size(Y,1)];
%     end
 %%   figure plot
   subplot(plotr,plotc,pos*plotc+inc);inc=inc+1;

    hold on, trisurf(tri,a(:,1),a(:,2),a(:,3),'FaceColor','interp','FaceAlpha',1), hold off %,'FaceColor','interp'
        axis equal;axis tight;
     text(a(1:max(max(tri)),1),a(1:max(max(tri)),2),a(1:max(max(tri)),3), labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
                           'HorizontalAlignment','right');
    title(['3D Mesh' num2str(k) ]);
    [E] = exterior_edges(tri);
    s = a(E(:,1),:); e = a(E(:,2),:);
    line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
    hold off;
    axis equal; axis tight; %axis off; 
    ori_a=a;

%% plot
    if has_har
        subplot(plotr,plotc,pos*plotc+inc);inc=inc+1;  
        z=[1:max(max(tri))]';
        hold on, trisurf(tri,b_har(1:max(max(tri)),1),b_har(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.6,'EdgeColor','interp','EdgeAlpha',0.2), 
        hold off
        %hold on, trisurf(tri,b(1:max(max(tri)),1),b(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,

        %          text(a(1:max(max(tri)),1), a(1:max(max(tri)),2),z, labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
    %                                   'HorizontalAlignment','right');
            title(['Harmonic' ]);         % set(gca,'Color','k');
    %  subplot(plotr,plotc,pos*plotc+3);  
         axis equal;axis tight;
    end
%     grid on;
%     hold on, trisurf(tri,b_har(1:max(max(tri)),1),b_har(1:max(max(tri)),2),(b(:,1)-b_har(:,1)),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,
% 
%      subplot(plotr,plotc,pos*plotc+4);  
%     axis equal;axis tight;
%     hold on, trisurf(tri,b_har(1:max(max(tri)),1),b_har(1:max(max(tri)),2),(b(:,2)-b_har(:,2)),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,



    [color_map2,area2]=area_ratio(a,b,tri);
    if has_har
    [color_map,area]=area_ratio(a,b_har,tri);
   subplot(plotr,plotc,pos*plotc+inc);inc=inc+1;
    z=[1:max(max(tri))]';
    z=ones(size(b,1),1);   

    hold on, trisurf(tri,b_har(1:max(max(tri)),1),b_har(1:max(max(tri)),2),flipud(z),'FaceVertexCData',color_map,'EdgeColor','none'), hold off %,'FaceColor','interp'  ,'FaceAlpha',0.6   
    axis equal;axis tight; title(['area:' num2str(area)]);
    grid on;    
    end
    b=b(1:max(max(tri)),:);
    area1=0.5*doublearea(b,tri);
    flip_num=sum(area1<0);
    flip_area=sum((area1<0).*area1);
%    close all;inc=1;
   subplot(plotr,plotc,pos*plotc+inc);inc=inc+1;  
    axis equal;axis tight;
    grid on;
    z=[1:max(max(tri))]';
    if smooth==0
    hold on, trisurf(tri,b(1:max(max(tri)),1),b(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,
          text(b(1:max(max(tri)),1), b(1:max(max(tri)),2),z, labels(1:max(max(tri))), 'VerticalAlignment','bottom', ...
                                   'HorizontalAlignment','right');
             ti=sprintf('Flip num: %d/%d,',flip_num,size(tri,1) );    title(ti);         set(gca,'Color','k');
    else
    hold on, trisurf(tri,b(:,1),b(:,2),b(:,3),'FaceColor','interp','FaceAlpha',1), hold off %,'FaceColor','interp'
    s = b(E(:,1),:); e = b(E(:,2),:);
    line([s(:,1)';e(:,1)'], [s(:,2)';e(:,2)'], [s(:,3)';e(:,3)'], 'Color', 'r', 'LineWidth', 3);
    end

 %   title(['Ellipsoid Map #Point: ' num2str(num_points) ', #Face: ' num2str(num_face) ])


    axis equal; axis tight; %axis off; 
    if smooth==0
   subplot(plotr,plotc,pos*plotc+inc);inc=inc+1;
%     area2=doublearea(b2,tri);
%     flip_num2=sum(area2<0);
%     flip_area2=sum((area2<0).*area2);    
    hold on, trisurf(tri,b(1:max(max(tri)),1),b(1:max(max(tri)),2),flipud(z),'FaceVertexCData',color_map2,'EdgeColor','none'), hold off %,'FaceColor','interp'  ,'FaceAlpha',0.6
    ti=sprintf('ori/map area: %.3f/%.2f',flip_area,area2 );    title(ti);      
%     hold on, trisurf(tri,b2(1:max(max(tri)),1),b2(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.5,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,
%     ti2=sprintf('2nd map: %d(%.3f)',flip_num2,flip_area2 );    title(ti2);  set(gca,'Color','k');
    axis equal;axis tight;        
    end
%     axis equal;
%     grid on;
%     a=b;
%     z=[1:max(max(tri))]';
%     hold on, trisurf(tri,a(1:max(max(tri)),1)*100,a(1:max(max(tri)),2)*100,zeros(max(max(tri)),1),'FaceColor','interp','FaceAlpha',0.3), hold off
    axis equal; axis tight; %axis off; 
    if map2==1
       subplot(plotr,plotc,pos*plotc+inc);inc=inc+1; 
             
%    hold on, trisurf(tri,b2(1:max(max(tri)),1),b2(1:max(max(tri)),2),b2(:,3),'FaceColor','interp'), hold off %,'FaceColor','interp'  ,'FaceAlpha',0.6
    hold on, trisurf(tri,b2(1:max(max(tri)),1),b2(1:max(max(tri)),2),ori_a(:,3),'FaceColor','interp','FaceAlpha',0.6,'EdgeColor','interp','EdgeAlpha',0.2), hold off %,'FaceColor','interp'  ,'FaceAlpha',0.6
    area1=0.5*doublearea(b2,tri);
    flip_num=sum(area1<0);
    flip_area=sum((area1<0).*area1);
    ti=sprintf('Flip num: %d/%d,',flip_num,size(tri,1) );    title(ti);     
    axis equal;axis tight;    
    
    [color_map2,area2]=area_ratio(b2,b,tri);        
   subplot(plotr,plotc,pos*plotc+inc);inc=inc+1;
    hold on, trisurf(tri,b2(1:max(max(tri)),1),b2(1:max(max(tri)),2),flipud(z),'FaceVertexCData',color_map2,'EdgeColor','none'), hold off %,'FaceColor','interp'  ,'FaceAlpha',0.6
    ti=sprintf('smoothed/map: %.3f/%.2f',flip_area,area2 );    title(ti);      
    axis equal;axis tight;  
    end
if pos==plotr-1   %  fail25,41,45,52  hard47
        pos=-1;
        suptitle(strrep(pre,'_',' '));
        w = waitforbuttonpress;
        if w == 0
            disp('Button click')
            else
            disp('Key press')
        end
        fig = figure;u = fig.Units;fig.OuterPosition=(winsize); 
   end
end
end

function area=triangle_area(pts,tri)

if size(pts,2)==2
    pts=[pts,zeros(size(pts,1),1)];
end
vin1=pts(tri(:,1),:);
vin2=pts(tri(:,2),:);
vin3=pts(tri(:,3),:);
s12=vin2-vin1 ;
s13=vin3-vin1 ; 
crossp=cross(s12', s13');
area=0.5*sum(crossp.^2,1);
    
area=area';
end
% calculate area change
function [color_map,area]=area_ratio(a,b,tri)

    in_area=0.5*doublearea(a,tri);
    out_area=0.5*doublearea(b,tri);
    ratio=in_area./out_area;
    [~,I]=sort(ratio);
    lib=jet(size(tri,1));
    color_map=zeros(size(tri));
    color_map(I,:)=lib;
    area=sum((out_area));
end

%%%
%% cubicinterp
% results/test_input_20170713_102559.txt
% dirname='../results/';
% time='20170713_102559';
% pre='test';
%trained by 3d_ellipsoid
% dirname='../../3d_ellipsoid/results/20170713_134849/';
% time='1';
% dirname='../results/20170713_153202/';
% time='20170713_155618';
% pre='test2';
% pre='cubicinterptest';
% % add mse of area +marginarea+cos
%     0.35z:first good result
% % dirname='../../cubicinterp/results/20170713_225654/';
% % time='20170713_232848';
% % pre='test_area';
%     randhalf dataset :good result  flatter the surface from inside the
%     sllipsoid(conclusion from >5 images),zero input has the same output
%     for each mesh
% dirname='../results/20170714_094322/';
% time='20170714_101639';
% pre='randhalf_area';
%     1z:mostly failed  perhaps need to panish small area,need to find the
%     boundary between failure and success in terms of lambda in lambda*z
%dirname='../results/20170714_095204/';
% time='20170714_102320';
% pre='cubicinterp1z_area';
% trained by v1map inputtrans_noglobal_3d_ellipsoid model bad
% dirname='../../3d_ellipsoid/results/';
% time='20170714_093259';
% pre='randhalfcos_marginarea';
% observation: squash along y (y+ or y-)

%partialy failed, only those parallel( or slightly tilt) to the x-y plane are the most
%successful ones =>perhaps the input transform is not good enough
% dirname='../results/20170714_153813/';  
% time='20170714_163057';
% pre='0.6z_area';

% dirname='../results/20170714_131129/';  %failure 25,29(the hole boundary is the most distorted part,but without inputtransform,this 29 looks just as normal as other notransform failure case,
%perhaps its failure is due to inputtransform's failure)
%perhaps they have learned the ability of keeping the consistensy for the
%same hole shape with different hole directions.(bu observing the mapped hole boundary's point order)
% time='20170714_131131';
% pre='over_0.8_area';

%trained with over_0.8,but tested with 0.6z noinputtransform, partially
%developed
% dirname='../results/20170714_214444/';  
% time='20170714_214446';
% pre='over_0.8_area';

%notrans    squash from the front view&from x- to x+;successful cases are
%mostly from this direction=>the input transform module learns to rotate
%the surface so that the hole is on the left
% dirname='../results/20170714_234706/';  
% time='20170714_234709';
% pre='over_0.8_area';
%   o2i1 mean fail25,41,45,52  hard47, has fixed case 29 in the last successful
%   one
% dirname='../results/20170717_105349/';  
% time='20170717_105350';
% pre='over_0.8_weightedarea_mean';
%o2i1 sum
% dirname='../results/20170717_105554/';
% time='20170717_105556';
% pre='over_0.8_weightedarea_mean';

% dirname='../../real10k/results/20170717_154803/';
% time='20170717_154804';
% pre='over_0.8_weightedarea_sum_real10k';

% dirname='../../real10k/results/20170717_193958/';
% time='20170717_200147';
% pre='bumpy_weightedarea_sum';
% 
% dirname='../../real10k/results/20170720_113734/';
% time='20170720_113736';
% pre='bumpy_weightedarea_sum';

%% trained by real 1k
%noweightedarea  most fail
% dirname='../pointnet.pytorch/real/results/20170731_235109/';
% time='20170801_023401';
% pre='real_weightedarea_sum_margin';
%most fail
% dirname='../pointnet.pytorch/real/results/20170801_001635/';
% time='20170801_072949';
% pre='real_weighted1area_sum_margin_cos0';

% same mesh, different cut: little difference
% 2 slightly different mesh: greatly different
% dirname='../pointnet.pytorch/real/results/20170801_002636/';
% time='20170801_030044';
% pre='real_weighted0area_sum_margin_cos1';


%% half real
% overlapped
% dirname='../pointnet.pytorch/half_real/results/20170802_172124/';
% time='20170802_174631';
% pre='real_weighted0_summargin1barrier0area_cos10';

% => barrier better than margin loss (at least the former is trying to flatten it), but some tends to be a line
% dirname='../pointnet.pytorch/half_real/results/20170802_172946/';
% time='20170802_175836';
% pre='real_weighted0_summargin0barrier1area_cos10';

%cos larger cos better than not, but the easiest failed(cylinder worse than
%yuanzhui),human is good
% dirname='../pointnet.pytorch/half_real/results/20170802_173311/';
% time='20170802_180148';
% pre='real_weighted0_summargin0barrier1area_cos100000';

% add barrier loss:become one line ==>cosine loss not functioning
% mostly not flipped
% dirname='../pointnet.pytorch/half_real/results/20170802_213816/';
% time='20170802_234513';
% pre='real_20weighted1000000000_summargin0barrier2000area_cos100000';
%barrier v2 is better,at least not like a line, but most human case
%failed(folded)
% dirname='../pointnet.pytorch/half_real/results/20170802_231636/';
% time='20170803_005227';
% pre='real_20weightedv21000000000_summargin0barrier2000area_cos100000';


% wrong normalization method
% (dirname='../pointnet.pytorch/half_real/results/20170803_023135/';
% time='20170803_045419';
% pre='norm_20weightedv21000000000_summargin0barrier2000area_cos100000';
% %results//_face_
% 
% dirname='../pointnet.pytorch/half_real/results/20170803_025551/';
% time='20170803_051445';
% pre='norm_20000weightedv21000000_summargin0wbarrier2000area_cos100000';)

% next step: barrier+cos?


%+normalize

%not good :has kept the abs(area), but many flipped
% time='20170803_131055';
% dirname='../pointnet.pytorch/half_real/results/20170803_110406/';
% pre='norm_2000weightedv21000000_summargin0wbarrier20000area_cos100000';

%good eg 121&131(hand) no flipping on the boundary
% time='20170803_132339';
% dirname='../pointnet.pytorch/half_real/results/20170803_110806/';
% pre='norm_20weightedv21000000000_summargin0wbarrier2000area_cos100000';



%+residual
%better than the next one in that less flip at guaidian, but worse in that
%the area in guaidian is quite small(but on the other hand, this particular 
%area property on guaidian preserves the property of curvature somehow)
% eg 41&61 to compare these two model
%convergence is very fast
% dirname='../pointnet.pytorch/half_real/results/20170803_110853/';
% time='20170803_132414';
% pre='res_norm_20weightedv21000000000_summargin0wbarrier2000area_cos100000';

%#1 xyz area (this one is actually with res(but forgot to include 'res' in its name))   Good! (best till now)
%failure 1:cannot let the face too big at guaidian, even for boundary  eg.two, three
% when flipping happens, decrease the area loss; greater difference for
% xyz;
% failure 2 overlap of two non-nerighboring parts  eg.630
% time='20170803_183844';
% dirname='../pointnet.pytorch/half_real/results/20170803_170951/';
% pre='xyz21000000000.0boundary20barrier20002000.0boundary0cos100000';
%xyz(no res)
%almost as good as the last one=>residual makes no difference,res only set
%the mapping's x-y coordinate among all possible rotations
%by zooming in, we can observe that the overlapped area is quite small
% dirname='../pointnet.pytorch/half_real/results/20170803_211504/';
% time='20170803_224312';
% pre='res0xyz21000000000boundary20_barrier20002000boundary0cos100000';


%assume the input faces are mostly the same size, differences among faces
%in the mapping will also indicate lots of info w.r.t. the area change of
%each face
% area of the highest faces are larger than #1
% dirname='../pointnet.pytorch/half_real/results/20170804_025604/';
% time='20170804_051902';
% pre='res0xyz41000000000boundary15_barrier20002000boundary0cos100000';

%compared to the last one, this one has larger infinite for BarrierArea (previous's basepara is 100000(not 2000))
%=> larger barrier loss (vs area loss is better)
%it preserves more original info of the points farthest to origin  e.g.
%four && larger (high face area)/(low face area) ratio
%less twirled e.g. case 111 (compared to last version)
%
%SSolution to do: boundary to curved ranking as input or supervision (perpoint)
% set the ranking fixed in a range of points
%SSolution: larger barrier loss is better
%dirname='../pointnet.pytorch/half_real/results/20170804_030107/';
% time='20170804_052355';
% pre='res0xyz41000000000boundary15_barrier200000002000boundary0cos100000';
% 
%till now typical failure: 201,one,two,three,four
% 
% dirname='../pointnet.pytorch/siamese/results/20170804_015620/';
% time='20170804_053155';
% pre='test_res0xyz21000000000boundary20_barrier20002000boundary0cos100000sia100000';

%%cos1 vs cos100000 makes little difference
% dirname='../pointnet.pytorch/half_real/results/20170805_134937/';
% time='20170805_155020';
% pre='res0xyz41000000000boundary15_barrier200000002000boundary0cos1';

%the following two failed
% (dirname='../pointnet.pytorch/half_real/results/20170806_172224/';
% time='20170806_213507';
% pre='signfirst0xyz41000000000boundary15_margin2000002000000000000000N0cos100000';
% 
% dirname='../pointnet.pytorch/half_real/results/20170806_171106/';
% pre='signfirst0xyz41000000000boundary15_margin2000002000000000000000boundary0cos100000';
% time='20170806_213749';)

%% 3d_ellipsoid
% folded
% time='20170708_154521'
% pre='mse_area'
% squeeze the ellipsoid to a flat surface
% marginarea_cos_input_20170708_204814.txt

% dirname='../results/'
% time='20170709_120330'
% pre='mse_cos'
% dirname='../results/'
% time='20170708_144725'
% pre='mse_cos'

%%%

