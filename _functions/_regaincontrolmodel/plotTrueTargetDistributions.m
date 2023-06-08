addpath(genpath('~/chris-lab/code_general/'));

% load a stimulus set
root = '~/chris-lab/projects/gain_behavior/_stimuli';
t1 = load(fullfile(root,'testingWide_200406hiloChord-bluerec1-dual.mat'),'DB','params');
t2 = load(fullfile(root,'testingWide_200406lohiChord-bluerec1-dual.mat'),'DB','params');

% target info
tSamps = t1.params.noiseD / t1.params.chordDuration;
tInd = 2;

% for each volume, get the target volume over all unique times and trials
for i = 1:size(t1.DB,4)
   
    cnt = 1;
   
    for j = 1:size(t1.DB,2)
        for k = 1:size(t1.DB,3)
                       
            noise1(cnt) = mean(t1.DB{1,j,k,i}(:,tSamps(j)+1),1);
            noise2(cnt) = mean(t2.DB{1,j,k,i}(:,tSamps(j)+1),1);
            vol1(i,cnt) = mean(t1.DB{tInd,j,k,i}(:,tSamps(j)+1),1);
            vol2(i,cnt) = mean(t2.DB{tInd,j,k,i}(:,tSamps(j)+1),1);
            cnt = cnt + 1;
           
        end
    end
   
end

edges = linspace(45,65,100);
lStr{1} = 'noise';
for i = 1:size(t1.DB,4)
    lStr{i+1} = sprintf('+%d dB',t1.params.targetDBShift(i));
end
xf = linspace(edges(1),edges(end),1000);

f1 = figure(1); clf;
subplot(2,1,1);
hold on
h(1) = histogram(noise1,edges,'Normalization','pdf');
h(1).FaceColor = [.5 .5 .5];
%h(1).LineWidth = 1.5;
for i = 1:size(t1.DB,4)
    h(i+1) = histogram(vol1(i,:),edges,'Normalization','pdf');
    h(i+1).FaceAlpha = .1*i;
    h(i+1).FaceColor = 'b';
    pd1 = fitdist(vol1(i,:)','Normal');
    pp = plot(xf,pdf(pd1,xf'),'k','LineWidth',1);
    %pp.Color(4) = 1;
end
text(.01,.9,sprintf('sigma = %3.2f\ntrue width = %d',pd1.sigma,t1.params.sd(2)),...
                    'units','normalized');
legend(h,lStr,'location','northeastoutside');
title('Low Contrast Targets')
ylabel('Trial Count')
axis tight; plotPrefs;

subplot(2,1,2);
hold on
h(1) = histogram(noise2,edges,'Normalization','pdf');
%h(1).LineWidth = 1.5;
h(1).FaceColor = [.5 .5 .5];
for i = 1:size(t1.DB,4)
    h(i+1) = histogram(vol2(i,:),edges,'Normalization','pdf');
    h(i+1).FaceAlpha = .1*i;
    h(i+1).FaceColor = 'r';
    pd2 = fitdist(vol2(i,:)','Normal');
    plot(xf,pdf(pd2,xf'),'k','LineWidth',1);
end
text(.01,.9,sprintf('sigma = %3.2f\ntrue width = %d',pd2.sigma,t2.params.sd(2)),...
     'units','normalized');
text(.01,.7,sprintf('sigma ratio = %3.2f\ntrue ratio = %d',pd2.sigma/pd1.sigma,...
                    t2.params.sd(2)/t1.params.sd(2)),'units','normalized');
legend(h,lStr,'location','northeastoutside');
title('High Contrast Targets')
xlabel('Mean Target Volume (dB)')
ylabel('Trial Count')
axis tight; plotPrefs;

saveFigPDF(f1,[600 300],'./target_distributions.pdf');