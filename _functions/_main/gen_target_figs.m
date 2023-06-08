addpath(genpath('~/chris-lab/code_general/'))

% resample to 1000Hz
r = 200;

%% HILO psych
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/testingWide_200406hiloChord-bluerec1-dual.mat';
load(fn)

fs = params.fs / r;

colormap('default');

f1 = figure(1);
for i = 1:6
    
    subplot(2,1,1)
    clm = params.mu + max(params.sd).*[-1 1];
    t = (0:params.chordDuration:(length(DB{2,2,1,i})*params.chordDuration)) ...
        - 3;
    imagesc(t,1:length(params.freqs),DB{2,2,1,i},clm); colorbar;
    octs = find(mod(params.freqs,1000)==0);
    set(gca,'ydir','normal'); % flip it to match frequency order
    set(gca,'ytick',octs);
    set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    plotPrefs; box off;
    xlim([-3 1.2])
    
    subplot(2,1,2)
    t = 1/fs:1/fs:4.2;
    I = 1:4.2*params.fs;
    y = stimf{2,2,1,i}(I);
    plot(t,y(1:r:end),'k') ;
    xlim([0 4.2]); ylim([-.25 .25]);
    colorbar;
    plotPrefs; box off;


    
    saveFigPDF(f1,[500 250],sprintf('./hilo_vol%d.pdf',i));
    
end

f1 = figure(1); clf;
subplot(2,1,1)
clm = params.mu + max(params.sd).*[-1 1];
t = (0:params.chordDuration:(length(DB{1,2,1,1})*params.chordDuration)) ...
    - 3;
imagesc(t,1:length(params.freqs),DB{1,2,1,1},clm); colorbar;
octs = find(mod(params.freqs,1000)==0);
set(gca,'ydir','normal'); % flip it to match frequency order
set(gca,'ytick',octs);
set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
plotPrefs; box off;
xlim([-3 1.2])

subplot(2,1,2)
t = 1/fs:1/fs:4.2;
I = 1:4.2*params.fs;
y = stimf{1,2,1,i}(I);
plot(t,y(1:r:end),'k') ;
xlim([0 4.2]); ylim([-.25 .25]);
colorbar;
plotPrefs; box off;

saveFigPDF(f1,[500 250],sprintf('./hilo_noiseOnly.pdf',i));


%% LOHI psych
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/testingWide_200406lohiChord-bluerec1-dual.mat';
load(fn)

f1 = figure(1); clf;
for i = 1:6
        
    subplot(2,1,1)
    clm = params.mu + max(params.sd).*[-1 1];
    t = (0:params.chordDuration:(length(DB{2,2,1,i})*params.chordDuration)) ...
        - 3;
    imagesc(t,1:length(params.freqs),DB{2,2,1,i},clm); colorbar;
    octs = find(mod(params.freqs,1000)==0);
    set(gca,'ydir','normal'); % flip it to match frequency order
    set(gca,'ytick',octs);
    set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    plotPrefs; box off;
    xlim([-3 1.2])
    
    subplot(2,1,2)
    t = 1/fs:1/fs:4.2;
    I = 1:4.2*params.fs;
    y = stimf{2,2,1,i}(I);
    plot(t,y(1:r:end),'k') ;
    xlim([0 4.2]); ylim([-.25 .25]);
    colorbar;
    plotPrefs; box off;
    
    saveFigPDF(f1,[500 250],sprintf('./lohi_vol%d.pdf',i));
    
end

f1 = figure(1); clf;
subplot(2,1,1)
clm = params.mu + max(params.sd).*[-1 1];
t = (0:params.chordDuration:(length(DB{1,2,1,1})*params.chordDuration)) ...
    - 3;
imagesc(t,1:length(params.freqs),DB{1,2,1,1},clm); colorbar;
octs = find(mod(params.freqs,1000)==0);
set(gca,'ydir','normal'); % flip it to match frequency order
set(gca,'ytick',octs);
set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
plotPrefs; box off;
xlim([-3 1.2])

subplot(2,1,2)
t = 1/fs:1/fs:4.2;
I = 1:4.2*params.fs;
y = stimf{1,2,1,i}(I);
plot(t,y(1:r:end),'k') ;
xlim([0 4.2]); ylim([-.25 .25]);
colorbar;
plotPrefs; box off;

saveFigPDF(f1,[500 250],sprintf('./lohi_noiseOnly.pdf',i));








%% HILO OFFSET
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/200406_offsetshiloChord-bluerec1-thresh9.783890e+00-dual.mat';
load(fn)

f1 = figure(1); clf;
for i = 1:size(DB,2)
    
    subplot(2,1,1)
    clm = params.mu + max(params.sd).*[-1 1];
    t = (0:params.chordDuration:(length(DB{2,i,1,2})*params.chordDuration)) ...
        - 3;
    imagesc(t,1:length(params.freqs),DB{2,i,1,2},clm); colorbar;
    octs = find(mod(params.freqs,1000)==0);
    set(gca,'ydir','normal'); % flip it to match frequency order
    set(gca,'ytick',octs);
    set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    plotPrefs; box off;
    xlim([-3 1.2])
    
    subplot(2,1,2)
    t = 1/fs:1/fs:min(4.2,length(stimf{2,i,1,2})/params.fs);
    I = 1:min(4.2*params.fs,length(stimf{2,i,1,2}));
    y = stimf{2,i,1,2}(I);
    plot(t,y(1:r:end),'k') ;
    xlim([0 4.2]); ylim([-.25 .25]);
    colorbar;
    plotPrefs; box off;
    
    saveFigPDF(f1,[500 250],sprintf('./hilo_time%d.pdf',i));
    
end


%% LOHI OFFSET
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/200305_offsetslohiChord-bluerec1-thresh13.4985-dual.mat';
load(fn)

f1 = figure(1);
for i = 1:5
        
    subplot(2,1,1)
    clm = params.mu + max(params.sd).*[-1 1];
    t = (0:params.chordDuration:(length(DB{2,i,1,2})*params.chordDuration)) ...
        - 3;
    imagesc(t,1:length(params.freqs),DB{2,i,1,2},clm); colorbar;
    octs = find(mod(params.freqs,1000)==0);
    set(gca,'ydir','normal'); % flip it to match frequency order
    set(gca,'ytick',octs);
    set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
    xlabel('Time (s)'); ylabel('Frequency (kHz)');
    plotPrefs; box off;
    xlim([-3 1.2])
    
    subplot(2,1,2)
    t = 1/fs:1/fs:min(4.2,length(stimf{2,i,1,2})/params.fs);
    I = 1:min(4.2*params.fs,length(stimf{2,i,1,2}));
    y = stimf{2,i,1,2}(I);
    plot(t,y(1:r:end),'k') ;
    xlim([0 4.2]); ylim([-.25 .25]);
    colorbar;
    plotPrefs; box off;
    
    saveFigPDF(f1,[500 250],sprintf('./lohi_time%d.pdf',i));
    
end



%% TARGETS ONLY - PSYCH
% plot just the sound wave, but with each different target overlaid
% in a different color (plot noise only first, then loudest to
% softest)
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/testingWide_200406hiloChord-bluerec1-dual.mat';
load(fn)
f1 = figure(1); clf; hold on;
t = 1/fs:1/fs:4.2;
I = 1:4.2*params.fs;
s0 = stimf{1,2,1,1}(I);
y0 = s0(1:r:end);
plot(t,y0,'k') ;

%c = cool;
%colors = c(round(linspace(1,64,6)),:);
colors = cool(6);
colormap(colors);

for i = 6:-1:1
    
    % find index of diff (ie. target points)
    si = stimf{2,2,1,i}(I);
    yi = si(1:r:end)
    di = find((yi - y0) ~= 0);
    
    % plot it with a different color
    plot(t(di),yi(di),'color',colors(i,:),'linewidth',1);
    
end
plot([3 3],ylim,'r--','linewidth',1)
ch = colorbar;
ch.Ticks = linspace(1/12,1-1/12,6);
ch.TickLabels = num2str(params.targetDBShift');
ch.TickDirection = 'out'
ct = get(ch,'Label');
set(ct,'String','Target Volume (dB SNR)')
xlim([2.75 4.2]); ylim([-.25 .25]);
plotPrefs; box off;

saveFigPDF(f1,[250 150],'./hilo_psych.pdf');



% plot just the sound wave, but with each different target overlaid
% in a different color (plot noise only first, each time)
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/200406_offsetshiloChord-bluerec1-thresh9.783890e+00-dual.mat';
load(fn)

nOff = size(stimf,2);

f1 = figure(1); clf; hold on;
t = 1/fs:1/fs:length(stimf{2,nOff,1,2})/params.fs;
I = 1:length(stimf{2,nOff,1,2});
s0 = stimf{1,nOff,1,2}(I);
y0 = s0(1:r:end);
plot(t,y0,'k') ;

color = [ones(nOff,1) linspace(.8,0,nOff)' ones(nOff,1)];
colormap(color);

for i = 1:nOff
    
    % new index for stim length
    si = stimf{2,i,1,2};
    yi = si(1:r:end);
    I = 1:length(yi)-.2*fs;
    
    % find index of diff (ie. target points)
    di = find(abs((yi(I) - y0(I))) > 0.001);
    
    % plot it with a different color
    plot(t(di),yi(di),'color',color(i,:),'linewidth',1);
    
end
plot([3 3],ylim,'r--','LineWidth',1)
ch = colorbar;
ch.Ticks = linspace(1/(2*nOff),1-1/(2*nOff),nOff)
ch.TickLabels = num2str((params.noiseD'-3)*1000);
ch.TickDirection = 'out'
ct = get(ch,'Label');
set(ct,'String','Target Time (ms)')
set(gca,'xtick',params.noiseD)
%xtickangle(45)
xlim([2.75 4.2]); ylim([-.25 .25]);
plotPrefs; box off;

saveFigPDF(f1,[250 150],'./hilo_offset.pdf');



%% generate figures for background vs no background
fn = '/Users/chris/chris-lab/projects/gain_behavior/_stimuli/testingWide_200406lohiChord-bluerec1-dual.mat';
load(fn)

f1 = figure(1); clf;
i = 6

subplot(6,1,[1 2])
clm = params.mu + max(params.sd).*[-1 1];
t = (0:params.chordDuration:(length(DB{2,2,1,i})*params.chordDuration)) ...
    - 3;
imagesc(t,1:length(params.freqs),DB{2,2,1,i},clm); colorbar;
octs = find(mod(params.freqs,1000)==0);
set(gca,'ydir','normal'); % flip it to match frequency order
set(gca,'ytick',octs);
set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
plotPrefs; box off;
xlim([-1 1]);

subplot(6,1,3)
t = 1/fs:1/fs:4.2;
I = 1:4.2*params.fs;
y = stimf{2,2,1,i}(I);
plot(t,y(1:r:end),'k') ;
xlim([2 4]); ylim([-.3 .3]);
colorbar;
plotPrefs; box off;


subplot(6,1,[4 5])
amps = DB{2,2,1,i};
t = (0:params.chordDuration:(length(amps)*params.chordDuration)) ...
    - 3;
amps(:,t >= 0 & t< .5) = 0;
amps(:,find(t >= .525)-1:end) = 0;
imagesc(t,1:length(params.freqs),amps,clm); colorbar;
octs = find(mod(params.freqs,1000)==0);
set(gca,'ydir','normal'); % flip it to match frequency order
set(gca,'ytick',octs);
set(gca,'yticklabels',num2str(params.freqs(octs)'/1000));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
plotPrefs; box off;
xlim([-1 1]);

subplot(6,1,6)
t = 1/fs:1/fs:4.2;
I = 1:4.2*params.fs;
y = stimf{2,2,1,i}(I);
y = y(1:r:end);
y(t >= 3 & t< 3.5) = 0;
y(t >= 3.525) = 0;
plot(t,y,'k') ;
xlim([2 4]); ylim([-.3 .3]);
colorbar;
plotPrefs; box off;


    saveFigPDF(f1,[250 250],'./lohi_silence_control.pdf');