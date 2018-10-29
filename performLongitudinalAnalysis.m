function performLongitudinalAnalysis

% Parameters
NSUR = 10000;
MIN_ASSESSMENTS = [3:18];

TARGET = [6,5,7];
FEATURES = 8:75;
TITLES = {'PHQ9','GAD7','Suicide'};

% Load model trained on group-level data
load('model')

% Read data from table
data = readtable('groundtruth_individualdata.csv');
targets = table2array(data(:,TARGET));
features = table2array(data(:,FEATURES));

% Compute number of assessments each participant completed
pid = unique(data.pid);
for p = 1:length(pid)
    num_assessments(p) = sum(data.pid == pid(p));
end

% Fit model to individual data
for nt = 1:length(MIN_ASSESSMENTS)
    fprintf('.')
    
    % Find participant that completed the minimal number of assessments
    pid1 = pid(num_assessments>=MIN_ASSESSMENTS(nt));
    npid = length(pid1);

    for p = 1:npid
        
        % select rows
        i = find(data.pid == pid1(p));
        
        for t = 1:length(TARGET)
            for n = 1:NSUR
                i = i(randperm(length(i),MIN_ASSESSMENTS(nt)));
                
                % select data
                x = zscore(features(i,:));
                target = zscore(targets(i,t));
                
                % estimate predicted scores using regression model
                yfit = x*model(t).beta;
                
                % correlate with measured scores
                R=corrcoef(yfit,target);
                
                r(n) = R(2,1);
                mse0(n) = immse(ones(size(target))*mean(target),target);
                mse1(n) = immse(yfit,target);
            end
            results.r(t,p,nt) = median(r);
            results.mse0(t,p,nt) = median(mse0);
            results.mse1(t,p,nt) = median(mse1);
        end
    end
end
fprintf('\n')


% Plot results
colours = cbrewer('qual','Set1',3);

set(gcf,'units','centimeters','position',[0 30 8.9 11])

results.r(results.r==0) = NaN;
for x = 1:3
    subplot('position',[0.09 1.055-x*0.325,0.86,0.23])
    plot([0.3 length(MIN_ASSESSMENTS)+0.7],[0 0],'k--')
    hold on
    UnivarScatter(squeeze(results.r(x,:,:)),'PointSize',14,'MarkerFaceColor',colours(x,:),...
        'SEMColor',[0.85, 0.85 0.85],'StdColor',[0.85, 0.85 0.85],'Width',0.5);
    hold off
    box off
    set(gca,'fontsize',8,'xlim',[0.3 length(MIN_ASSESSMENTS)+0.7],'xticklabels',MIN_ASSESSMENTS,'ytick',[-1,0,1])
    title(TITLES{x},'fontsize',9)
    ylabel('r','fontsize',9)
    if x == 3
        xlabel('number of assessments','fontsize',9)
    end
end