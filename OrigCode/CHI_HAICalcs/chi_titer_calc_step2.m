%%
load('CHI_170213.mat')

%% studies to use
fnames = {'CHI'};
cntr = {'CHI'}; % center
lbl = strcat(fnames,'_',cntr);
study = 'CHI';
%%
fnames = [fnames, strcat(fnames,'_y')];
cntr = repmat(cntr,1,2);
lbl = strcat(fnames,'_',cntr);
%%
ageth = 36;
idx_y = CHI.age < ageth;
idx_o = CHI.age > ageth;
idx_a = true(size(CHI.age));

%% split young and old
study = 'CHI';
for yo = 'y'
idx = eval(strcat('idx_',yo));
newstudy = strcat(study,'_',yo);
ts = CHI;
ts.subject = ts.subject(idx);
ts.age = ts.age(idx);
ts.hai = ts.hai(idx,:,:);
ts.fc = ts.fc(idx,:);
ts.d0 = ts.d0(idx,:);
ts.fc_norm = ts.fc_norm(idx,:);
ts.d0_norm = ts.d0_norm(idx,:);
% ts.fc_std_norm = ts.fc_std_norm(idx,:);
% ts.d0_std_norm = ts.d0_std_norm(idx,:);
ts.fc_norm_max = ts.fc_norm_max(idx,:);
ts.fc_norm_max_ivt = ts.fc_norm_max_ivt(idx,:);
ts.d0_norm_max = ts.d0_norm_max(idx,:);
ts.d0_max = ts.d0_max(idx,:);
ts.fc_max = ts.fc_max(idx,:);
ts.fc_4fc = ts.fc_4fc(idx,:);
ts.fc_norm_max_d20 = ts.fc_norm_max_d20(idx,:);
ts.fc_norm_max_d30 = ts.fc_norm_max_d30(idx,:);
% ts = rmfield(ts,{'table','fc_mad_norm','d0_mad_norm'});
% ts = rmfield(ts,{'table'});
eval(strcat(newstudy,' = ts;'));
end

%% fc_res_max
d0_bin{1} = {[-10 1], [1.1 50]};
d0_bin{2} = {[-10 1], [1.1 50]};

[CHI.fc_res_max2, CHI.bininfo] = decorrelate_by_bin_idx(CHI.d0_norm_max, CHI.fc_norm_max_ivt, d0_bin{1});
[CHI_y.fc_res_max2, CHI_y.bininfo] = decorrelate_by_bin_idx(CHI_y.d0_norm_max, CHI_y.fc_norm_max_ivt, d0_bin{2});

%%
CHI_y.fc_norm_max_d20 = discretize(CHI_y.fc_norm_max,20,20);
CHI_y.fc_norm_max_d30 = discretize(CHI_y.fc_norm_max,30,30);
CHI_y.d0_bin = {[-10 1], [1.1 50]};
bin=20;
rt = decorrelate_by_bin(CHI_y.d0_norm_max, CHI_y.fc_norm_max_ivt, bin, bin, CHI_y.d0_bin);
CHI_y.fc_res_max = rt.fc_res;
CHI_y.fc_res_max_d20 = discretize(CHI_y.fc_res_max,20,20);
CHI_y.fc_res_max_d30 = discretize(CHI_y.fc_res_max,30,30);

%% checking before decorrelation
s = 2;
% study = fnames{s};
% idx = idx_a;
ts = CHI_y;
idx = true(size(ts.age));
    if min(ts.age(idx)) >= ageth
        agelbl = 'old';
%         xlim(xl)
    elseif max(ts.age(idx)) < ageth
        agelbl = 'young';
%         xlim(xl)
    else
        agelbl = 'all';
%         xl = xlim;
    end
x = ts.d0_norm_max(idx);
y = ts.fc_norm_max_ivt(idx);
[cc, pv] = corr( x, y, 'type', 'Spearman', 'rows', 'pairwise');
clf
densitybubbleplot( x, y, 500, 'b', 0 ); 
title([untex(lbl{s}),', ',agelbl, ', n=',num2str(sum(idx)),', spearman corr=' num2str(cc,'%.2f'), ', p=', num2str(pv,'%.1e')])
    xlabel('d0\_norm\_max')
    ylabel('fc\_norm\_max\_ivt')
    ylim([-2.5 2.5])
%     xlim([-1 8])
hv = vline(unique(cell2mat(ts.d0_bin)));
set(hv,'LineWidth',2)
textonaxis(gca, sprintf('bin %d: n=%d, p=%.3f\n',[1:numel(ts.bininfo.binN);ts.bininfo.binN;ts.bininfo.pv]),'se');
% xlim(xl)

%% and save
set(gcf,'color','w')
filename = sprintf('%s_fc_max_ivt_vs_d0_norm_max_before',study);
export_fig(filename,'-png')
%% checking after decorrelation
% idx = idx_a;
idx = true(size(ts.age));
    if min(ts.age(idx)) > ageth
        agelbl = 'old';
%         xlim(xl)
    elseif max(ts.age(idx)) < ageth
        agelbl = 'young';
%         xlim(xl)
    else
        agelbl = 'all';
%         xl = xlim;
    end
x = ts.d0_norm_max(idx);
y = ts.fc_res_max(idx);
[cc, pv] = corr( x, y, 'type', 'Spearman', 'rows', 'pairwise');
clf
densitybubbleplot( x, y, 500, 'b', 0 ); 
title([untex(lbl{s}),', ',agelbl, ', n=',num2str(sum(idx)),', spearman corr=' num2str(cc), ', p=', num2str(pv)])
    xlabel('d0\_norm\_max')
    ylabel('fc\_res\_max')
%     ylim([-2.5 2.5])
%     xlim([-1 7])
vline(unique(cell2mat(ts.d0_bin)))
%% and save
set(gcf,'color','w')
filename = sprintf('%s_fc_max_ivt_vs_d0_norm_max_after',study);
export_fig(filename,'-png')

%% discretize fc_res_max
bin = [20 30];
    for b = bin
        study = fnames{s};
        CHI_y.(sprintf('fc_res_max_d%d',b)) = discretize(CHI_y.fc_res_max,b,b);
    end

%% output study to single file
fn = {'age','fc','d0','d0_norm_max','fc_norm_max','fc_norm_max_ivt','d0_max','fc_max',...
    'fc_4fc','fc_norm_max_d20','fc_norm_max_d30','fc_res_max','fc_res_max_d20','fc_res_max_d30'};

% for s = 1:numel(fnames)
    study = 'SDY80_young';
    t = CHI_y;
    out = table(t.subject,'VariableNames',{'subject'});
    for f = 1:numel(fn)
        if size(t.(fn{f}),2) > 1
            a = num2cell(t.(fn{f}),1);
            avar = strcat(fn{f},'_',t.virus);
            if verLessThan('matlab','8.3.0')
                avar = regexprep(avar,'\W','_');
            else
                avar = matlab.lang.makeValidName(avar);
            end
            out = [out, table( a{:}, 'VariableNames',{avar{:}} )];
        else
            out = [out, table( t.(fn{f}), 'VariableNames',fn(f) )];
        end
    end
    out = sortrows(out, 'subject');
    filename = sprintf('%s_hai_titer_table.txt',study);
    writetable(out,filename,'delimiter','\t')
% end

