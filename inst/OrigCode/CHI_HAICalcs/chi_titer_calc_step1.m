demo_import = importdata('chi_data/flu-demographics-db-export-coh2.txt');
demo.gender = demo_import.textdata(2:end,1);
demo.race = demo_import.textdata(2:end,2);
demo.sample_id = demo_import.data(:,3);
demo.primary_id = demo_import.data(:,2);
demo.age = demo_import.data(:,1);
demo.gender_int = ones( size(demo.gender) );
demo.gender_int( strcmp( demo.gender, 'F') ) = 2;
%%
d70 = importdata('chi_data/day70_d0_d70match_08_2013_raw.txt');
d0 = importdata('chi_data/day0_d0_d70match_08_2013_raw.txt');
assert(sum(d70.data(:,1)-d0.data(:,1))==0);
%%
primary_ids = demo.primary_id;

mask = ismember(d70.data(:,1), primary_ids);

d0_raw.uru = d0.data(mask,4);
d0_raw.bri = d0.data(mask,3);
d0_raw.swi = d0.data(mask,2);
d0_raw.b_bri = d0.data(mask,5);

d70_raw.uru = d70.data(mask,4);
d70_raw.bri = d70.data(mask,3);
d70_raw.swi = d70.data(mask,2);
d70_raw.b_bri = d70.data(mask,5);

% get the primary ids of the subjects
primary_ids = d70.data(mask,1);

%% get fc_max
use_log = 0;

r = compute_titer_fc(d0_raw.uru, d70_raw.uru, use_log);
r2 = compute_titer_fc(d0_raw.bri, d70_raw.bri, use_log);
r3 = compute_titer_fc(d0_raw.swi, d70_raw.swi, use_log);
r4 = compute_titer_fc(d0_raw.b_bri, d70_raw.b_bri, use_log);

% standardize using median and mad
r.fc_mad_norm = mad_standardize(r.fc);
r2.fc_mad_norm = mad_standardize(r2.fc);
r3.fc_mad_norm = mad_standardize(r3.fc);
r4.fc_mad_norm = mad_standardize(r4.fc);

% set infinity points to NaN if needed
r.fc_mad_norm( isinf(r.fc_mad_norm) ) = NaN;
r2.fc_mad_norm( isinf(r2.fc_mad_norm) ) = NaN;
r3.fc_mad_norm( isinf(r3.fc_mad_norm) ) = NaN;
r4.fc_mad_norm( isinf(r4.fc_mad_norm) ) = NaN;

% set fc_max - here fc_max is the max of the mad/median normalized numbers
% across the four viruses
fc_max = nanmax( [r.fc_mad_norm r2.fc_mad_norm r3.fc_mad_norm r4.fc_mad_norm], [], 2);
% inverse_normalize fc_max (so that only the rank counts)
inv_norm_fc_max=inv_normal_transform2_tied(fc_max);

% process day0 titer; only use the seasonal 
uru_d0 = mad_standardize(d0_raw.uru);
bri_d0 = mad_standardize(d0_raw.bri);
b_bri_d0 = mad_standardize(d0_raw.b_bri);
swi_d0 = mad_standardize(d0_raw.swi);
ub_max = max( [ uru_d0 bri_d0 b_bri_d0 ], [], 2 );

%% after running john\compute_titer_1114.m (2 first cells)
CHI.subject = primary_ids;
% [CHI.subject, sidx] = sort(CHI.subject);
CHI.timeptn = [0;70];
CHI.virus = {'A/Uruguay/716/07_H3N2','A/Brisbane/59/07_H1N1','A/California/07/2009_H1N1','B/Brisbane/60/2008'}';
d70 = cell2mat(struct2cell(d70_raw)');
d0 = cell2mat(struct2cell(d0_raw)');
CHI.hai = permute(cat(3,d0,d70),[1 3 2]);
CHI.age = [];
CHI.fc = d70./d0;
CHI.d0 = d0;
CHI.fc_norm = [r.fc_mad_norm, r2.fc_mad_norm, r3.fc_mad_norm, r4.fc_mad_norm];
CHI.d0_norm = [uru_d0, bri_d0, swi_d0, b_bri_d0];
CHI.d0_max = nanmax(CHI.d0, [], 2);
CHI.fc_max = nanmax(CHI.fc, [], 2);
CHI.fc_4fc = ones(size(CHI.fc_max));
CHI.fc_4fc(all(CHI.fc<4,2)) = 0;
CHI.fc_4fc( sum(CHI.fc>=4,2) > floor(numel(CHI.virus)*2/3) ) = 2;
CHI.d0_norm_max = ub_max;
CHI.fc_norm_max = fc_max;
CHI.fc_norm_max_ivt = inv_norm_fc_max;

%% get flow
fdata = importdata('chi_data\flow_coh2_day0.log.agrcorr.txt');
flow_subject = fdata.data(1,:);
[cidx, fidx] = ismember(CHI.subject, flow_subject);
fidx(fidx==0) = [];

%% filter subjects by d0 flow
fnames = fieldnames(CHI);
nsubj = numel(CHI.subject);
for f=1:numel(fnames)
    if size(CHI.(fnames{f}),1) ~= nsubj
        continue
    end
    CHI.(fnames{f}) = CHI.(fnames{f})(cidx,:,:);
end

%% get age
[cidx,didx] = ismember(CHI.subject,demo.primary_id);
assert(all(cidx),'Missing subject(s)')
CHI.age = demo.age(didx);

%% save
save CHI_170213
