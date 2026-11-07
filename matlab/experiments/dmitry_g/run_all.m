function run_all(limit,selection,bg1K)
% Три моделирования: одна точка Eb/N0 на 12 процессов через parfor.
% run_all — расчёт 1e6; run_all(1e3) — пробный расчёт.
% run_all(1e6,2) — только второе моделирование.
% run_all(1e6,1:3,528) — BG1 K=528 вместо K=3872.
% run_all(1e6,1,528) — только сравнение BG1 K=528 и BG2 K=520.
if nargin<1, limit=1e6; end
if nargin<2, selection=1:3; end
if nargin<3, bg1K=3872; end
assert(isnumeric(bg1K) && isscalar(bg1K) && ismember(bg1K,[528 3872]), ...
    'Для BG1 поддерживаются K=528 и K=3872.');
if bg1K==528, bg1Zc=24; else, bg1Zc=176; end
assert(isnumeric(limit) && isscalar(limit) && isfinite(limit) && limit>=1 && limit==floor(limit));
assert(isnumeric(selection) && ~isempty(selection) && all(ismember(selection,1:3)));
selection=unique(selection(:).');
base=fileparts(mfilename('fullpath'));
core=fullfile(base,'..','..');
addpath(base,core);
assert(exist('ldpc_simulate','file')==2 && exist('ldpc_default_params','file')==2, ...
    'Не найдены общие функции проекта в папке matlab.');
set(groot,'defaultFigureVisible','off');
if limit==1e6, folderName='results'; else, folderName=sprintf('results_%.0f',limit); end
out=cell(1,3);
for model=selection
    modelFolder=folderName;
    % Для K=528 отдельная папка только первого моделирования.
    % Моделирования 2 и 3 используют прежние общие результаты BG2.
    if model==1 && bg1K==528, modelFolder=[folderName '_k528']; end
    out{model}=fullfile(base,sprintf('моделирование %d',model),modelFolder);
    if ~isfolder(out{model}), mkdir(out{model}); end
end
p=ldpc_default_params();p.maxFrames=limit;
p.stopAtZeroFer=false; % Рассчитываем весь диапазон, даже после точки без ошибок.
sets={make_algorithms(p,1,bg1Zc),make_algorithms(p,2,52), ...
    make_iterations(p),make_termination(p)};
owners=[1 1 2 3];
labelsA={'BP (flooding)','Layered BP','NMS \alpha=0.75','OMS \beta=0.5','Bit-flip'};
labelsI=arrayfun(@(n) sprintf('I_{max}=%d',n),[5 10 20 50 100],'UniformOutput',false);
labelsT={'early (синдром = 0)','max (всегда 50 итераций)'};
% Используем приложенные совместимые MAT; недостающие точки досчитываются.
inputNames={{'e1b_bg1_bitflip.mat','bg1_k3872_checkpoint.mat'}, ...
    {'e1b_bg2_bitflip.mat'},{'e2_bp_iterations.mat'},{'e2_termination.mat'}};
if bg1K==528, inputNames{1}={'e1b_bg1_bitflip.mat'}; end
for group=1:4
    model=owners(group);
    if ~ismember(model,selection),continue;end
    for fileIndex=1:numel(inputNames{group})
        file=fullfile(out{model},inputNames{group}{fileIndex});
        if ~isfile(file),continue;end
        C=load(file,'results');
        assert(isfield(C,'results'),'Нет results в MAT: %s',file);
        sets{group}=import_series(sets{group},C.results,limit);
    end
end

% Возобновляем только результаты этой версии и этой папки. Старые допы не нужны.
for model=selection
    file=fullfile(out{model},'checkpoint.mat');
    if ~isfile(file), continue; end
    C=load(file,'state','limit','schemeVersion');
    assert(isfield(C,'schemeVersion') && C.schemeVersion==1 && C.limit==limit, ...
        'Несовместимый checkpoint: %s. Переименуйте папку результатов.',file);
    ids=find(owners==model);assert(numel(C.state)==numel(ids));
    for j=1:numel(ids)
        expected=sets{ids(j)};cached=C.state{j};
        assert(numel(cached)==numel(expected));
        for i=1:numel(expected)
            assert(isequaln(cached(i).params,expected(i).params) && ...
                cached(i).K==expected(i).K && isequal(cached(i).EbN0dB,expected(i).EbN0dB), ...
                'Параметры изменились. Переименуйте папку результатов: %s',out{model});
            check_result(cached(i),limit);
        end
        sets{ids(j)}=cached;
    end
end
jobs=struct('set',{},'curve',{},'point',{},'params',{});
total=0;done=0;
for n=1:4
    if ~ismember(owners(n),selection), continue; end
    for i=1:numel(sets{n})
        r=sets{n}(i);total=total+numel(r.EbN0dB);done=done+nnz(r.numFrames>0);
        for k=1:numel(r.EbN0dB)
            if r.numFrames(k)>0, continue; end
            q=r.params;q.EbN0dB=r.EbN0dB(k);q.seed=q.seed+k-1;q.verbose=false;
            jobs(end+1)=struct('set',n,'curve',i,'point',k,'params',q); %#ok<AGROW>
        end
    end
end
fprintf('Готово %d/%d точек. Осталось %d. Максимум кадров на точку: %g.\n', ...
    done,total,numel(jobs),limit);
if ~isempty(jobs)
    assert(license('test','Distrib_Computing_Toolbox'),'Нужен Parallel Computing Toolbox.');
    workers=12;
    pool=gcp('nocreate');
    if isempty(pool)
        cluster=parcluster('Processes');cluster.NumWorkers=workers;
        pool=parpool(cluster,workers);
    else
        assert(isa(pool,'parallel.ProcessPool') && pool.NumWorkers==workers, ...
            ['Нужен пул Processes из 12 процессов. Текущий пул автоматически не закрывается. ' ...
             'Завершите его вручную и повторите запуск.']);
    end
    attachments={which(mfilename),which('ldpc_simulate'), ...
        which('ldpc_bitflip_decode'),which('nr_ldpc_pcm')};
    assert(all(~cellfun(@isempty,attachments)),'Не найдены необходимые файлы декодера.');
    missing=attachments(~ismember(attachments,pool.AttachedFiles));
    if ~isempty(missing),addAttachedFiles(pool,missing);end
    setup=parfevalOnAll(pool,@addpath,0,base,core);
    fetchOutputs(setup);
    for j=1:numel(jobs)
        job=jobs(j);model=owners(job.set);
        fprintf('Моделирование %d: BG%d K=%d, %s, maxIter=%d, %.2f дБ.\n', ...
            model,job.params.bgn,sets{job.set}(job.curve).K,job.params.algorithm, ...
            job.params.maxIter,job.params.EbN0dB);
        one=compute_point(pool,job.params,fullfile(out{model},'point_in_progress.mat'));
        r=sets{job.set}(job.curve);
        assert(one.K==r.K && one.numFrames>0 && one.numFrames<=limit);
        for key={'ber','fer','avgIter','numFrames','numBitErrors','numFrameErrors','timePerFrame'}
            r.(key{1})(job.point)=one.(key{1});
        end
        check_result(r,limit);
        sets{job.set}(job.curve)=r;
        save_model_checkpoint(out{model},sets(owners==model),limit);
        fprintf('[%d/%d] Готово, BER=%.3e, FER=%.3e.\n', ...
            done+j,total,one.ber,one.fer);
    end
end
for n=1:4
    if ismember(owners(n),selection)
        assert(all(arrayfun(@(r) all(r.numFrames>0),sets{n})));
    end
end
for model=selection
    save_model_checkpoint(out{model},sets(owners==model),limit);
end
samplingInfo='Imported points retain original samples; new points use frame chunks across 12 processes and wave-level stopping';
if ismember(1,selection)
    A=sets{1}; B=sets{2}; labels=labelsA;
    results=A; save(fullfile(out{1},'e1b_bg1_bitflip.mat'),'results','labels','samplingInfo');
    results=B; save(fullfile(out{1},'e1b_bg2_bitflip.mat'),'results','labels','samplingInfo');
    results=A(1:4); labels=labelsA(1:4);
    save(fullfile(out{1},'e1_bg1_algorithms.mat'),'results','labels','samplingInfo');
    results=B(1:4); save(fullfile(out{1},'e1_bg2_algorithms.mat'),'results','labels','samplingInfo');
    plot_metrics(A(1:4),labelsA(1:4),sprintf('BG1 K=%d: мягкие декодеры',bg1K),fullfile(out{1},'e1_bg1_algorithms.png'));
    plot_metrics(A,labelsA,sprintf('BG1 K=%d',bg1K),fullfile(out{1},'e1_bg1_all.png'));
    plot_metrics(B(1:4),labelsA(1:4),'BG2 K=520: мягкие декодеры',fullfile(out{1},'e1_bg2_algorithms.png'));
    plot_metrics(B,labelsA,'BG2 K=520',fullfile(out{1},'e1_bg2_all.png'));
    for metric={'ber','fer'}
        plot_pair(A,B,labelsA,metric{1},fullfile(out{1},['bg1_vs_bg2_' metric{1} '.png']));
    end
    write_raw({A,B},fullfile(out{1},'bg_raw_measurements.csv'));
    write_thresholds({A,B},labelsA,fullfile(out{1},'bg_comparison.csv'));
end
if ismember(2,selection)
    results=sets{3}; labels=labelsI; iterList=[5 10 20 50 100];
    save(fullfile(out{2},'e2_bp_iterations.mat'),'results','labels','iterList','samplingInfo');
    plot_metrics(results,labels,'BP flooding: Eb/N0 до 2,5 дБ',fullfile(out{2},'e2_bp_iterations.png'));
    plot_one(results,labels,'ber','BP flooding: влияние maxIter',fullfile(out{2},'e2_bp_iterations_ber.png'));
    write_raw({results},fullfile(out{2},'e2_bp_iterations.csv'));
end
if ismember(3,selection)
    results=sets{4}; labels=labelsT;
    save(fullfile(out{3},'e2_termination.mat'),'results','labels','samplingInfo');
    plot_metrics(results,labels,'BP flooding: early и max',fullfile(out{3},'e2_termination.png'));
    write_raw({results},fullfile(out{3},'e2_termination.csv'));
end
fprintf('Готово. Результаты сохранены в папках выбранных моделирований; BG1 K=%d.\n',bg1K);
end

function rr=import_series(rr,cached,limit)
assert(numel(rr)==numel(cached),'Неполный набор кривых во входном MAT.');
for i=1:numel(rr)
    p=rr(i).params;q=cached(i).params;
    ignored={'verbose','stopAtZeroFer','EbN0dB'};
    assert(rr(i).K==cached(i).K && isequal(rr(i).EbN0dB,cached(i).EbN0dB), ...
        'Другая длина блока или сетка Eb/N0 во входном MAT.');
    assert(isequaln(rmfield(p,ignored),rmfield(q,ignored)), ...
        'Параметры входного MAT несовместимы с текущим расчётом.');
    check_result(cached(i),limit);
    for key={'ber','fer','avgIter','numFrames','numBitErrors','numFrameErrors','timePerFrame'}
        rr(i).(key{1})=cached(i).(key{1});
    end
end
end

function save_model_checkpoint(folder,state,limit)
file=fullfile(folder,'checkpoint.mat');tmp=[file '.tmp.mat'];schemeVersion=1;
samplingInfo='Frame-parallel Monte Carlo, chunk-derived seeds, aggregate wave-level stopping';
save(tmp,'state','limit','schemeVersion','samplingInfo');
[ok,msg]=movefile(tmp,file,'f');assert(ok,'Не удалось сохранить checkpoint: %s',msg);
end

function check_result(r,limit)
ok=r.numFrames>0;
assert(all(r.numFrames>=0 & r.numFrames<=limit));
assert(all(isfinite(r.ber(ok))) && all(isfinite(r.fer(ok))));
assert(all(r.numFrameErrors(ok)>=0 & r.numFrameErrors(ok)<=r.numFrames(ok)));
assert(all(r.numBitErrors(ok)>=0 & r.numBitErrors(ok)<=r.numFrames(ok)*r.K));
assert(all(abs(r.ber(ok)-r.numBitErrors(ok)./(r.numFrames(ok)*r.K))<1e-12));
assert(all(abs(r.fer(ok)-r.numFrameErrors(ok)./r.numFrames(ok))<1e-12));
assert(all(r.avgIter(ok)>=0 & r.avgIter(ok)<=r.params.maxIter));
end

function C=read_optional(file)
C=struct();if ~isfile(file),return;end
% Частичный прогресс этой версии не игнорируем молча при повреждении MAT.
C=load(file);
end

function r=compute_point(pool,p,file)
% Все процессы считаются для одной точки. Прогресс сохраняется после каждой волны.
version=1;
C=read_optional(file);
if isfield(C,'progress') && isfield(C.progress,'version') && ...
        C.progress.version==version && isequaln(C.progress.params,p)
    progress=C.progress;
else
    progress=struct('version',version,'params',p,'numFrames',0, ...
        'numBitErrors',0,'numFrameErrors',0,'sumIter',0,'decoderSeconds',0, ...
        'nextChunk',1,'startedAt',char(datetime('now')));
end
workers=pool.NumWorkers;
while progress.numFrames<p.maxFrames && progress.numFrameErrors<p.minFrameErrors
    % В начале малыми порциями оцениваем частоту ошибок, затем берём до 512 кадров.
    chunkFrames=p.batchSize;
    if progress.numFrames>0
        if progress.numFrameErrors==0
            chunkFrames=512;
        else
            need=(p.minFrameErrors-progress.numFrameErrors) / ...
                (workers*progress.numFrameErrors/progress.numFrames);
            chunkFrames=max(p.batchSize,min(512,ceil(need/p.batchSize)*p.batchSize));
        end
    end
    remaining=p.maxFrames-progress.numFrames;
    n=min(workers,ceil(remaining/chunkFrames));
    counts=min(chunkFrames,remaining-(0:n-1)*chunkFrames);
    assert(all(counts>0) && sum(counts)<=remaining);
    ids=progress.nextChunk+(0:n-1);
    chunks=cell(1,n);
    parfor w=1:n
        chunks{w}=simulate_chunk(p,counts(w),ids(w));
    end
    % Складываем в порядке номеров порций, независимо от порядка завершения процессов.
    for w=1:n
        one=chunks{w};
        assert(one.numFrames==counts(w));
        progress.numFrames=progress.numFrames+one.numFrames;
        progress.numBitErrors=progress.numBitErrors+one.numBitErrors;
        progress.numFrameErrors=progress.numFrameErrors+one.numFrameErrors;
        progress.sumIter=progress.sumIter+round(one.avgIter*one.numFrames);
        progress.decoderSeconds=progress.decoderSeconds+one.timePerFrame*one.numFrames;
    end
    progress.nextChunk=progress.nextChunk+n;
    assert(progress.numFrames<=p.maxFrames);
    tmp=[file '.tmp.mat']; save(tmp,'progress');
    [ok,msg]=movefile(tmp,file,'f');assert(ok,'Не удалось сохранить прогресс: %s',msg);
    fprintf('  Кадров %d/%d, ошибочных кадров %d/%g, процессов в волне %d.\n', ...
        progress.numFrames,p.maxFrames,progress.numFrameErrors,p.minFrameErrors,n);
end
r=empty_result(p);
r.numFrames=progress.numFrames;
r.numBitErrors=progress.numBitErrors;
r.numFrameErrors=progress.numFrameErrors;
r.ber=r.numBitErrors/(r.numFrames*r.K);
r.fer=r.numFrameErrors/r.numFrames;
r.avgIter=progress.sumIter/r.numFrames;
r.timePerFrame=progress.decoderSeconds/r.numFrames;
end

function r=simulate_chunk(p,count,id)
% У порций разные детерминированные seed; одинаковая порция не дублируется на процессах.
maxNumCompThreads(1);
q=p; q.maxFrames=count; q.minFrameErrors=Inf; q.stopAtZeroFer=false; q.verbose=false;
% Новая схема выборки: она не воспроизводит побитово прежний последовательный поток.
q.seed=mod(double(p.seed)*1000003+double(id)*9176,2^32-2);
r=ldpc_simulate(q);
end

function rr=make_algorithms(p,bgn,Zc)
algs={'bp','layered-bp','norm-min-sum','offset-min-sum','bit-flip'}; rr=[];
for i=1:5
    q=p; q.bgn=bgn; q.Zc=Zc; q.maxIter=50; q.termination='early';
    q.algorithm=algs{i}; q.EbN0dB=0:0.25:2;
    if i==5, q.EbN0dB=3:0.5:9; end
    rr=[rr,empty_result(q)]; %#ok<AGROW>
end
end

function rr=make_iterations(p)
rr=[];
for it=[5 10 20 50 100]
    q=p; q.bgn=2; q.Zc=52; q.algorithm='bp'; q.maxIter=it;
    q.termination='early'; q.EbN0dB=0:0.25:2.5;
    rr=[rr,empty_result(q)]; %#ok<AGROW>
end
end

function rr=make_termination(p)
rr=[]; terms={'early','max'};
for i=1:2
    q=p; q.bgn=2; q.Zc=52; q.algorithm='bp'; q.maxIter=50;
    q.termination=terms{i}; q.EbN0dB=0:0.25:2;
    rr=[rr,empty_result(q)]; %#ok<AGROW>
end
end

function r=empty_result(p)
if p.bgn==1, K=22*p.Zc; N=66*p.Zc; else, K=10*p.Zc; N=50*p.Zc; end
n=numel(p.EbN0dB); E=round(K/p.rate);
r.EbN0dB=p.EbN0dB(:).';
r.ber=nan(1,n); r.fer=nan(1,n); r.avgIter=nan(1,n);
r.numFrames=zeros(1,n); r.numBitErrors=zeros(1,n); r.numFrameErrors=zeros(1,n);
r.timePerFrame=nan(1,n);
r.codeName=sprintf('5G NR BG%d Zc=%d K=%d R=%g',p.bgn,p.Zc,K,K/E);
r.K=K; r.N=N; r.E=E; r.params=p;
end

function draw_curve(ax,r,metric,label,i)
markers={'o','s','d','^','v'};
colors=lines(5);
v=r.(metric); ok=r.numFrames>0 & isfinite(v);
if ~strcmp(metric,'avgIter'), v(v==0)=NaN; end
plot(ax,r.EbN0dB(ok),v(ok),['-' markers{mod(i-1,5)+1}], ...
    'Color',colors(mod(i-1,5)+1,:),'LineWidth',1.6,'MarkerSize',6,'DisplayName',label);
end

function plot_metrics(rr,labels,ttl,file)
fig=figure('Color','w','Position',[70 70 1600 650]);
cleanup=onCleanup(@() close(fig));
tl=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
metrics={'ber','fer','avgIter'}; yl={'BER','FER','Среднее число итераций'};
axesList=gobjects(1,3);
for m=1:3
    ax=nexttile(tl); axesList(m)=ax;
    hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.FontSize=12;
    if m<3, ax.YScale='log'; end
    for i=1:numel(rr), draw_curve(ax,rr(i),metrics{m},labels{i},i); end
    xlabel(ax,'E_b/N_0, дБ'); ylabel(ax,yl{m});
end
title(tl,ttl,'FontSize',17);
lg=legend(axesList(1),'Location','southoutside','Orientation','horizontal', ...
    'NumColumns',numel(labels),'FontSize',12);
lg.Layout.Tile='south';
exportgraphics(fig,file,'Resolution',180);
end

function plot_one(rr,labels,metric,ttl,file)
fig=figure('Color','w','Position',[100 100 900 680]);
cleanup=onCleanup(@() close(fig));
tl=tiledlayout(fig,1,1,'Padding','compact');
ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.YScale='log'; ax.FontSize=13;
for i=1:numel(rr), draw_curve(ax,rr(i),metric,labels{i},i); end
xlabel(ax,'E_b/N_0, дБ'); ylabel(ax,upper(metric)); title(ax,ttl);
xlim(ax,[min(rr(1).EbN0dB),max(rr(1).EbN0dB)]);
lg=legend(ax,'Location','southoutside','Orientation','horizontal','NumColumns',numel(labels));
lg.Layout.Tile='south';
exportgraphics(fig,file,'Resolution',180);
end

function plot_pair(left,right,labels,metric,file)
fig=figure('Color','w','Position',[80 80 1600 680]);
cleanup=onCleanup(@() close(fig));
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
first=[];
for j=1:2
    if j==1, rr=left; else, rr=right; end
    ax=nexttile(tl); if j==1, first=ax; end
    hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.YScale='log'; ax.FontSize=13;
    for i=1:numel(rr), draw_curve(ax,rr(i),metric,labels{i},i); end
    xlim(ax,[0 9]); ylim(ax,[1e-6 1]);
    xlabel(ax,'E_b/N_0, дБ'); ylabel(ax,upper(metric));
    title(ax,sprintf('BG%d, K=%d, Zc=%d',rr(1).params.bgn,rr(1).K,rr(1).params.Zc));
end
title(tl,sprintf('BPSK / AWGN, R=1/2, maxIter=50, maxFrames=%g',left(1).params.maxFrames),'FontSize',16);
lg=legend(first,'Location','southoutside','Orientation','horizontal','NumColumns',5,'FontSize',13);
lg.Layout.Tile='south';
exportgraphics(fig,file,'Resolution',180);
end

function write_raw(sets,file)
rows=[];
for c=1:numel(sets)
    rr=sets{c};
    for i=1:numel(rr)
        r=rr(i); n=numel(r.EbN0dB);
        rows=[rows; repmat([r.params.bgn,r.K,i],n,1),r.EbN0dB(:), ...
            r.ber(:),r.fer(:),r.avgIter(:),r.numFrames(:), ...
            r.numBitErrors(:),r.numFrameErrors(:)]; %#ok<AGROW>
    end
end
writetable(array2table(rows,'VariableNames', ...
    {'BG','K','algorithm','EbN0dB','BER','FER','avgIter','numFrames','numBitErrors','numFrameErrors'}),file);
end

function write_thresholds(sets,labels,file)
metric={}; target=[]; bg=[]; K=[]; algorithm={}; threshold=[];
targets={'ber',1e-2;'ber',1e-3;'ber',1e-4;'fer',1e-2};
for c=1:numel(sets)
    rr=sets{c};
    for t=1:size(targets,1)
        for i=1:numel(rr)
            metric{end+1,1}=targets{t,1}; target(end+1,1)=targets{t,2}; %#ok<AGROW>
            bg(end+1,1)=rr(i).params.bgn; K(end+1,1)=rr(i).K; %#ok<AGROW>
            algorithm{end+1,1}=labels{i}; %#ok<AGROW>
            threshold(end+1,1)=ebno_at_target(rr(i),targets{t,1},targets{t,2}); %#ok<AGROW>
        end
    end
end
writetable(table(metric,target,bg,K,algorithm,threshold, ...
    'VariableNames',{'metric','target','BG','K','algorithm','EbN0dB'}),file);
end

function eb = ebno_at_target(r, metric, target)
% Первое пересечение заданного уровня сверху вниз между соседними рассчитанными точками; логарифмическая интерполяция.
v = r.(metric); eb = NaN;
for k = 1:numel(v)-1
    if all(r.numFrames(k:k+1)>0) && all(isfinite(v(k:k+1))) && ...
            all(v(k:k+1)>0) && v(k)>=target && v(k+1)<=target
        if v(k)==v(k+1), eb=r.EbN0dB(k); return; end
        f=(log10(target)-log10(v(k)))/(log10(v(k+1))-log10(v(k)));
        eb=r.EbN0dB(k)+f*(r.EbN0dB(k+1)-r.EbN0dB(k)); return
    end
end
end
