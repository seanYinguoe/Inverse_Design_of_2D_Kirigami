function result = genetic_search(evaluate,x0,lower,upper,cfg,checkpoint)
%GENETIC_SEARCH Bounded real-valued GA: tournament selection and one elite.
% A serial implementation needs no Global Optimization/Parallel Toolbox.
if nargin<6,checkpoint=@(~) [];end
validateattributes(cfg.population,{'numeric'},{'scalar','integer','>=',2});
validateattributes(cfg.generations,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(cfg.crossover,{'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(cfg.mutation,{'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(cfg.stall_generations,{'numeric'},{'scalar','integer','positive'});
assert(all(isfinite(x0))&&all(lower<=x0)&&all(x0<=upper),'Invalid seed or bounds.');
N=cfg.population;D=numel(x0);
pop=lower+rand(N,D).*(upper-lower);pop(1,:)=x0;
fitness=inf(N,1);details=cell(N,1);failures={};
% The seed must solve. Do not hide connection/licence/configuration failures.
[fitness(1),details{1}]=evaluate(x0);
assert(isfinite(fitness(1))&&isreal(fitness(1))&&fitness(1)>=0,'Invalid baseline fitness.');
if fitness(1)>cfg.tolerance
    for k=2:N,[fitness(k),details{k}]=candidate(pop(k,:));end
end
history=[];stall=0;previous=inf;
for generation=0:cfg.generations
    [best,order]=sort(fitness);b=order(1);
    history(end+1,:)=[generation best(1) mean(fitness(isfinite(fitness)))]; %#ok<AGROW>
    result=struct('x',pop(b,:),'fitness',best(1),'evaluation',details{b}, ...
        'history',history,'failures',{failures},'generation',generation,'stop_reason','');
    fprintf('Generation %d: best shape error %.6g\n',generation,best(1));
    if previous-best(1)>1e-8,stall=0;else,stall=stall+1;end
    previous=best(1);
    if best(1)<=cfg.tolerance,result.stop_reason='target tolerance';
    elseif generation>=cfg.generations,result.stop_reason='generation limit';
    elseif stall>=cfg.stall_generations,result.stop_reason='stalled';end
    checkpoint(result);
    if ~isempty(result.stop_reason),break;end
    next=zeros(size(pop));next(1,:)=pop(b,:);
    for k=2:N
        a=tournament();c=tournament();child=pop(a,:);
        if rand<cfg.crossover
            mix=rand(1,D);child=mix.*child+(1-mix).*pop(c,:);
        end
        mutate=rand(1,D)<cfg.mutation;
        scale=0.2*(1-generation/(cfg.generations+1));
        child=child+mutate.*randn(1,D).*scale.*(upper-lower);
        next(k,:)=max(lower,min(upper,child));
    end
    pop=next;fitness(1)=best(1);details{1}=result.evaluation;
    for k=2:N,[fitness(k),details{k}]=candidate(pop(k,:));end
end
    function index=tournament
        pair=randi(N,1,2);[~,winner]=min(fitness(pair));index=pair(winner);
    end
    function [f,detail]=candidate(x)
        try
            [f,detail]=evaluate(x);
            assert(isfinite(f)&&isreal(f)&&f>=0,'Invalid candidate fitness.');
        catch e
            f=inf;detail=[];
            failures{end+1}=struct('x',x,'identifier',e.identifier,'message',e.message);
            message=strsplit(e.message,newline);
            fprintf('  Candidate rejected: %s\n',message{1});
        end
    end
end
