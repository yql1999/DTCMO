classdef DTCMO < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Evolutionary multitasking-based constrained multiobjective optimization
%------------------------------- Reference --------------------------------
% EMCMO

    methods
        function main(Algorithm,Problem)
            Population{1} = Problem.Initialization(); % 考逐个加入约束
            Population{2} = Problem.Initialization(); % 不考虑约束
            Population{3} = Problem.Initialization(); % 约束动态变化，e松弛约束，先慢后快
            Fitness{1}    = CalFitness(Population{1}.objs,Population{1}.cons);
            Fitness{2}    = CalFitness(Population{2}.objs);
            Fitness{3}    = CalFitness(Population{3}.objs);
            transfer_state=0;
            generation=1; % 迭代次数
            
            isUPF = 0; % P2是否到达UPF

            % 切比雪夫距离
%             [W,Problem.N] = UniformPoint(Problem.N,Problem.M); 

            % DE_current_to_other_pbest_1
            p = 0.1;
            
            processcon = 0;
            processed = []; % 按约束优先级排序的约束
            index = 1;
            totalcon = size(Population{1}(1,1).con,2); % 约束数目
            firstTotalcon = 0;
            currTotalcon = 0;
            currIndex = 1;
            Pops = [];
            change_threshold = 1e-3;
            change_threshold2   = 1e-2;
            change_rate2        = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            isAdd1 = 0;
            isAdd2 = 0;
            for i = 1:totalcon
                Pops{1,i} = Population{1};
                cons = Population{1}.cons;
                PopCon = max(0,cons(:,i));
                Pops{2,i} = CalFitness(Population{1}.objs,sum(PopCon,2));
                Pops{3,i} = i;
                Pops{4,i} = 0; % 保存只有i约束时最小Pareto支配等级
            end
            Objvalues(1) = sum(sum(Population{1}.objs,1));
            
            % 根据 MSCEA 设置约束边界epsn动态缩小，先慢后快
            [~, nCon]               = size(Population{3}.cons);
            [initialE3, ~]          = max(max(0,Population{3}.cons), [], 1);
            initialE3(initialE3==0) = 1;      
            epsn3                   = initialE3;
            cp = Algorithm.ParameterSet(5);
            
            % 方法二：MTCMO 文献[31]
            consVAR0 = [Population{3}.cons];
            consVAR0(consVAR0<0) = 0;
            VAR0 = max(sum(consVAR0,2));
            if VAR0 == 0
                VAR0 = 1;
            end
            
            %% Optimization
            while Algorithm.NotTerminated(Population{1})
                generation =generation+1;
                
                % 判断约束优先级
                while isAdd2 && processcon <=totalcon
                    AllPops = [];
                    for j = 1:size(Pops,2)
                        AllPops = [AllPops,Pops{1,j}];
                    end
                    [FrontNo,~] = NDSort(AllPops.objs,inf);
                    if processcon == 0
                        ranks = [];
                        for j = 1:size(Pops,2)
                            minF = min(FrontNo((j-1)*Problem.N+1:j*Problem.N));
                            ranks = [ranks,minF];
                        end
                        [~,seq] = sort(ranks);
                        seq = fliplr(seq); % 将数组从左向右翻转
			processed = seq;
                    else
                        processed = [processed,processcon];
                        Pops{4,processcon} = 1;
                        Minindex = min(FrontNo((processcon-1)*Problem.N+1:processcon*Problem.N));
                        for i = 1:size(Pops,2)
                            if i ~= processcon
                                maxindex = max(FrontNo((i-1)*Problem.N+1:i*Problem.N));
                                if maxindex <= Minindex
                                    Pops{4,i} = 1; % 保存只有i约束时最小Pareto支配等级
                                end
                            end
                        end
                    end
                    unpro = 0;
                    for i = 1:size(Pops,2)
                        unpro = unpro+Pops{4,i};
                    end
                    if unpro < totalcon
                        while Pops{4,seq(index)} == 1
                            index = index + 1;
                        end
                        processcon = seq(index);    % 依次取优先级等级最高的约束
                        if ~firstTotalcon
                            firstTotalcon = processcon; % 找到优先级最高的约束
                            currTotalcon = firstTotalcon;
                        end
                    else
                        processcon = totalcon + 1;
                    end
                end                
                
                % 如果P1的解变化小于阈值，则P1进入下一阶段，从约束优先级队列中新增下一个约束
                % 方法一：C3M P1所有目标上的总值变化小于阈值
%                 Objvalues(generation) = sum(sum(abs(Population{1}.objs),1));
%                 isAdd1 = is_stable(Objvalues,generation,Population{1},Problem.N,change_threshold,Problem.M);
                % 方法二：MSCMO P1在每个目标上的平均值变化小于阈值
                change_rate2 = Normalization(Population{1},change_rate2,ceil(Problem.FE/Problem.N));
                isAdd1 = Convertion(change_rate2,generation,100,change_threshold2);
                % 方法三：MOEA-D-DAE P1可行率大于0.9
                 P1Cons = Population{1}.cons;
                 CV    = sum(max(P1Cons(:,currTotalcon),0),2);
                fr    = length(find(CV<=0))/Problem.N;
                isAdd1 = (fr > 0.9) && Convertion(change_rate2,generation,100,change_threshold2);
                
                if currIndex < totalcon && isAdd1
                    currIndex = currIndex + 1;
                    currTotalcon = [currTotalcon,processed(currIndex)];
                    isAdd1 = 0;
                end
                if  Problem.FE/Problem.maxFE >=0.8 && currIndex < totalcon % 如果到最后阶段约束还没加完，直接考虑所有约束
                    currTotalcon = processed;
                    currIndex = totalcon;
                    isAdd1 = 1;
                end
                
                % 如果P3的解都在e-constraint边界内，则P3进入下一阶段，epsn3缩小
                PopCon3   = max(0,Population{3}.cons);
                if sum(sum(PopCon3<=epsn3,2)==nCon) == length(Population{3})
                    % 方法一：MSCEA 文献[46] 模拟退火，减小速度先慢后快
                    epsn3 = ReduceBoundary(initialE3,ceil(Problem.FE/Problem.N),ceil(Problem.maxFE/Problem.N)-1,cp);
                    % 方法二：MTCMO 文献[31]
%                     MTCMO_cp  = (-log(VAR0)-6)/log(1-0.5);
%                     epsn3 = VAR0*(1-(Problem.FE/Problem.maxFE))^MTCMO_cp;                    
                end
                
                if isUPF == 0
                    std_obj(generation,:) = std(Population{2}.objs,[],1);
                    if generation>100
                        if  sum(std(std_obj(generation-100:generation,:),[],1)<0.5) == Problem.M
                            isUPF = 1; % P2 到达UPF，停止进化
                            isAdd2 = 1;
                        end
                    end
                end
                    
                % 早期转移所有子代种群
                if   transfer_state == 0
                    if currIndex == totalcon && isAdd1 && isAdd2
                        transfer_state = 1; % 阶段切换，早期转后期
                    end

                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness{1});
                    valOffspring{1}(1:Problem.N/2) = OperatorGAhalf(Problem,Population{1}(MatingPool1));
%                     valOffspring{1}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{1}, Population{2}, Problem.N/2);
                    valOffspring{1}( 1+Problem.N/2 : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{1},Problem.N/2,Fitness{2},Population{2},p);

                    if isUPF == 0
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
                        valOffspring{2}(1:Problem.N/2) = OperatorGAhalf(Problem,Population{2}(MatingPool2));
%                         valOffspring{2}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{2}, Population{3}, Problem.N/2);
                        valOffspring{2}( 1+Problem.N/2 : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{2},Problem.N/2,Fitness{3},Population{3},p);
                    else
                        valOffspring{2} = [];
%                         MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
%                         valOffspring{2} = OperatorGAhalf(Problem,Population{2}(MatingPool2));
                    end
                    
                    MatingPool3 = TournamentSelection(2,Problem.N,sum(max(0,Population{3}.cons-epsn3),2));
                	valOffspring{3}(1:Problem.N/2) = OperatorGAhalf(Problem,Population{3}(MatingPool3));
%                     valOffspring{3}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{3}, Population{2}, Problem.N/2);
                    valOffspring{3}( 1+Problem.N/2 : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{3},Problem.N/2,Fitness{2},Population{2},p);
                    
                    [Population{1},Fitness{1},~] = EnvironmentalSelection2( [Population{1},valOffspring{1:3}],Problem.N,currTotalcon); % 如果为1则考虑所有约束
                    [Population{2},Fitness{2},~] = EnvironmentalSelection2( [Population{2},valOffspring{1:3}],Problem.N,0); % 不考虑任何约束
                    [Population{3},Fitness{3},~] = EnvironmentalSelection3( [Population{3},valOffspring{1:3}],Problem.N,epsn3);
                  
%                     % 阶段切换，早期转后期
%                     if Problem.FE/Problem.maxFE >=0.2
%                         transfer_state = 1;
%                     end
                  
                % 后期根据辅助任务上的成功率确定转移种群
                else
                    if isUPF == 0
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
                        valOffspring{2}(1:Problem.N/2) = OperatorGAhalf(Problem,Population{2}(MatingPool2));
%                         valOffspring{2}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{2}, Population{1}, Problem.N/2);
                        valOffspring{2}( 1+Problem.N/2 : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{2},Problem.N/2,Fitness{3},Population{3},p);
                    else
                        valOffspring{2} = [];
%                         MatingPool2 = TournamentSelection(2,Problem.N,Fitness{2});
%                         valOffspring{2} = OperatorGAhalf(Problem,Population{2}(MatingPool2));
                    end
                    
                    MatingPool1 = TournamentSelection(2,Problem.N,Fitness{1});
                    valOffspring{1}(1:Problem.N/2) = OperatorGAhalf(Problem,Population{1}(MatingPool1));
%                     valOffspring{1}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{1}, Population{3}, Problem.N/2);
                    valOffspring{1}( 1+Problem.N/2 : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{1},Problem.N/2,Fitness{3},Population{3},p);

                    MatingPool3 = TournamentSelection(2,Problem.N,sum(max(0,Population{3}.cons-epsn3),2));
                	valOffspring{3}(1:Problem.N/2) = OperatorGAhalf(Problem,Population{3}(MatingPool3));
%                     valOffspring{3}( 1+Problem.N/2 : Problem.N) = DE_transfer(Problem, Population{3}, Population{1}, Problem.N/2);
                    valOffspring{3}( 1+Problem.N/2 : Problem.N) = DE_current_to_other_pbest_1(Problem,Population{3},Problem.N/2,Fitness{1},Population{1},p);
                    
                    % 父代中存活到下一代的比率和子代中存活到下一代的比率
%                     [~,~,Next] = EnvironmentalSelection2( [Population{2},valOffspring{2}],Problem.N,2);
%                     succ_rate(1,generation) =  (sum(Next(1:Problem.N))/Problem.N) - (sum(Next(Problem.N+1:end))/Problem.N/2);
                    [~,~,Next] = EnvironmentalSelection2( [Population{1},valOffspring{1}],Problem.N,-1);
                    succ_rate(3,generation) =  (sum(Next(1:Problem.N))/Problem.N) - (sum(Next(Problem.N+1:end))/Problem.N/2);
                    [~,~,Next] = EnvironmentalSelection3( [Population{3},valOffspring{3}],Problem.N,epsn3);
                    succ_rate(1,generation) =  (sum(Next(1:Problem.N))/Problem.N) - (sum(Next(Problem.N+1:end))/Problem.N/2);
                    
                    % P2到达UPF之后，P2停止进化，判断P2是否要参与P2和P3的后续环境选择
                    
                    if   succ_rate(1,generation) >0
                        rand_number = randperm(Problem.N);
                        [Population{1},Fitness{1},~] = EnvironmentalSelection2( [Population{1},valOffspring{1},valOffspring{2},Population{3}(rand_number(1:Problem.N/2))],Problem.N,-1);
                    else
                        [Population{1},Fitness{1},~] = EnvironmentalSelection2( [Population{1},valOffspring{1},valOffspring{2},valOffspring{3}],Problem.N,-1);
                    end
                    if   succ_rate(3,generation) >0
                        rand_number = randperm(Problem.N);
                        [Population{3},Fitness{3},~] = EnvironmentalSelection3( [Population{3},valOffspring{3},valOffspring{2},Population{1}(rand_number(1:Problem.N/2))],Problem.N,epsn3);
                    else
                        [Population{3},Fitness{3},~] = EnvironmentalSelection3( [Population{3},valOffspring{3},valOffspring{2},valOffspring{1}],Problem.N,epsn3);
                    end
                end
            end
        end
    end
end

% C3M
function result = is_stable(Objvalues,gen,Population,N,change_threshold,M)
    result = 0;
    [FrontNo,~]=NDSort(Population.objs,size(Population.objs,1));
    NC=size(find(FrontNo==1),2);
    max_change = abs(Objvalues(gen)-Objvalues(gen-1));
    if NC == N
        change_threshold = change_threshold * abs(((Objvalues(gen) / N))/(M))*10^(M-2);
        if max_change <= change_threshold
            result = 1;
        end
    end
end
% MSCMO
function outcome = Convertion(Average,G,last_gen,change_threshold)
    outcome = 0;
    if G > last_gen
        k = max(abs(Average(G,:)-Average(G-last_gen,:)));
        if k<=change_threshold
            outcome = 1;
        end
    end
end