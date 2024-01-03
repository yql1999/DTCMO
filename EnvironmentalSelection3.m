function [Population,Fitness, Next] = EnvironmentalSelection3(Population,N,epsn)
% The environmental selection1 of MSCEA
   
    [~, nCon] = size(Population.cons);
    PopCon    = max(0,Population.cons);
    if sum(sum(PopCon<=epsn, 2)==nCon) > N       % 如果父代和子代中满足松弛约束的解>N，在这些满足松弛约束的解中进行适应度选择
        tmp        = sum(PopCon<=epsn, 2)==nCon;
        Population = Population(1:end, tmp);
        CV         = sum(max(0,Population.cons),2);
        Fitness    = CalFitness(Population.objs,CV);
        Next = Fitness < 1;
        if sum(Next) < N
            [~,Rank] = sort(Fitness);
            Next(Rank(1:N)) = true;
        elseif sum(Next) > N
            Del  = Truncation(Population(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
    else       % 如果父代和子代中满足松弛约束的解<N，直接根据CV选择
        CV         = sum(max(0,Population.cons),2);
        [~, rank]  = sort(CV);
        Next       = rank(1:N);
       Fitness    = CalFitness(Population.objs,CV);
    end   
    Population = Population(Next);
    Fitness = Fitness(Next);
end

function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end