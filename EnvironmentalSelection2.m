function [Population,Fitness,Next] = EnvironmentalSelection2(Population,N,currTotalcon)
% The environmental selection of SPEA2
% -1所有约束；0无约束；其他某一些约束
    %% Calculate the fitness of each solution
    if currTotalcon==-1
        CV = sum(max(0,Population.cons),2);
        Fitness = CalFitness(Population.objs,CV);
    elseif currTotalcon==0
        Fitness = CalFitness(Population.objs);
    else
        Cons = Population.cons;
        CV = Cons(:,currTotalcon);
        CV    = sum(max(0,CV),2);
        Fitness = CalFitness(Population.objs,CV);
    end

    %% Environmental selection
    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = Truncation(Population(Next).objs,sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end

    % Population for next generation
    Population = Population(Next);
    Fitness    = Fitness(Next);
    % Sort the population
    [Fitness,rank] = sort(Fitness);
    Population = Population(rank);
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