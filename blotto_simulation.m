function strategy_dead = blotto_simulation
% blotto simultation
% pawns allocated to n+1 must be greater than pawns at n

clc; 
strategy_live = cell(1); % has pawns left -at the start all strategies are live
strategy_live{1}.pawns_alloc = [zeros(1,10),100]; %the last is a dummy
strategy_live{1}.pawns_left = 100;

global num_dead strategy_dead
strategy_dead = cell(10000); %no pawns left
num_dead = 0;


for castle = 10:-1:1
    
    num_live = length(strategy_live);    
     %the children inheriting from the current live set
    strategy_live_new = {};
    new_lv_ctr = 0;
    
    for lv = 1 : num_live
        % take a live strategy and spawn
        curr_strategy = strategy_live{lv};
        
        % starting the counter with 
        % this value enforces the constraint on the allocation while
        % making sure that there are enough pawns left
        start_value = min( curr_strategy.pawns_alloc(castle+1), ...
                    curr_strategy.pawns_left);
        if start_value == 0
            error('fuckup in making em dead');
        end
        for pawns = start_value :-1: 1
            % check if this allocation strategy is viable, i.e. all 
            % pawns can be used
            if pawns * castle < curr_strategy.pawns_left 
                % this allocation strategy and by extension all 
                % lower allocations will not use all pawns
                % fuhgedabouddit and move on to the next live strategy
                break;
            end    % if pawns         

            %inherit your parent's allocation
            new_strategy.pawns_alloc = curr_strategy.pawns_alloc;
            new_strategy.pawns_alloc(castle) = pawns;
            new_strategy.pawns_left = curr_strategy.pawns_left - pawns;
            if new_strategy.pawns_left == 0
                num_dead = num_dead + 1;
                strategy_dead{num_dead} = new_strategy; 
                strategy_dead{num_dead}.wins = 0;
                process_dead_list;
            elseif new_strategy.pawns_left > 0
                new_lv_ctr = new_lv_ctr + 1;
                strategy_live_new{new_lv_ctr} = new_strategy;  
            elseif new_strategy.pawns_left < 0
                % not enough for this allocation fuhgeddaboudit
                error('fuckup in pawn counting'); 
            end %if new_strategy
        end %for pawns
    end %for lv
    % done processing all live ones for this castle
    strategy_live = strategy_live_new;
    display( [' num_live=',num2str(num_live), ...
              ' num_dead=',num2str(num_dead), ...
              ' castle=',num2str(castle)]);
    save 'blotto_sim.mat' strategy_dead num_dead num_live castle
end % for castle 
return

function process_dead_list
    % uses the strategy_dead_array to update win-loss count
    global num_dead strategy_dead    
    if num_dead > 1
        for dd = 1 : num_dead-1;
            result = compare( strategy_dead{dd}, strategy_dead{num_dead} );
            if 1 == result %strategy dd
                strategy_dead{dd}.wins = strategy_dead{dd}.wins+1;
            elseif 2 == result %strat num_dead
                strategy_dead{num_dead}.wins = strategy_dead{num_dead}.wins+1;
            elseif 3 == result %tie
                strategy_dead{dd}.wins = strategy_dead{dd}.wins+0.5;
                strategy_dead{num_dead}.wins = strategy_dead{num_dead}.wins+0.5;
            else
                error
            end % if results
        end % for dd            
    end % if num_dead
    return ;

function result = compare(strategy_1, strategy_2)
    castles_1 = find( strategy_1.pawns_alloc(1:10) > strategy_2.pawns_alloc(1:10));
    castles_2 = find( strategy_1.pawns_alloc(1:10) < strategy_2.pawns_alloc(1:10));
    castles_tied = find( strategy_1.pawns_alloc(1:10) == strategy_2.pawns_alloc(1:10));
    
    if sum(castles_1) + sum(castles_2) + sum(castles_tied) ~= 55
        error ('fuckup in compare');
    end
    
    if sum(castles_1) > sum(castles_2) 
        result = 1;
    elseif sum(castles_1) > sum(castles_2) 
        result = 2;
    else
        result = 3;
    end
    return;
    
        
    

    