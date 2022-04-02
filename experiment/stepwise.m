function [lower, higher] = stepwise(f, lower, higher, step)
    node = lower;
    
    while node < higher
        if f(node)*f(node+step) < 0
            lower = node;
            higher = node + step;
            break
        end
        
        node = node + step;
    end
end