function [res, err] = bisection(f, lower, higher, err)
    while f(lower)*f(higher)<0
        mid = (lower+higher)/2;
        
        if f(lower)*f(mid) < 0
            higher = mid;
        elseif f(higher)*f(mid) < 0
            lower = mid;
        end
        
        
        if  abs(f(mid)) < err
            break;
        end
    end
    
    err = abs(f(mid));
    res = mid;
end