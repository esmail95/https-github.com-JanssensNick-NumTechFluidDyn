function stop = stoppingCriterion(u_corr,v_corr,p_corr,tol)
%STOPPINGCRITERION Summary of this function goes here
%   Detailed explanation goes here
    stop = false;
    u_corr = norm(u_corr);
    v_corr = norm(v_corr);
    p_corr = norm(p_corr);
    Unorm = max(u_corr,v_corr);
    if max(Unorm,p_corr) < tol
        stop = true;
    end
end

