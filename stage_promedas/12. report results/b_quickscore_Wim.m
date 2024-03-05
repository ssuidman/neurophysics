function [pfplus, P_joint, posterior, dt, dt_array] = b_quickscore_Wim(previn,pfmin,pfminneg)
    
    if(~isempty(pfmin))
        [m,n] = size(pfmin);     % m number of pos tests. n number of diags.
        if(n ~= length(previn))
            error('inconsistent n');
        end
    else
        m=0;                     % 0  pos tests.
        n = length(previn) ;
    end
    
    pfplus=zeros(n+1,1);      % will be used for  P(F+,F-), P(F+,F-|d_1=1),... P(F+,F-|d_n=1)  (note F- will be absorbed below)
    prev = repmat(previn, n+1,1);  % copies of prevalences, but with (i+1,i) entry 1 to condition on d_i=1,  etc.
    for i = 1:n; prev(i+1,i)=1; end;

    % absorb negative findings in prevminneg  (so these are p(F-|d_1=1)p(d_1=1), ... p(F-|d_n=1)p(d_n=1) )
    % needed in heckerman eqn  11 (compare with heckerman 10, there is a
    % constant difference which is incorporated in the block below)
    if(~isempty(pfminneg))
        prevminneg = prev.*pfminneg;
    else
        prevminneg = prev;
    end
    
    t1 = tic;
    % pfplus_matrix=zeros(2^m,n+1,1);
    dt_array = zeros(m+1,1);
    for i=0:(2^m-1)
        v = de2bi(i,m);
        myset = find(v==1);
        if(isempty(myset))
            pfplus = pfplus+exp(sum(log(1e-50+prevminneg + (1-prev)),2));   % I use the exp(sum(log( variant of prod. Small number is added to prevent log(0) - maybe there is a better solution
    %         pfplus_matrix(i+1,:) = prod(1e-50+prevminneg + (1-prev),2);
%             pfplus = pfplus+prod(1e-50+prevminneg + (1-prev),2);   % I use the exp(sum(log( variant of prod. Small number is added to prevent log(0) - maybe there is a better solution
        elseif(length(myset)==1)
            pfplus = pfplus + ((-1)^length(myset))*exp(sum(log(1e-50+pfmin(myset,:).*prevminneg + (1-prev)),2));
    %         pfplus_matrix(i+1,:) = ((-1)^length(myset))*prod(1e-50+ (pfmin(myset,:).*prevminneg + (1-prev)),2);
%             pfplus = pfplus + ((-1)^length(myset))*prod(1e-50+ (pfmin(myset,:).*prevminneg + (1-prev)),2);
            dt = toc(t1);
            dt_array(myset(1),1) = dt; 
    %         dt_array = [dt_array,dt];
            fprintf('%d:\t[ ',myset(1));
            fprintf('%d ', v);
            fprintf(']\t--> %5f\n',dt);
        else
            pfplus = pfplus + ((-1)^length(myset))*exp(sum(log(1e-50+ (prod(pfmin(myset,:)).*prevminneg + (1-prev))),2));
    %         pfplus_matrix(i+1,:) = ((-1)^length(myset))*prod(1e-50+ (prod(pfmin(myset,:)).*prevminneg + (1-prev)),2);
%             pfplus = pfplus + ((-1)^length(myset))*prod(1e-50+ (prod(pfmin(myset,:)).*prevminneg + (1-prev)),2);
        end
    end
    dt = toc(t1);
    dt_array(end,1) = dt; 
    dt_array = dt_array(2:end); 
    P_joint = pfplus(2:end)' .* previn; 
    posterior = P_joint / pfplus(1); 
    
    fprintf('Running time: %f\n',dt(end))
end
