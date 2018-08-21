function xxx = clustergrid(lo, mid, hi, n1, n2, c1, c2)
%@assert lo < mid < hi "lo, mid, hi must be ordered in that way"
    if lo >= 0  % all are positive
        g2 = linspace(mid^c2, hi^c2, n2).^(1/c2);
        g1 = mid - linspace(mid^c1, lo^c1, n1+1).^(1/c1);
        xxx = [g1(1:end-1), g2];
    else  % lo is negative
        %shift all lo, mid, hi up by abs(lo) and form grid when all are positive
        % then shift back down by abs(lo) to re-align boundaries
        alo = abs(lo);
        xxx = clustergrid(0, mid+alo, hi+alo, n1, n2, c1, c2) - alo;
    end
end