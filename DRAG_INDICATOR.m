function drag_indie = DRAG_INDICATOR(x, R)
    F_CW = 15000;
    s1 = x(1:21);
    s2 = x(22:42);
    sl = x(43:end);
    Radii_cuntribution = F_CW.*sl./R;
    Stag_cuntribution_pos = F_CW*2*(s1-s2)./sl;
    Stag_cuntribution_neg = -Stag_cuntribution_pos;
    drag_pos = Stag_cuntribution_pos + Radii_cuntribution;
    drag_neg = Stag_cuntribution_neg + Radii_cuntribution;
    drag_pos= abs(drag_pos);
    drag_neg = abs(drag_neg);
    VIOLATION_pos = find(drag_pos > 900 | drag_pos < 70);
   	VIOLATION_neg = find(drag_neg > 900 | drag_neg < 70);
    if ~isempty(VIOLATION_pos) | ~isempty(VIOLATION_neg)
        drag_indie = 1;
    else
        drag_indie = 0;
    end
end