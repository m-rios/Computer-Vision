function [goodness] = ad(window1, window2)
    window1 = window1(:);
    window2 = window2(:);
    goodness = sum(abs(window1 - window2));
end

