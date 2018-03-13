function [goodness] = sd(window1, window2)
    window1 = window1(:);
    window2 = window2(:);
    goodness = sum((window1 - window2).^2);
end

