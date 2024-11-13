function [Sw2Plot,SwSD2Plot] = GetSwTime4Plots(TCond)
        %CycleTime,StanceTime

        mSwTime = mean(TCond.CycleTime - TCond.StanceTime);
        sdSwTime = std(TCond.CycleTime - TCond.StanceTime);
        mCycleTime = mean(TCond.CycleTime);% SwTime,CycleTimeSwSt,CycleTimeStSw

        Sw2Plot = round(100*mSwTime/mCycleTime);
        SwSD2Plot = round(100*sdSwTime/mCycleTime);
