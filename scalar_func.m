%% scalarizing functions for decomposition methods
function max_fun = scalar_func(y_obj,idealpoint,namda)
%% // Tchebycheff approach
    max_fun = -1.0e+30;
    for n=1:2
        diff = abs(y_obj(n) - idealpoint(n) );
        if namda(n)==0
            feval = 0.00001*diff;
        else
            feval = diff*namda(n);
        end
        if feval>max_fun
            max_fun = feval;
        end
    end
end
