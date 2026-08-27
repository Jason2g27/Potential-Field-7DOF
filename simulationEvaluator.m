function [ind1, ind2] = simulationEvaluator(humantorque,inducedtorque,impairment,performance, debug)
    % check direction -> whether the self-interactive characteristics are
    % fully used
    mild = humantorque(find(impairment==2),:).*inducedtorque(find(impairment==2),:);
    severe = humantorque(find(impairment==1),:).*inducedtorque(find(impairment==1),:);
    num_tot = length(humantorque(1,:));
    num_mild = length(find(mild<=0));
    num_severe = length(find(severe>=0));
    num_mild_true = length(find(mild<0));
    num_severe_true = length(find(severe>0));
    percentage_mild = num_mild/num_tot;
    percentage_severe = num_severe/num_tot;
    percentage_mild_true = num_mild_true/num_tot;
    percentage_severe_true = num_severe_true/num_tot;
    % check transmission ratio -> overal difficulty
    trans = abs(sum(inducedtorque(find(impairment==2),find(mild<0)))) / abs(sum(inducedtorque(find(impairment==1),find(severe>0))));

    if debug
        fprintf("Portion for mildly impaired joints: %.2f\n", percentage_mild);
        fprintf("Portion for severely impaired joints: %.2f\n", percentage_severe);
        fprintf("True portion for mildly impaired joints: %.2f\n", percentage_mild_true);
        fprintf("True portion for severely impaired joints: %.2f\n", percentage_severe_true);
        fprintf("Transmission ratio: %.2f\n", trans);
    end
    if percentage_mild >= performance.overalpct && percentage_severe >= performance.overalpct && percentage_mild_true >= performance.truepct && percentage_severe_true >= performance.truepct
        ind1 = 1;
    else
        ind1 = 0;
    end
    if trans >= performance.ratio
        ind2 = 1;
    else
        ind2 = 0;
    end

end