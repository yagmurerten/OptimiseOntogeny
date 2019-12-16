function changed_str=get_trait(focal_str,other_str,trait)
    trait_index=0;
    if (trait==string('prob_a'))
        trait_index=1;
    elseif (trait==string('prob_d'))
        trait_index=2;
    elseif (trait==string('tel'))
        trait_index=3;
    elseif (trait==string('onco_thr'))
        trait_index=4;
    elseif (trait==string('onco_prob'))
        trait_index=5;
    elseif (trait==string('dif_l'))
        trait_index=6;
    elseif (trait==string('div_prop'))
        trait_index=7;
    end
    changed_str=focal_str;
    changed_str(trait_index)=other_str(trait_index);
end