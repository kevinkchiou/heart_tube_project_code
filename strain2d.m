function s2d = strain2d(E_c, E_s, s0, E_t)
    model = mphload('cell_on_gel');
    model.param.set('E_c', strcat(num2str(E_c), '[kPa]'));
    model.param.set('E_s', strcat(num2str(E_s), '[kPa]'));
    model.param.set('s0', b(1));
    
    s2d = zeros(size(E_t));
    for i=1:length(E_t)
        E = E_t(i);
        model.param.set('E_t', strcat(num2str(E), '[kPa]'));
        model.study('std1').run;
        s2d(i) = abs(mphmean(model, 'strain2d', 'surface', 'selection', 5));
    end
end