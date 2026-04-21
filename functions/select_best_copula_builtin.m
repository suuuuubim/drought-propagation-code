function [best_family, best_params, best_aic] = select_best_copula_builtin(uv_data, candidates)
    best_aic = inf; best_family = ''; best_params = {}; 
    for i = 1:length(candidates)
        family = candidates{i};
        try
            params_cell = {}; num_params = 1;
            if strcmpi(family, 't')
                [rho_matrix, nu] = copulafit(family, uv_data);
                params_cell = {rho_matrix, nu}; num_params = 2;
            else
                param = copulafit(family, uv_data);
                params_cell = {param};
            end
            pdf_vals = copulapdf(family, uv_data, params_cell{:});
            logL = sum(log(pdf_vals + 1e-10)); 
            aic = -2*logL + 2*num_params;
            if ~isnan(aic) && aic < best_aic
                best_aic = aic; best_family = family; best_params = params_cell;
            end
        catch
        end
    end
end