%% Function to estimate regression statistics given linear and quadratic
%models fit to the same data.

% Edited July 10th, 2026, RJLIII


function [Lin_NRMSE,Lin_p,Lin_VarExplained,Quad_NRMSE,Quad_p,Quad_VarExplained] ...
    = Regression_Metrics(X, Y, Lin_model, Quad_model, model_type)


if model_type == "polyfit"
    try
        n = length(X);
        
        %% Normalized Root Mean Square Error Estimate
        %fprintf("If you get an error evaluating the model, check your model type.\n")
        
        Lin_pred = polyval(Lin_model, X); % Evaluate the linear model
        Quad_pred = polyval(Quad_model, X); % Evaluate the quadratic model
        
        Lin_resid = Y - Lin_pred;
        Quad_resid = Y - Quad_pred;
        
        Lin_RMSE = sqrt(mean(Lin_resid.^2));
        Quad_RMSE = sqrt(mean(Quad_resid.^2));
        
        Lin_NRMSE = Lin_RMSE / std(Y);
        Quad_NRMSE = Quad_RMSE / std(Y);
        
        %% Anova Test
        % Define the degrees of freedom
        k_quad = 3;
        % Difference in degrees of freedom
        df_err_quad = n - k_quad;
        df_x2 = 1; 
        
        % sum of squares (SS)
        SS_total = sum((Y - mean(Y)).^2);
        SS_err_lin = sum(Lin_resid.^2);
        SS_err_quad = sum(Quad_resid.^2);
        SS_reg_lin = SS_total - SS_err_lin;
        SS_reg_quad = SS_total - SS_err_quad;
        SS_x2 = SS_err_lin - SS_err_quad; 
        % Defines the variance explained by adding the quadratic term
        
        % Mean squares (MS)
        MS_x_lin = SS_reg_lin /1;
        MS_x2 = SS_x2 / df_x2;
        MS_err_quad = SS_err_quad / df_err_quad;
        
        % F-stats
        F_x = MS_x_lin / MS_err_quad;
        F_x2 = MS_x2 / MS_err_quad;
        
        % p-values
        Lin_p = 1 - fcdf(F_x, 1, df_err_quad);
        Quad_p = 1 - fcdf(F_x2, 1, df_err_quad);
        
        Lin_VarExplained = 100*(SS_reg_lin / (SS_reg_lin + SS_reg_quad));
        Quad_VarExplained = 100*(SS_reg_quad / (SS_reg_lin + SS_reg_quad));
        fprintf("Regression stats done with polyfit model.\n")

    catch err
        fprintf("Error occurred: check your model type.")
    end
elseif model_type == "fitlm"
   try
        %fprintf("If you get an error evaluating the model, check your model type.\n")
        
        Lin_NRMSE = (sqrt(mean(Lin_model.Residuals.Raw.^2)) / std(Y));
        Quad_NRMSE = (sqrt(mean(Quad_model.Residuals.Raw.^2)) / std(Y));
        Lin_p = anova(Lin_model).pValue(1);
        Quad_p = Quad_model.ModelFitVsNullModel.Pvalue;
        Lin_VarExplained = 100*(anova(Quad_model).SumSq(1) / (anova(Quad_model).SumSq(1) + anova(Quad_model).SumSq(2)));
        Quad_VarExplained = 100*(anova(Quad_model).SumSq(2) / (anova(Quad_model).SumSq(1) + anova(Quad_model).SumSq(2)));
        fprintf("Regression stats done with fitlm model.\n")
   catch err
        fprintf("Error occurred: check your model type.")
   end
else 
    fprintf("Model type not recognized. Acceptable inputs: 'fitlm' or 'polyfit'\n")
end