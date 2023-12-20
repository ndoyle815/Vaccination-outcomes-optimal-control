% basic cost function to evaluate

% weights: (w1,w2) determine how we weigh disease burden and lockdown
% stringency component
% para: ODE parameter set for main simulation
% NB: thresholds defined as para.Imin, para.Imax in para
% Hospital capacity defined as para.Hmax in para

function F = CostFunction(weights, para, burden, stringency, peak_hospital, apply_discounting)

MAXBURDEN = 33600;%33655;    % informed by Rapid Control strategy under worst-case scenario (no vaccine)
MAXSTRINGENCY = 275;%270;  % informed by Suppression strategy under worst-case scenario (no vaccine)

if apply_discounting == 1
    d = 0.035;
    discounting = 1./(1 + d).^([para.t0:para.maxtime]./365)';

    Burden = round(sum(discounting.*burden))/MAXBURDEN;
    Stringency = sum(discounting.*stringency)/MAXSTRINGENCY;
else
    Burden = round(sum(burden))/MAXBURDEN;
    Stringency = sum(stringency)/MAXSTRINGENCY;
end

F = weights(1)*Burden + (1-weights(1))*Stringency + exp(weights(2)*(peak_hospital-para.Hmax));