%% PDE-Coefficient c
function [c] = coeff_c(location,state,obj)
    % COEFF_C sets diffusive coefficients of the PDE (the second spacial
    % derivative)
    
    % Rename for easier readability
    % state.u: line_1 to line_numSpeciesInlet represents mole flow
    % state.u: line_numSpeciesInlet+1 represents Temperature of reaction,
    %          if aplicable
    % state.ux: their respective first derivative in x-direction

    % nCurr contains only the lines regarding mole flow
    nCurr = state.u(1:obj.numSpeciesInlet,:);
    % TCurrRct contains only the last line regarding Temperature of
    % reaction
    if obj.TRctSwitch
        TCurrRct = state.u(obj.numSpeciesInlet+1,:);
    else
        TCurrRct = obj.TConst*ones(1,length(location.x));
    end
    
    
    % nSum is current total mole flow
    nSum = sum(nCurr(1:obj.numSpeciesInlet,:));
    % Set pCurr to pInlet since isobaric is assumed
    pCurr = obj.pInlet;
    
    % Preallocate c
    % numX: number of values in x direction
    numX = length(location.x);
    c = zeros(obj.numEquations,numX);
    
    % Set Diffusion coefficient to constant value since method is not
    % implemented yet
    D = obj.diff(pCurr,TCurrRct,nCurr,nSum);
    
    % Workaround if numSpeciesInlet is not equal numSpeciesGasPhase
    if obj.numSpeciesInlet ~= obj.numSpeciesGasPhase
        D(obj.numSpeciesGasPhase+1:obj.numSpeciesInlet,:) = ones(obj.numSpeciesInlet-obj.numSpeciesGasPhase,1) * D(obj.numSpeciesGasPhase,:);
    end
    
    if ~isreal(D)
    disp('D is complex!')
    end
    % Calculate parameters dependent on current values
    cPRct = obj.specHeatCap(TCurrRct,nCurr,nSum);
    lambdaRct = obj.thermCondGas(pCurr,TCurrRct,nCurr,nSum);
    lambdaWall = obj.solidPhase.thermCondS('T',TCurrRct);
    cWall = obj.solidPhase.specHeatCapS('T',TCurrRct);
    rhoRct = obj.densityGas(pCurr,TCurrRct,nCurr,nSum);

    %% PDE Coefficients
    % assuming m independent variables (mole flow species 1 ...,
    % Temperatures zone 1...)
    % and k locations, an m x k matrix f linewise from 1 x k line vectors
    % m and k can be derived from size of state.u/state.ux
    
    % Mass balance
    c(1:obj.numSpeciesInlet,:) = obj.DFactor * D;
    
    % Enthalpy balance
    if obj.TRctSwitch
        c(obj.numSpeciesInlet+1,:) = obj.aRctFactor * ((lambdaRct)./(cPRct.*rhoRct));
    end
    
    if obj.TWallSwitch
        c(obj.numSpeciesInlet+2,:) = obj.aWallFactor * lambdaWall./(cWall.*obj.rhoWall);
    end
    
    % saves coefficient c as class property for each step an array for c is
    % created
    if numX > 2
        % sort location.x and get indices
        [obj.locXSort,obj.idxSort] = sort(location.x);
        % sort parameters with the help of indices
        obj.DFinal = D(1:end,obj.idxSort);
        obj.lambdaRctFinal = lambdaRct(1,obj.idxSort);
        obj.cPRctFinal = cPRct(1,obj.idxSort);
        obj.rhoRctFinal = rhoRct(1,obj.idxSort);
        if obj.TRctSwitch
            obj.aRctFinal = c(obj.numSpeciesInlet+1,obj.idxSort);
        end
    end

end