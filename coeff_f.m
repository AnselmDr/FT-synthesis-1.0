%% PDE-Coefficient f
function [f] = coeff_f(location,state,obj)
    % COEFF_F Coefficient f of the PDE, arbitrary term, in this case source and first
    % spacial derivative term (convective term)
    
    % Rename for easier readability
    % state.u: line_1 to line_numSpeciesInlet represents mole flow
    % state.u: line_numSpeciesInlet+1 represents Temperature, if applicable
    % state.ux: their respective first derivative in x-direction

    obj.counter = obj.counter +1;
    
    % nCurr contains only the lines regarding mole flow
    nCurr = state.u(1:obj.numSpeciesInlet,:);
    dnCurrdx = state.ux(1:obj.numSpeciesInlet,:);
    
    % TCurrRct contains only the last line regarding Temperature
    if obj.TRctSwitch && ~obj.TWallSwitch
        TCurrRct = state.u(obj.numSpeciesInlet+1,:);
        dTCurrRctdx = state.ux(obj.numSpeciesInlet+1,:);
        TCurrWall = obj.TConst*ones(1,length(location.x));
    elseif obj.TRctSwitch && obj.TWallSwitch
        TCurrRct = state.u(obj.numSpeciesInlet+1,:);
        dTCurrRctdx = state.ux(obj.numSpeciesInlet+1,:);
        TCurrWall = state.u(obj.numSpeciesInlet+2,:);
        TCool = obj.TCool*ones(1,length(location.x));
    else
        TCurrRct = obj.TConst*ones(1,length(location.x));
        fTCurr = obj.TConst*ones(1,length(location.x));
    end
    
    % nSum is current total mole flow
    nSum = sum(nCurr(1:obj.numSpeciesInlet,:));
    
    % Set pCurr to pInlet since isobaric is assumed
    pCurr = obj.pInlet;
    
    % numX: number of values in x direction
    numX = length(location.x);
    
    % Calculate parameters dependent on current values
    cPRct = obj.specHeatCap(TCurrRct,nCurr,nSum);
    V = obj.volumeFlow(pCurr,TCurrRct,nCurr);%obj.VInlet.*ones(1,numX)
    u = V./(obj.volumeRct/obj.lengthRct);%.*ones(1,numX);
    cWall = obj.solidPhase.specHeatCapS('T',TCurrRct);
    lambdaWall = obj.solidPhase.thermCondS('T',TCurrRct);
    rhoRct = obj.densityGas(pCurr,TCurrRct,nCurr,nSum);
    
    if obj.kModeSwitch == 1
        % Set constant k_i
        kRct = obj.kRct.*ones(1,numX);
        kCool = obj.kCool.*ones(1,numX);
    elseif obj.kModeSwitch == 2
        % Calculating parameters needed for Heat Law
        lambdaRct = obj.thermCondGas(pCurr,TCurrRct,nCurr,nSum);
        rho_g = densityGas(obj,pCurr,TCurrRct,nCurr,nSum);
        eta = dynViscGas(obj,pCurr,TCurrRct,nCurr,nSum);

        % Calculating dimensionless numbers
        Re = u*obj.dRct*rho_g/eta;
        Pr = eta.*cPRct./lambdaRct;

        % Calculation of heat transfer per unit length / W/m
        Nu = lawNusseltPipe(Re, Pr, obj.dRct, location.x);
        kRct = Nu*lambdaRct/obj.dRct;
        
        % No heat law for cooling side implemented yet
        kCool = obj.kCool.*ones(1,numX);
    elseif obj.kModeSwitch == 3
        % Calculate k_i from alpha_i
        kRct = 1./(1/obj.alphaRct + obj.thicknessWallmeanRct*obj.dRct./(lambdaWall*obj.dLogMeanRct));
        kCool = 1./(1/obj.alphaCool + obj.thicknessWallmeanCool*obj.dCool./(lambdaWall*obj.dLogMeanCool));
    elseif obj.kModeSwitch == 4
        kRct = 1./(1/obj.alphaRct + obj.thicknessWall./(2*lambdaWall));
        kCool = 1./(1/obj.alphaCool + obj.thicknessWall./(2*lambdaWall));
    end
    
    % Rate Law
    partial = nCurr./nSum * pCurr;
    Rate = obj.rateFactor*rateLawKwackViaKirsch(obj.R,TCurrRct,partial,obj.numSpeciesInlet);
    
    if ~isreal(cPRct)
        disp('cPRct is complex!')
    end
    if ~isreal(V)
        disp('V is complex!')
    end
    if ~isreal(partial)
        disp('partial is complex!')
    end
    if ~isreal(Rate)
        disp('Rate is complex!')
    end    
    
    %% PDE Coefficients
    % assuming m independent variables (mole flow species 1 ...,
    % Temperatures zone 1...)
    % and k locations, an m x k matrix f linewise from 1 x k line vectors
    % m and k can be derived from size of state.u/state.ux
    
    %% Mass balance
    fnCurr = zeros(obj.numSpeciesInlet,numX);
    for i = 1:obj.numSpeciesInlet % number of species / rows
        for j = 1:numX % number of locations in x direction / columns
                fnCurr(i,j) = -u(j) * dnCurrdx(i,j) ...                               % Convective term
                        + V(j).*obj.rhoCat .* Rate(i,j) ...                               % Source term - Reaction
                        ;
        end
    end
                
                
    %% Enthalpy balance Rct
    if obj.TRctSwitch
        fTCurr = zeros(1,numX);
        for j = 1:numX % number of locations in x direction / columns || .*molarWeightRct(j)
                fTCurr(j) = -u(j) * dTCurrRctdx(j) ...                                                                         % Convective term 
                    + (obj.rhoCat) ./ (cPRct(j).*rhoRct(j)) * obj.DeltaRctH * Rate(2,j) ...                                    % Source term - Reaction
                    - (kRct(j).*obj.URct)./(obj.areaCrossRct.*cPRct(j).*rhoRct(j)) * (TCurrRct(j) - TCurrWall(j)) ...   % Source term - Heat transfer wall
                    ;
        end
    end
    
    %% Enthalpy balance Wall
    if obj.TWallSwitch
        fTWallCurr = zeros(1,numX);
        for j = 1:numX
            fTWallCurr(j) = (kRct(j).*obj.URct) ./ (obj.areaCrossWall.*obj.rhoWall.*cWall(j)) * (TCurrRct(j) - TCurrWall(j)) ...
                - (kCool(j).*obj.UCool)./(obj.areaCrossWall.*obj.rhoWall.*cWall(j)) * (TCurrWall(j) - TCool(j)) ...
                ;
        end
    end
    
    %% Concatenate coeff_f array
    if obj.TRctSwitch && ~obj.TWallSwitch
        f = [fnCurr; fTCurr];
    elseif obj.TRctSwitch && obj.TRctSwitch
        f = [fnCurr; fTCurr; fTWallCurr];
    else
        f = fnCurr;
    end
    
    if numX > 2
        % sort location.x and get indices
        [obj.locXSort,obj.idxSort] = sort(location.x);
        % sort parameters with the help of indices
        obj.kFinal = [kRct(1,obj.idxSort); kCool(1,obj.idxSort)];
        obj.uFinal = u(1,obj.idxSort);
        obj.VFinal = V(1,obj.idxSort);
        obj.RateFinal = Rate(1:end,obj.idxSort);
        obj.location = location;
        
%         obj.kFinal = [kRct;kCool];
%         obj.uFinal = u;
%         obj.VFinal = V;
%         obj.RateFinal = Rate;
%         obj.location = location;
    end
end