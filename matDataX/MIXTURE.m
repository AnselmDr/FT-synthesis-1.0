classdef MIXTURE < handle
    %MIXTURE can be used to describe thermodynamic properties of mixtures
    %   Instances of this class can be used to describe thermodynamic
    %   properties of mixtures. It comprises methods that calculate
    %   thermodynamic properties of mixtures based on correlations from the
    %   VDI Heat Atlas. Use curly brackets to list material names. This
    %   equals a cell array.
    %   Provide arguments as name value pairs.
    %   
    %   For molarWeight: Provide argument 'y' as column vector. Provide as
    %   many columns for 'y' as you wish. Result will be a line vector
    
    properties
        
        names               % cell line vector: designations of all materials
        N                   % integer: number of all species
        materials           % array of instances of xMATERIAL: encoding 
                            % material properties of the components of the mixture
        constMIXTURESwitch  % boolean: true, if MIXTURE properties shall be calculated with constant pressure and temperature
        pConst              % pressure for pressure-and-temperature-constant property in bar
        TConst              % temperature for pressure-and-temperature-constant property in K
        diffBinArray        % array: array for binary diffusion coefficients 
                            % not dependent on T or y 
                    
    end
    
    methods
        
        %% Constructor
        function obj = MIXTURE(varargin)
            %MIXTURE Construct an instance of this class
            %   Input arguments
            %   - one dimensional cell vector with designations of all materials in
            %     the mixture, e.g. {'H2', 'CO', ..} or semicolon instead of
            %     comma
            %   Output arguments
            %   - instance of class
            
            % Create input parser
            parser = inputParser;
            
            % Define input parser elements
            addParameter (parser,'names',[]);                       % array of strings: species names
            addParameter (parser,'constMIXTURESwitch',false);   % boolean: true, if MIXTURE properties shall be calculated with constant pressure and temperature
            addParameter (parser,'pConst',1.013);                   % scalar or line vector of double: for pressure-and-temperature-constant property
            addParameter (parser,'TConst',273.15);                  % scalar or line vector of double: for pressure-and-temperature-constant property
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Determine length of names (number of species)
            obj.names = parser.Results.names;
            obj.N = numel(obj.names);
            obj.constMIXTURESwitch = parser.Results.constMIXTURESwitch;
            obj.pConst = parser.Results.pConst;
            obj.TConst = parser.Results.TConst;
            
            % Go through all species and construct instances of xMATERIAL
            % employing appending; also calculate values of
            % constant-p-T-properties cp, lambda, eta
            for i = 1:1:obj.N
                obj.materials = [obj.materials; xMATERIAL('name', parser.Results.names{i})];
                obj.materials(i).cPConst = obj.materials(i).specHeatCap('T',obj.TConst);
                obj.materials(i).lambdaConst = obj.materials(i).thermCond('T',obj.TConst);
                obj.materials(i).etaConst = obj.materials(i).dynVisc('T',obj.TConst);
            end
        
        end
        
        %% Method checking input validity
        function checkInput(obj,varargin)
            %CHECKINPUT checks input of class methods for validity
            %   The method checks the molar fraction handed over the class
            %   methods for validity (sum must be unity etc.)
            
            % Create input parser
            parser = inputParser;
            
            % Define input parser elements
            addParameter (parser,'T',273.15);    % double: Temperature in K
            addParameter (parser,'y',[]);        % array of double: Composition (molar fractions)

            % Parse for input
            parse (parser,varargin{:});
            
            T = parser.Results.T;
            y = parser.Results.y;
            
            % Check if implausible temperatures occur
            % lower bound set to 273.15K
            % only warn if actual temperature arrays are passed to MIXTURE
            % and not just single column vectors filled with ones
            % The temperatures get reset in the method
            if any(T<273.15) && ~exist('T')
                warning('TCurrRct < 273.15K!')
            end
            
            % Check input for validity
            % Non negative molar fractions
            for i=1:1:numel(y)
                if y(i) < -0.01
                    warning('Molar fractions smaller than -0.001 occured')
                    break
                end
            end
            
            %sum of molar fraction = 1
            %exception: y is empty. Empty y is passed to checkInput by
            %diffBin which does not use y
            %without the exception the error would be caused unnecessarily
            if (sum(y) < 0.99 | sum(y) > 1.01)  & ~isempty(y)
                warning('Sum of molar fractions unequal 1')
            end
            
            %length of y must equal N
            if size(y,1) ~= obj.N & ~isempty(y)
                error('Number of molar fractions does not match number of materials')
            end
            
        end
        
        %% Method returning the molar weight of the mixture in kg/mol
        function molarW = molarWeight(obj, varargin)
            %MOLARWEIGHT returns the molar weight of the mixture
            %   The method expects an input encoding the composition (molar fractions)
            %   Source:
            
            % Create input parser
            parser = inputParser;
            
            % Define input parser elements
            addParameter (parser,'y',[]);        % array of double: Composition (molar fractions)
            addParameter (parser,'switchCheckInput',true); %boolean: if true check input
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            y = parser.Results.y;
            switchCheckInput = parser.Results.switchCheckInput;
            
            % Set negative y to zero
            y(y<0) = 0;
            % Check input
            if switchCheckInput
                obj.checkInput('y',y)
            end
            
            % Identify number of columns of y (different compositions of
            % same mixture)
            numCols = size(y,2);
            
            % Pre-allocate molarWeight for all compositions
            molarW = zeros(1,numCols);
            
            % Calculate molarWeight
            for j = 1:numCols % Go through all compositions
                for i = 1:obj.N % Go through all species

                    % Add up specific molar weights weighed with their mole fraction
                    molarW(j) = molarW(j) + y(i,j)*obj.materials(i).molarWeight;

                end
            end
        end
        
        %% Method returning the specific heat capacity of the mixture in J/(kg*K)
        function c_p = specHeatCap(obj, varargin)
            %SPECHEATCAP returns the specific heat capacity of the mixture
            %   The method expects an input encoding the
            %   temperature and the composition (molar fractions)
            %   Source: VDI Heat Atlas D1 p. 158
            
            % Create input parser
            parser = inputParser;
            
            % Define input parser elements
            addParameter (parser,'T',273.15);    % double: Temperature in K
            addParameter (parser,'y',[]);        % array of double: Composition (molar fractions)
            addParameter (parser,'switchCheckInput',true); %boolean: if true check input
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            y = parser.Results.y;
            switchCheckInput = parser.Results.switchCheckInput;

            % Set negative y to zero
            y(y<0) = 0;
            % Check input
            if switchCheckInput
                obj.checkInput('T',T,'y',y)
            end

            % lower temperatures, especially negative ones, can lead to
            % errors in the methods and are therefore set to 273.15K
            T(T<273.15) = 273.15;
            
            % Go through all species
            c_p_molar = 0;
            for i = 1 : 1 : obj.N
                if ~obj.constMIXTURESwitch
                    % Add up specific molar heat capacities of the components
                    % weighed with their mole fraction. specHeatCap is
                    % calculated depending on current temperature
                    c_p_molar = c_p_molar + y(i)*obj.materials(i).specHeatCap('T',T)*obj.materials(i).molarWeight;
                else
                    % specHeatCap is not being calculated! value at
                    % constant Temperature is used
                    c_p_molar = c_p_molar + y(i)*obj.materials(i).cPConst*obj.materials(i).molarWeight;
                end
                
            end
            
            % Convert molar to gravimetric specific heat capacity of the
            % mixture
            if switchCheckInput
                c_p = c_p_molar./obj.molarWeight('y',y);
            else
                c_p = c_p_molar./obj.molarWeight('y',y,'switchCheckInput',false); 
            end
        end
        
        %% Method returning the dynamic viscosity of the mixture in Pa*s
        function eta = dynVisc(obj, varargin)
            %DYNVISC returns the dynamic viscosity of the mixture
            %   The method expects an input encoding the pressure, the
            %   temperature and the the composition (molar fractions)
            %   Source: VDI Heat Atlas D1 p. 166
            
            % Create input parser
            parser = inputParser;
            
            % Define input parser elements
            addParameter (parser,'p',1);      % double: Pressure in bar
            addParameter (parser,'T',273.15);    % double: Temperature in K
            addParameter (parser,'y',[]);        % array of double: Composition (molar fractions)
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            p = parser.Results.p;
            T = parser.Results.T;
            y = parser.Results.y;
            
            % Check input
            obj.checkInput('T',T,'y',y)
%             % Set negative y to zero
%             y(y<0) = 0;
            % lower temperatures, especially negative ones, can lead to
            % errors in the methods and are therefore set to 273.15K
            T(T<273.15) = 273.15;
            
            % Calculate auxiliary values
            F = zeros(obj.N);
            for i = 1 : 1 : obj.N
                for j = 1 : 1 : obj.N
                    if ~obj.constMIXTURESwitch
                        % Calculate auxiliary values
                        F(i,j) = (1+sqrt(obj.materials(i).dynVisc('T',T)/...
                                obj.materials(j).dynVisc('T',T)*...
                                sqrt(obj.materials(j).molarWeight/obj.materials(i).molarWeight)))^2/...
                                sqrt(8*(1+obj.materials(i).molarWeight/obj.materials(j).molarWeight));
                    else
                        F(i,j) = (1+sqrt(obj.materials(i).etaConst/...
                                obj.materials(j).etaConst*...
                                sqrt(obj.materials(j).molarWeight/obj.materials(i).molarWeight)))^2/...
                                sqrt(8*(1+obj.materials(i).molarWeight/obj.materials(j).molarWeight));
                    end
                end
            end
            
            % Calculate mixture dynamic viscosity
            eta = 0.0;
            for i = 1 : 1 : obj.N
                den = 0.0;
                for j = 1 : 1 : obj.N
                    den = den + y(j)*F(i, j);
                end
                if ~obj.constMIXTURESwitch
                    eta = eta + y(i)*obj.materials(i).dynVisc('T',T)/den;
                else
                    eta = eta + y(i)*obj.materials(i).etaConst/den;
                end
            end
            
        end
        
        %% Method returning the thermal conductivity of the mixture in W/m*K)
        function lambda = thermCond(obj, varargin)
            %THERMCOND returns the thermal conductivity of the mixture
            %   The method expects an input encoding the pressure, the
            %   temperature and the the composition (molar fractions)
            %   Source: VDI Heat Atlas D1 p. 169
            
            % Create input parser
            parser = inputParser;
            
            % Define input parser elements
            addParameter (parser,'p',1);      % double: Pressure in bar
            addParameter (parser,'T',273.15);    % double: Temperature in K
            addParameter (parser,'y',[]);        % array of double: Composition (molar fractions)
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            p = parser.Results.p;
            T = parser.Results.T;
            y = parser.Results.y;
            
            % Set negative y to zero
            y(y<0) = 0;
            % Check input
            obj.checkInput('T',T,'y',y)
            
            % lower temperatures, especially negative ones, can lead to
            % errors in the methods and are therefore set to 273.15K
            T(T<273.15) = 273.15;
            
            % Calculate auxiliary values
            % lambda is temperature dependent
            % therefore calculate F-array for each NxN component at each
            % temperature (point in reactor / x direction)
            
            if ~obj.constMIXTURESwitch
                F = zeros(obj.N,obj.N,length(T));
            else
                F = zeros(obj.N,obj.N);
            end
            for i = 1 : 1 : obj.N
                for j = 1 : 1 : obj.N
                    % Calculate auxiliary values
                    if ~obj.constMIXTURESwitch
                        F(i,j,:) = (1+sqrt(obj.materials(i).dynVisc('T',T)./...
                                obj.materials(j).dynVisc('T',T).*...
                                sqrt(obj.materials(j).molarWeight/obj.materials(i).molarWeight))).^2./...
                                sqrt(8*(1+obj.materials(i).molarWeight/obj.materials(j).molarWeight));
                    else
                        F(i,j) = (1+sqrt(obj.materials(i).etaConst/...
                             obj.materials(j).etaConst*...
                             sqrt(obj.materials(j).molarWeight/obj.materials(i).molarWeight)))^2/...
                             sqrt(8*(1+obj.materials(i).molarWeight/obj.materials(j).molarWeight));
                    end
                end
            end
            
            % Set breakpoint at F to stop when temperature vector and full
            % array of y is being considered
            if length(T) > 2
                if any(y(5,:)>1e-4)
                    F;
                end
            end
            
            % Calculate mixture thermal conductivity
            lambda = zeros(1,length(T));
            for k = 1:length(T)
                if ~obj.constMIXTURESwitch
                    for i = 1 : 1 : obj.N
                        den = 0.0;
                        % calculate denominator by summing up ober all
                        % y_j*F_ij
                        for j = 1 : 1 : obj.N
                            den = den + y(j,k).*F(i, j, k);
                        end
                        % Calculate lambda by summing up over all
                        % y_i*lambda_i
                        lambda(1,k) = lambda(1,k) + y(i,k).*obj.materials(i).thermCond('T',T(k))./den;
                    end
                else
                    for i = 1 : 1 : obj.N
                        den = 0.0;
                        for j = 1 : 1 : obj.N
                            den = den + y(j,k).*F(i, j, 1);
                        end
    %                     lambda(1,:) = (lambda + y(i,:)*obj.materials(i).lambdaConst./(repmat(den(1),1,size(y,2))));%.*ones(1,length(T));
                        lambda(1,k) = lambda(1,k) + y(i,k)*obj.materials(i).lambdaConst./den;
                    end
                end
            end
            
        end

        %% Method returning the binary diffusivites of all components in the mixture
        function diffBin = diffusivityBin(obj, varargin)
            %DIFFUSIVITYBIN returns the binary diffusivites of all components of the mixture in
            %   each component in the mixture in m^2/s at a given pressure
            %   and temeprature.
            %   The method expects an input encoding temperature, pressure,
            %   index of species A
            %   Source: Equation by Fuller et. al. according to
            %           Reid, The Properties of Gases and Liquids,p.587-588
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'p',1);         % double: pressure in bar
            addParameter (parser,'T',273.15);    % double: Temperature in K
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            p = parser.Results.p;
            T = parser.Results.T;
            
            % Check input
            obj.checkInput('T',T)
%             % Set negative y to zero
%             y(y<0) = 0;
            % lower temperatures, especially negative ones, can lead to
            % errors in the methods and are therefore set to 273.15K
            T(T<273.15) = 273.15;
            
            % Calculate & return binary diffusivities of each component of
            % mixture in every other component of mixture (=binary
            % diffusivity in other components, as well as self diffusivity)
            
            % if block for computation time reduction:
            % except for parameter T.^1.75 the binary diffusivitys are
            % independent of T (and y). Therefore calculate array with
            % loop only if obj.diffBinArray is empty and write it into MIXTURE properties. 
            % If its not empty (array exists) only multiply it with current Temperature T
            
% % % % % % % % % OLD COMPUTATION ROUTINE, acc. to Poling, The properties of gases and liquids, equ. 11-4.4            
            if isempty(obj.diffBinArray)
                obj.diffBinArray = zeros(obj.N, obj.N); % pre-allocating result array
                for i = 1:obj.N % rows of result array
                    for j = 1:obj.N % columns of result array
                        obj.diffBinArray(i,j) = (1/(1000*obj.materials(i).molarWeight) + 1/(1000*obj.materials(j).molarWeight))^0.5 / ...
                        (p * 2^0.5 * ...
                        (obj.materials(i).diffVol^(1/3) + obj.materials(j).diffVol^(1/3))^2);
                    end
                end
                diffBin = 0.00143*T.^1.75 * obj.diffBinArray * 1e-4; % calcualte diffBin with new diffBinArray from loop, Unit conversion to m^2/s
            else
                diffBin = 0.00143*T.^1.75 * obj.diffBinArray * 1e-4; % calculate diffBin with diffBinArray from first run of method, Unit conversion to m^2/s
            end
% % % % % % % % % 

% % % % % % % % % % % NEW COMPUTATION ROUTINE, acc. to Jess, Chemical Technology, p58, and Diss. Hannah Kirsch
% % % % % % % % % % % tabulated diffusion volumes differ a little from values in Database.xlsx
% %             if isempty(obj.diffBinArray)
% %                 obj.diffBinArray = zeros(obj.N, obj.N); % pre-allocating result array
% %                 for i = 1:obj.N % rows of result array
% %                     for j = 1:obj.N % columns of result array
% %                         obj.diffBinArray(i,j) = ((obj.materials(i).molarWeight + obj.materials(j).molarWeight)/(obj.materials(i).molarWeight * obj.materials(j).molarWeight))^0.5 / ...
% %                                                 (p*1e5 * ((1e-6*obj.materials(i).diffVol)^(1/3) + (1e-6*obj.materials(j).diffVol)^(1/3))^2);
% %                     end
% %                 end
% %                 diffBin = 3.16e-8*T.^1.75 * obj.diffBinArray; % calcualte diffBin with new diffBinArray from loop
% %             else
% %                 diffBin = 3.16e-8*T.^1.75 * obj.diffBinArray; % calculate diffBin with diffBinArray from first run of method
% %             end
% % % % % % % % % % 
            
        end
        
        %% Method returning the diffusivity of each component in the mixture in m^2/s
        function diffMix = diffusivityMix(obj, varargin)
            %DIFFUSIVITYMIX returns the diffusivity of each component of
            %the mixture in the mixture
            %   The method expects an input encoding temperature, pressure,
            %   molar fractions
            %   Source: Jess, Wasserscheid - Chemical Technology, p.58-59
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'p',1);         % double: pressure in bar
            addParameter (parser,'T',273.15);    % double or array of double: Temperature in K
            addParameter (parser,'y',[]);        % array of double: Composition (molar fractions)
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            p = parser.Results.p;
            T = parser.Results.T;
            y = parser.Results.y;
            
            % Set negative y to zero
            y(y<0) = 0;
            % Check input
            obj.checkInput('T',T,'y',y)
            
            % lower temperatures, especially negative ones, can lead to
            % errors in the methods and are therefore set to 273.15K
            if any(T<0)
                beep
                disp('T < 0K!')
            end
            T(T<273.15) = 273.15;
            
            % Calculate binary diffusivity, if Temperature for properties is constant
            if obj.constMIXTURESwitch
                diffBinConst = obj.diffusivityBin('p',obj.pConst,'T',obj.TConst);
            end
            
            % Calculate & return diffusivity in an NxX array
            % for each specie in the mixture N
            % for each temperature value X = for each node = set of
            % temperature and vector of molar fractions
            diffMix = zeros(obj.N,length(T));
            for j = 1:length(T) % index for the case that function is called with an array for temperature, e.g. a value for each node of a mesh for PDE calculation
                
                if ~obj.constMIXTURESwitch
                    diffBin = obj.diffusivityBin('p', p,'T',T(j)); % Calculate binary diffusivitys for Temperature T(j)
                else
                    diffBin = diffBinConst;
                end
                
                for k = 1:obj.N % indexing number of species
                        denomArray = zeros(1,obj.N); % Preallocate array of denominator
                        for i = 1:obj.N % indexing number of species
                            if k ~= i % prevent eigendiffusivity from being considered
                                denomArray(i) = (y(i,j) / diffBin(k,i)); % Array in denominator of the equation, each element represents one term of the sum
                            end
                        end
                    denomSum = sum(denomArray); % sum in denominator
                    diffMix(k,j) = (1 - y(k,j)) / denomSum; % result for diffusivity
                end
            end
        end % function
    end % methods
end % classdef

