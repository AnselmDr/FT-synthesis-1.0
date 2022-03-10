classdef xMATERIAL < handle
    %XMATERIAL can be employed to describe thermodynamical properties
    %   Instances of this class can be used to calculate thermodynamic
    %   properties of variours materials. It comprises methods that
    %   calculate thermodynamic properties using correlations from the VDI
    %   Heat atlas. The parameters needed for the calculations are extracted
    %   from the corresponding chapters in the VDI Heat Atlas and are
    %   stored in the instances properties.
    
    properties
        material        % string: material designation
        molarWeight     % double: molar weight in kg/mol
        R = 8.3144626182% double: universal gas constant in J/(mol*K)
        N_A = 6.02214e23% double: avogadro's constant / mol^-1
        k_B = 1.3806e-23% double: Boltzmann's constant = R/N_A / J/K = m2 kg s-2 K-1
        
        % Gaseous material properties
        cpParam         % line vector of double: paramters needed for the 
                        % calculation of the specfic heat capacity
        lambdaParam     % line vector of double: paramters needed for the 
                        % calculation of the thermal conductivity
        etaParam        % line vector of double: paramters needed for the 
                        % calculation of the dynamic viscosity
        diffVol         % double: diffusion volume [10^6*m^3/mol] needed for the
                        % calculation of the diffusivity
        henryParam      % line vector of double: parameters to calculate
                        % coefficient of gas in C30
                        
        % Liquid material properties
        HevapParam      % line vector of double: paramters needed for the 
                        % calculation of the evaporation enthalpy
        cParam          % line vector of double: paramters needed for the 
                        % calculation of the liquid specific heat capacity
        TCrit           % double: critical temperatur ein K
        pCrit           % double: critical pressure in bara
        rhoCrit         % double: critical density in kg/m^3
        pVapParam       % line vector of double: paramters needed for the 
                        % calculation of the vapor pressure
        areaColl        % double: collision area / nm^2
        rhoBp           % double: molar density at boiling point in mol/m^3
                        % to calculate diffusion of gas in C30
                        
        % Solid material properties
        Rp02Param       % line vector of double: paramters needed for the 
                        % calculation of the 0.2% proof strength
        EParam          % line vector of double: paramters needed for the 
                        % calculation of the modulus of elasticity
        lambdaSolidParam% line vector of double: paramters needed for the 
                        % calculation of the solid thermal conductivity
        cSolidParam     % line vector of double: paramters needed for the 
                        % calculation of the solid specific heat capacity
                        
        % Permeable sintered materials properties
        porosity        % double / -: porosity
        coeffsPerm      % double column vector length 2:
                        % - m^2: permeability coefficient of volume Psi_v according
                        %        DIN EN ISO 4022 alpha acc. GKN SIKA Porous Axial
                        % - m:   permeability coefficient of inertia according
                        %        DIN EN ISO 4022 beta acc. GKN SIKA Porous Axial
        poreSizeDist    % double column vector length 4 / µm: pore size distribution
                        % values acc. GKN SIKA Porous Axial: d_min, d_mean,
                        % d_max, d_eq
                        
        % Properties at constant T and if applicable p; scalar or line
        % vector double
        cPConst         % gaseous specific heat capacity in J/(kg*K)
        lambdaConst     % gaseous thermal conductivity in W/(m*K)
        etaConst        % gaseous dynamic viscosity in Pa*s = kg/(m*s)
        %%% reserved for liquid properties currently unused %%%
        Rp02Const       % 0.2% proof strength in MPa
        EConst          % modulus of elasticity in MPa?
        lambdaSConst    % solid thermal conductivity in W/(m*K)
        cSConst         % solid specific heat capacity in j/(kg*K)
                        
        % For displaying purposes
        valuest         % double line vector for typical values of
                        % temperature in °C from VDI Heat Atlas
        valuesT         % double line vector for typical values of
                        % temperature in K from VDI Heat Atlas
        % TPlot is a double line vector for creating figures of
        % thermophysical property derived from TLow, THigh, TReso.
        TLow
        THigh
        TResol
        TPlot
        % pPlot is a double line vector for creating figures of
        % thermophysical property derived from pLow, pHigh, pReso.
        pLow
        pHigh
        pNumber
        pPlot
    end
    
    methods
        %% Constructor
        function obj = xMATERIAL(varargin)
            %XMATERIAL Construct an instance of this class
            %   The constructor expects a string input that encodes the material
            %   that the constructed instance shall represent. The
            %   constructor reads a csv file containing all parameters
            %   needed for the calculations of material properties,
            %   searches for the entries that correspond to the instance
            %   material and stores them in the instance properties. If the
            %   designation is not found, an error is displayed.
           
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'name','noName'); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign property values
            obj.material = parser.Results.name;
            
            % Import data from file Database.csv
            delimiterIn = ';';
            headerlinesIn = 8;
            a = importdata('Database.csv', delimiterIn, headerlinesIn);
            %Size of array
            [lMax, ~] = size(a.data);
            
            % Search for material in the first column when found: set crit
            % to true and stop
            crit = false; % initialize crit to be false
            for i=1:1:lMax
                crit = strcmp(a.textdata(headerlinesIn+i,1), obj.material);
                if crit == true
                    break
                end
            end
            
            % Display error if the material was not found
            if crit == false
                % Return error
                error('Material designation not found in CSV file.')
            end
                      
            % Assign value of molar weight
            obj.molarWeight = a.data(i,1);
            
            % Assign parameters for cp for the respective material
            obj.cpParam = a.data(i,2:8);
            
            % Assign parameters for thermal conductivity for the respective material
            obj.lambdaParam = a.data(i,9:13);
            
            % Assign parameters for dynamic viscosity for the respective material
            obj.etaParam = a.data(i,14:18);
            
            % Assign parameters for evaporation enthalpy
            % VDI Heat Atlas D3.1 Tab. 4 S. 385
            obj.HevapParam = a.data(i,19:23);
            
            % Assign parameters for liquid specific heat capacity
            % VDI Heat Atlas D3.1 Tab. 5 S. 395
            obj.cParam = a.data(i,24:29);
            
            % Assign critical calues
                % Critical temperature in K
                obj.TCrit = a.data(i, 30);
                % Critical pressure in bara
                obj.pCrit = a.data(i, 31);
                % Critical density in kg/m^3
                obj.rhoCrit = a.data(i, 32);
            
            % Assign values for vapor pressure
            % VDI Heat Atlas D3.1 Tab. 3 on p. 375
            obj.pVapParam = a.data(i, 33:36);
            
            % Coefficients for calcutlation of 0,2% Dehngrenze
            % in N/mm^2 from T in K (valid interval 77,15K - 823,15 K)
            % with polynome of order 3 (A+B*T+…) (Source: WIAM online
            % database, mechanisch-technologische Eigenschaften)
            obj.Rp02Param = a.data(i, 37:41);
            
            % Coefficients for Elastizitätsmodul in kN/mm^2 from T in K
            % (valid interval 173,15 - 1273,15 K) with polynom of order 1
            % (straight line A+B*T) (Source: WIAM online database,
            % physikalische Eigenschaften (Source SEW 310:1992-08)
            obj.EParam = a.data(i, 42:43);
            
            % Coefficients for thermal conductivity in W/(m*K) from T in K
            % (valid interval 293,15 - 1073,15 K) with polynom of order 1
            % (straight line A+B*T) (Source: WIAM online database,
            % physikalische Eigenschaften (Source SEW 310:1992-08)
            obj.lambdaSolidParam = a.data(i, 48:49);
            
            % Coefficients for actual specific heat capacity in 10^3 J/(kg*K)
            % from T in K (valid interval 173,15 - 1273,15 K) with polynome
            % of order 3 (A+B*T+…) (Source: WIAM online database,
            % physikalische Eigenschaften (Source SEW 310:1992-08)
            obj.cSolidParam = a.data(i, 52:55);
            
            % Assign diffusivity volume for the respective material
            obj.diffVol = a.data(i,56);
            
            % Assign collision area for the respective material
            obj.areaColl = a.data(i,65);
            
            % Assign properties for the respective frit
            obj.porosity = a.data(i,66);
            obj.coeffsPerm = [a.data(i,67)*1e-12; a.data(i,68)*1e-7];
            obj.poreSizeDist = a.data(i,69:72)';
            
            % Coefficients for Henry Coefficient
            obj.henryParam = a.data(i,73:78);

            % molar densities at boiling point
            obj.rhoBp = a.data(i,79);
            
            % Assign parameters for new property for the respective material
            %obj.newProp = a.data(i,59:60);
            
            obj.valuest = [-200 -175 -150 -125 -100 -75 -50 -25 0 25 50 75 100 125 150 200 300 400 500 600 700 800 900 1000];
                            % double line vector for typical values of
                            % temperature in K from VDI Heat Atlas
            obj.valuesT = obj.valuest + 273.15;
            
            % TPlot is a double line vector for creating figures of thermophysical property
            obj.TLow = 0;
            obj.THigh = obj.TLow + 1000 + 273.15;
            obj.TResol = 100;
            obj.TPlot = obj.TLow:(obj.THigh-obj.TLow)/obj.TResol:obj.THigh;
            
            % pPlot is a double line vector for creating figures of thermophysical property
            obj.pLow = -3;
            obj.pHigh = 2;
            obj.pNumber = length(obj.TPlot);
            obj.pPlot = logspace(obj.pLow, obj.pHigh, obj.pNumber);
        end
        
        %% Method returning the density in kg/m^3
        function rho = density(obj, varargin)
            %DENSITY returns the density of the material according to ideal
            %gas law
            %   Input arguments
            %   - double scalar (or line vector): pressure in bar
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): density in kg/mol
            %   Source: VDI Heat Atlas D1 p. 150 Eq. 24
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'p',1.013);
            addParameter (parser,'T',273.15);
            addParameter (parser,'plot',false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            p = parser.Results.p;
            T = parser.Results.T;
            checkPlot = parser.Results.plot;
            
            % Calculate & return rho
            rho = p.*1e5.*obj.molarWeight./(obj.R.*T);
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.density('p',p,'T',obj.TPlot));
                title(obj.material);
                text(800, 0.7*max(obj.density('p',p,'T',obj.TPlot)), cat(2, 'p = ', num2str(p), ' bar'),'FontSize',16);
                xlabel('Temperature / K');
                ylabel('Density / kg/m^3');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.density('p',p,'T',obj.valuesT)]);
            end
            
        end
        
        %% Method returning the specific heat capacity in J/(kg*K)
        function c_p = specHeatCap(obj, varargin)
            %SPECHEATCAP returns the specific heat capacity of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): specific heat capactiy in
            %   J/kgK
            %   Source: VDI Heat Atlas D3 p. 358 Eq. 10 and D1, S. 157, Gl. 51a
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'plot',false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            checkPlot = parser.Results.plot;
            
            % Calculate auxiliary value
            A = obj.cpParam(1);
            B = obj.cpParam(2);
            C = obj.cpParam(3);
            D = obj.cpParam(4);
            E = obj.cpParam(5);
            F = obj.cpParam(6);
            G = obj.cpParam(7);
            
            % Calculate and return c_p gravimetric
            c_p_molar_perR = B + ...
                            ( C - B ) .* (T./(A+T)).^2 .* ...
                                ( 1 - A./(A+T) .* ( D + E.*T./(A+T) + F.*(T./(A+T)).^2 + G.*(T./(A+T)).^3) );
            c_p = c_p_molar_perR.*obj.R./obj.molarWeight;
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.specHeatCap('T',obj.TPlot));
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('Specific heat capacity / J/(kg*K)');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.specHeatCap('T',obj.valuesT)]);
            end
            
        end

        %% Method returning the thermal conductivity in W/(m*K)
        function lambda = thermCond(obj, varargin)
            %THERMCOND returns the thermal conductivity of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): thermal conductivity in W/(m*K)
            %   Source: VDI Heat Atlas D3 p. 357 Eq. 5 and D1, S. 167, Gl. 103
            %   Later the pressure correction might be implemented
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'plot',false);
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            checkPlot = parser.Results.plot;
            
            % Calculate auxiliary value
            A = obj.lambdaParam(1)*1e-3;
            B = obj.lambdaParam(2)*1e-3;
            C = obj.lambdaParam(3)*1e-6;
            D = obj.lambdaParam(4)*1e-9;
            E = obj.lambdaParam(5)*1e-12;
            
            % Calculate & return lambda
            lambda = A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.thermCond('T',obj.TPlot));
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('Thermal conductivity / W/(m*K)');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.thermCond('T',obj.valuesT)]);
            end
            
        end

        %% Method returning the dynamic viscosity in Pa*s
        function eta = dynVisc(obj, varargin)
            %DYNVISC returns the dynamic viscosity of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): dynamic viscosity in Pa*s
            %   Source: VDI Heat Atlas D3 p. 357 Eq. 3 COrrections in D1
            %   Later the pressure correction might be implemented
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'checkPlot',[]); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            checkPlot = parser.Results.checkPlot;
            
            % Calculate auxiliary value
            A = obj.etaParam(1)*1e-5;
            B = obj.etaParam(2)*1e-7;
            C = obj.etaParam(3)*1e-10;
            D = obj.etaParam(4)*1e-12;
            E = obj.etaParam(5)*1e-15;
            
            % Calculate & return eta
            eta = A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;
            
            % If applicable plot figure and value chart
            if checkPlot == true
                plot(obj.TPlot, obj.dynVisc('T',obj.TPlot));
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('Dynamic viscosity / Pa*s');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.dynVisc('T',obj.valuesT)]);
            end
            
        end

        %% Method returning the 0.2% proof strength in N/mm^2=MPa
        function Rp02 = pStrength02(obj, varargin)
            %PSTRENGTH02 returns the 0.2% proof strength of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): 0.2% proof strength in N/mm^2=MPa
            %   Source:  (valid interval 77,15K - 823,15 K)
            %   with polynome of order 4 (A+B*T+…) (Source: WIAM online
            %   database, mechanisch-technologische Eigenschaften)
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'plot',false); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            checkPlot = parser.Results.plot;
            
            % Calculate auxiliary value
            A = obj.Rp02Param(1);
            B = obj.Rp02Param(2);
            C = obj.Rp02Param(3);
            D = obj.Rp02Param(4);
            E = obj.Rp02Param(5);
            
            % Calculate & return new property
            Rp02 = A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.pStrength02('T',obj.TPlot));
                title(obj.material);
                axis([77 823 0 400]);
                xlabel('Temperature / K');
                ylabel('0.2% proof strength / MPa');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.pStrength02('T',obj.valuesT)]);
            end
            
        end
        
        %% Method returning the solid thermal conductivity in W/mK
        function lambdaS = thermCondS(obj, varargin)
            %THERMCONDS returns the solid thermal conductivity of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): solid thermal conductivity in W/mK
            %   (valid interval 293,15 - 1073,15 K) with polynom of order 1
            %   (straight line A+B*T) (Source: WIAM online database,
            %   physikalische Eigenschaften (Source SEW 310:1992-08)
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'plot',false); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            checkPlot = parser.Results.plot;
            
            % Calculate auxiliary value
            A = obj.lambdaSolidParam(1);
            B = obj.lambdaSolidParam(2);
            
            % Calculate & return new property
            lambdaS = A+ B.*T;
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.thermCondS('T',obj.TPlot));
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('Thermal conductivity of solid / W/(m*K)');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.thermCondS('T',obj.valuesT)]);
            end
            
        end

        %% Method returning the solid specific heat capacity in J/(kg*K)
        function cS = specHeatCapS(obj, varargin)
            %SPECHEATCAPS returns the solid specific heat capacity of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   Output arguments
            %   - double scalar (or line vector): actual specific heat
            %   capacity in J/(kg*K) from T in K (valid interval 173,15 -
            %   1273,15 K) with polynome of order 3 (A+B*T+…) (Source:
            %   WIAM online database, physikalische Eigenschaften
            %   (Source SEW 310:1992-08)
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'plot',false); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            checkPlot = parser.Results.plot;
            
            % Calculate auxiliary value
            A = obj.cSolidParam(1);
            B = obj.cSolidParam(2);
            C = obj.cSolidParam(3);
            D = obj.cSolidParam(4);
            
            % Calculate & return new property
            cS = A + B.*T + C.*T.^2 + D.*T.^3;

            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.specHeatCapS('T',obj.TPlot));
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('solid specific heat capacity / J/(kg*K)');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.specHeatCapS('T',obj.valuesT)]);
            end
            
        end
        
         %% Method returning the mean free path of a gas
        function MFP = meanFreePath(obj, varargin)
            %MEANFREEPATH returns the mean free path of the material
            %   Input arguments
            %   - double scalar (or line vector): temperature in K
            %   - double scalar (or line vector): pressure in bar
            %   Output arguments
            %   - double scalar (or line vector): mean free path in m
            %   Source: Atkins, Peter W.; Paula, Julio de; Keeler, James J.
            %   (2022): Physikalische Chemie. Sechste Auflage. Weinheim, Germany: Wiley-VCH. 
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'p',1.0);
            addParameter (parser,'plot',false); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            p = parser.Results.p;
            checkPlot = parser.Results.plot;
           
            % Calculate & return new property
            MFP = obj.k_B.*T./(sqrtm(2.0).*obj.areaColl.*1e-18.*p.*1e5);
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, values);
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('Mean free path / m');
                
                format shortG;
                %disp([obj.valuest;obj.valuesT;obj.meanFreePath('T',obj.valuesT)]);
            end
            
        end
        
         %% Method returning the mean free path of a gas
        function meanFreePathPlot(obj, varargin)
            %MEANFREEPATHPLOT see mean free path.
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'save',false); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            saveSwitch = parser.Results.save;
            
            colorMap = parula(5);
            pressure = [0.1 0.5 1 5 10 50];

            for i = 1:5
                MFP(i,:) = obj.meanFreePath('T',obj.TPlot, 'p', pressure(i));
                hold on
                plot(obj.TPlot, MFP(i,:),'Color',colorMap(i,:),'Linewidth',2)
            end
            hold off
            xlabel('Temperature T / K')
            ylabel('Mean free path / m')
            title(strcat('Mean free path of ',obj.material))
            set(gca,'linewidth',2)
            lgd = legend(string(pressure));
            lgd.Title.String = 'Pressure / bar';
            legend('Location','best');
            
            if saveSwitch
                if ~isfolder(strcat('results/',char(obj.sDate)))
                    mkdir('results/',obj.sDate)
                end
                filename = strcat('results/',obj.sDate,'/',obj.sDateTime, '_meanFreePath.png');
                saveas(gcf, filename);
            end
            
        end
        
        
        %% Method returning the Henry Constant in C30 in bar
        function HiL = henryConst(obj, varargin)
            %HENRYCOEFF returns the the henry coefficient HiL for a gaseous
            %   component in a liquid alkane
            %   Input arguments: Temperature in K (double or line vector)
            %   and number of Carbons in the alkane nC (double)
            %   Output arguments: double or linevector 
            %   Source: 10.1016/S0378-3812(97)00166-0
            %   Check for validity intervals of alkanes and temperatures
            %   ATTENTION: ONLY PARAMETERS FOR H2 and CO in Database.csv
            
            % Create input parser
            parser = inputParser;

            % Define input parser elements
            addParameter (parser,'T',273.15);
            addParameter (parser,'nC',30);
            addParameter (parser,'plot',false); 
            
            % Parse for input
            parse (parser,varargin{:});
            
            % Assign values
            T = parser.Results.T;
            % number of Carbons in C30
            nC = parser.Results.nC; 
            checkPlot = parser.Results.plot;
            
            % rewrite parameters for readability
            A = obj.henryParam(1);
            B = obj.henryParam(2);
            C = obj.henryParam(3);
            D = obj.henryParam(4);
            E = obj.henryParam(5);
            DeltaH = obj.henryParam(6);
            
            % Calculate Hi0
            Hi0 = A + B./T + C.*log(T) + D.*T.^2 + E./(T.^2);
            
            % Calculate Henry Constant
            HiL = exp(Hi0 - nC * DeltaH);
            
            % If applicable plot figure and value chart
            if checkPlot
                plot(obj.TPlot, obj.newProp('T',obj.TPlot));
                title(obj.material);
                xlabel('Temperature / K');
                ylabel('HiC30 / bar');
                
                format shortG;
                disp([obj.valuest;obj.valuesT;obj.newProp('T',obj.valuesT)]);
            end
            
        end
         %% Method returning the new property in new unit
%         function nP = newProp(obj, varargin)
%             %NEWPROP returns the new property of the material
%             %   Input arguments - double scalar (or line vector):
%             %   temeprature in K Output arguments - double scalar (or line
%             %   vector): new property in new unit Source: Anything e.g. VDI
%             %   Heat Atlas D3 p. 357 Eq. 3
%             
%             % Create input parser
%             parser = inputParser;
% 
%             % Define input parser elements
%             addParameter (parser,'T',273.15);
%             addParameter (parser,'plot',false); 
%             
%             % Parse for input
%             parse (parser,varargin{:});
%             
%             % Assign values
%             T = parser.Results.T;
%             checkPlot = parser.Results.plot;
%             
%             % Calculate auxiliary value
%             A = obj.propParam(1)*1e-5;
%             B = obj.propParam(2)*1e-7;
%             C = obj.propParam(3)*1e-10;
%             D = obj.propParam(4)*1e-12;
%             E = obj.propParam(5)*1e-15;
%             
%             % Calculate & return new property
%             nP = A + B.*T + C.*T.^2 + D.*T.^3 + E.*T.^4;
%             
%             % If applicable plot figure and value chart
%             if checkPlot
%                 plot(obj.TPlot, obj.newProp('T',obj.TPlot));
%                 title(obj.material);
%                 xlabel('Temperature / K');
%                 ylabel('New property / new unit');
%                 
%                 format shortG;
%                 disp([obj.valuest;obj.valuesT;obj.newProp('T',obj.valuesT)]);
%             end
%             
%         end
        
    end
    
end

