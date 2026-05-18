function price = priceGEOSatellite(satellite)
    % Function to estimate price of GEO satellites based on various heuristics
    % Input: satellite structure with all fields
    % Output: price in dollars (or NaN if cannot be priced)
    
    % Initialize price as NaN
    price = NaN;
    
    % Extract relevant fields for easier access
    mass = satellite.LaunchMass_kg;
    purpose = satellite.Purpose;
    operator = satellite.Operator;
    users = satellite.Users;
    name = satellite.Name;
    
    %% Case 1: Military/Intelligence satellites
    if contains(users, 'Military', 'IgnoreCase', true) || ...
       contains(operator, 'Defense', 'IgnoreCase', true) || ...
       contains(operator, 'NRO', 'IgnoreCase', true) || ...
       contains(name, 'USA', 'IgnoreCase', true)
        
        % High-end military satellites
        if contains(purpose, 'Electronic Intelligence', 'IgnoreCase', true) || ...
           contains(purpose, 'Surveillance', 'IgnoreCase', true)
            price = 1.5e9;  % $1.5 billion for SIGINT/ELINT
            
        elseif contains(name, 'WGS', 'IgnoreCase', true) || ...
               contains(name, 'Wideband Global', 'IgnoreCase', true)
            price = 500e6;  % $500M for WGS satellites
            
        elseif contains(name, 'MUOS', 'IgnoreCase', true)
            price = 1.2e9;  % $1.2B for MUOS
            
        elseif contains(name, 'AEHF', 'IgnoreCase', true)
            price = 1.8e9;  % $1.8B for AEHF
            
        elseif contains(name, 'SBIRS', 'IgnoreCase', true)
            price = 1.5e9;  % $1.5B for SBIRS
            
        elseif contains(name, 'Milstar', 'IgnoreCase', true) || ...
               contains(name, 'DSCS', 'IgnoreCase', true)
            price = 800e6;  % $800M for older military comsats
            
        else
            % Generic military satellite pricing based on mass
            if ~isnan(mass)
                price = mass * 250000;  % $250k per kg for military
            else
                price = 600e6;  % Default $600M for unknown military
            end
        end
        return;
    end
    
    %% Case 2: Large commercial GEO comsats (by mass)
    if ~isnan(mass)
        if mass >= 6000
            % Very large commercial satellites
            price = 400e6;  % $400M
            
        elseif mass >= 5000
            % Large satellites (Intelsat, SES, etc.)
            price = 300e6;  % $300M
            
        elseif mass >= 4000
            price = 250e6;  % $250M
            
        elseif mass >= 3000
            price = 200e6;  % $200M
            
        elseif mass >= 2000
            price = 150e6;  % $150M
            
        elseif mass >= 1000
            price = 100e6;  % $100M
            
        else
            % Smaller satellites
            price = mass * 100000;  % $100k per kg
        end
        
        % Adjust for specific operators (premium or budget)
        if contains(operator, 'Intelsat', 'IgnoreCase', true) || ...
           contains(operator, 'SES', 'IgnoreCase', true) || ...
           contains(operator, 'EUTELSAT', 'IgnoreCase', true)
            price = price * 1.1;  % 10% premium for major operators
        end
        
        return;
    end
    
    %% Case 3: Purpose-based pricing (when mass is unknown)
    if contains(purpose, 'Communications', 'IgnoreCase', true)
        price = 250e6;  % Default $250M for GEO comsat
        
    elseif contains(purpose, 'Navigation', 'IgnoreCase', true) || ...
           contains(purpose, 'GPS', 'IgnoreCase', true)
        price = 500e6;  % $500M for navigation satellites
        
    elseif contains(purpose, 'Meteorology', 'IgnoreCase', true) || ...
           contains(purpose, 'Earth Observation', 'IgnoreCase', true)
        price = 400e6;  % $400M for weather satellites
        
    elseif contains(purpose, 'Technology Development', 'IgnoreCase', true)
        price = 150e6;  % $150M for tech demo
        
    elseif contains(purpose, 'Space Science', 'IgnoreCase', true)
        price = 300e6;  % $300M for science missions
        
    else
        % Unknown purpose
        price = 200e6;  % Default $200M
    end
    
    %% Case 4: Specific satellite families/contractors (override if needed)
    if contains(name, 'GOES', 'IgnoreCase', true)
        price = 550e6;  % GOES satellites are expensive
        
    elseif contains(name, 'Meteosat', 'IgnoreCase', true)
        price = 400e6;
        
    elseif contains(name, 'Beidou', 'IgnoreCase', true) || ...
           contains(name, 'IRNSS', 'IgnoreCase', true) || ...
           contains(name, 'QZS', 'IgnoreCase', true)
        price = 250e6;  % Regional navigation systems
        
    elseif contains(name, 'TDRS', 'IgnoreCase', true)
        price = 600e6;  % NASA data relay satellites
        
    elseif contains(name, 'Express', 'IgnoreCase', true) && ...
           contains(operator, 'Russia', 'IgnoreCase', true)
        if ~isnan(mass) && mass > 2000
            price = 150e6;  % Larger Russian comsats
        else
            price = 100e6;  % Smaller Russian comsats
        end
        
    elseif contains(name, 'Galaxy', 'IgnoreCase', true) || ...
           contains(name, 'AMC', 'IgnoreCase', true)
        price = 200e6;  % Intelsat/PanAmSat fleet
        
    elseif contains(name, 'GSAT', 'IgnoreCase', true)
        price = 100e6;  % Indian satellites tend to be cheaper
    end
    
end