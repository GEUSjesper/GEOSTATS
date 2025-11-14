function [m,layer_index,flag_vector] = prior_m_reals(INFO,z,flag_vector)


% Number of sections
N = INFO.SECTIONS.N_sections;


% Initialize lithology vector
m = randsample([INFO.SECTIONS.types{N} INFO.SECTIONS.types{N}],1)*ones(size(z));


% Initialize layer vector
layer_count = 1;
layer_index = layer_count*ones(size(z));
layer_count = layer_count + 1;
if N == 1
    return
end


% random  vector for occurence of layers
r = rand(1,N-1);


% Prelocate vectors
thick_sections = zeros(1,N);
N_layers = zeros(1,N-1);
types_layers = cell(1,N-1);
thick_layers = cell(1,N-1);
for i = 1:N-1
    if r(i) <= INFO.SECTIONS.frequency(i)
        thick_sections(i) = rand(1)*(INFO.SECTIONS.max_sections(i)-INFO.SECTIONS.min_sections(i))+INFO.SECTIONS.min_sections(i);              % Randomly draw thickness of sections
        N_layers(i) = randi(INFO.SECTIONS.max_layers(i) - INFO.SECTIONS.min_layers(i) + 1) + INFO.SECTIONS.min_layers(i) - 1;               % Randomly draw number of layers
        types_layers{i} = randsample(repmat(INFO.SECTIONS.types{i},1,N_layers(i)+1),N_layers(i));                                                  % Randomly draw lithology types
        thick_layers{i} = rand(1,N_layers(i)).*(INFO.TYPES.max_thick(types_layers{i}) - INFO.TYPES.min_thick(types_layers{i})) + INFO.TYPES.min_thick(types_layers{i});    % Randomly draw thicknesses of layers  
    elseif r(i) > INFO.SECTIONS.frequency(i)
        thick_sections(i) = 0;
        N_layers(i) = 0;
        types_layers{i} = [];
        thick_layers{i} = [];
    end
end                  


% Make sure thicknesses fit within their intervals
if N > 1
    for i = find(thick_sections ~= 0)
        thick_layers{i} = thick_layers{i}./(sum(thick_layers{i})/thick_sections(i));
    end
end


% Check if dimensions fit
tries = 1;
checksum = 0;
for i = find(thick_sections ~= 0)
    layers_max_check = sum(thick_layers{i} >= 1.1*INFO.TYPES.max_thick(types_layers{i}));
    layers_min_check = sum(thick_layers{i} <= (1/1.1)*INFO.TYPES.min_thick(types_layers{i})); 
    checksum = checksum + layers_max_check + layers_min_check;
end


% Do it all again if dimensions don't fit
while checksum > 0
    for i = 1:N-1
        if r(i) <= INFO.SECTIONS.frequency(i)
            % if tried enough drawn a different number of layers and write a warning
            if tries > 1000
                N_layers(i) = randi(INFO.SECTIONS.max_layers(i) - INFO.SECTIONS.min_layers(i) + 1) + INFO.SECTIONS.min_layers(i) - 1; 
                flag_vector(2) = 1;
            end
            thick_sections(i) = rand(1)*(INFO.SECTIONS.max_sections(i)-INFO.SECTIONS.min_sections(i))+INFO.SECTIONS.min_sections(i);              % Randomly draw thickness of sections            N_layers(i) = randi(INFO.SECTIONS.max_layers(i) - INFO.SECTIONS.min_layers(i) + 1) + INFO.SECTIONS.min_layers(i) - 1;               % Randomly draw number of layers
            types_layers{i} = randsample(repmat(INFO.SECTIONS.types{i},1,N_layers(i)+1),N_layers(i));                                                  % Randomly draw lithology types
            thick_layers{i} = rand(1,N_layers(i)).*(INFO.TYPES.max_thick(types_layers{i}) - INFO.TYPES.min_thick(types_layers{i})) + INFO.TYPES.min_thick(types_layers{i});    % Randomly draw thicknesses of layers  
        elseif r(i) > INFO.SECTIONS.frequency(i)
            thick_sections(i) = 0;
            N_layers(i) = 0;
            types_layers{i} = [];
            thick_layers{i} = [];
        end
    end          
    
    
    % Make sure thicknesses fit within their intervals
    if N > 1
        for i = find(thick_sections ~= 0)
            thick_layers{i} = thick_layers{i}./(sum(thick_layers{i})/thick_sections(i));
        end
    end
    
    
    % Check if dimensions fit
    checksum = 0;
    for i = find(thick_sections ~= 0)
        layers_max_check = sum(thick_layers{i} >= 1.1*INFO.TYPES.max_thick(types_layers{i}));
        layers_min_check = sum(thick_layers{i} <= (1/1.1)*INFO.TYPES.min_thick(types_layers{i})); 
        checksum = checksum + layers_max_check + layers_min_check;
    end


    % If 2000 loops continue with a warning
    tries = tries + 1;
    if tries > 2000
        flag_vector(1) = 1;
        break
    end
end


% Combine into arrays
Ts_all = cat(2,thick_layers{:});
types_all = cat(2,types_layers{:});


% Convert to depths
Ds = cumsum(Ts_all);


% Fill into m vector
for i = numel(types_all):-1:1
    m(z < Ds(i)) = types_all(i);
    layer_index(z < Ds(i)) = layer_count;
    layer_count = layer_count + 1;
end

