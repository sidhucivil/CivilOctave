% Program to formulate global stiffness matrix
clc
clear

%% Initialisation [PRE-PROCESSING]

% Follwing are the main variables used. No explaination is given as
%	names are self explanatory.
% Number_of_storeys
% Mass
% Stiffness_story
%	(Summed shear stiffness of all columns in a story)


load input.mat

%Variable initialisation
Stiffness_matrix(Number_of_storeys, Number_of_storeys) = 0;
Time_period(Number_of_storeys, Number_of_storeys) = 0;


%% Display input

Number_of_storeys

Mass

Stiffness_story

%% [PROCESSING]

% Following loop adds local matrices of each element to form Global
% Matrix.

for storey_i = 1:Number_of_storeys

	Stiffness_matrix(storey_i, storey_i) = ...
		Stiffness_story(storey_i);

	if (storey_i < Number_of_storeys )
		Stiffness_matrix(storey_i, storey_i) = ...
			Stiffness_matrix(storey_i, storey_i) + ...
			Stiffness_story(storey_i + 1);
		Stiffness_matrix(storey_i, storey_i + 1) = ...
			- Stiffness_story(storey_i + 1);
                Stiffness_matrix(storey_i + 1, storey_i) = ...
			Stiffness_matrix(storey_i, storey_i + 1);
 	endif
    
end

[Eigen_vector, Omega_square] = eig(Stiffness_matrix, Mass);
Omega = sqrt(Omega_square);
for storey_i = 1:Number_of_storeys
  Time_period(storey_i, storey_i) = 2 * pi() ...
    / sqrt(Omega_square(storey_i, storey_i)); 
end

for storey_i = 1:Number_of_storeys
 Frequency(storey_i,1) = Omega(storey_i, storey_i);
end

for storey_i = 1:Number_of_storeys
 Time_periods(storey_i,1) = Time_period(storey_i, storey_i);
end

%% [POST-PROCESSING]

%% Echo input data

Number_of_storeys
Stiffness_story

%% Echo processed data

Stiffness_matrix

Eigen_vector
%Omega_square
Frequency
Time_periods

% End of file
