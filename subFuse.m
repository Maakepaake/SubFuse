function [A,B,C,D,data_to_use_final,model_order] = subFuse(DU,DY)
%SUBFUSE - LTI system identification with MOESP algorithm. 
%
% Syntax:  [output1,output2,output3,output4,output5,output6] = SubFuse(input1,input2)
%
% Inputs:
%    input1 - Matrix of input deviation from nominal at each timestep. Each
%    row is to represent an input signal; each column is to represent a
%    sample.
%    input2 - Matrix of output deviation from nominal at each timestep. Each
%    row is to represent a measurement signal; each column is to represent a
%    sample.
%
% Outputs:
%    output1 - Identified A matrix of LTI system state-space presentation
%    output2 - Identified B matrix of LTI system state-space presentation
%    output3 - Identified C matrix of LTI system state-space presentation
%    output4 - Identified D matrix of LTI system state-space presentation
%    output5 - Scalar amount of measurements used for identification
%    output6 - Scalar identified dimension of LTI system
%
% Author: Markus Neuvonen
% University of Oulu
% email: markus.neuvonen@oulu.fi
% February 2021; Last revision: 20-July-2021
%------------- BEGIN CODE --------------
arguments
    DU double {mustBeFinite,mustBeReal}
    DY double {mustBeFinite,mustBeReal}
end
assert(size(DU,2) == size(DY,2),...
    'Inputs and outputs must contain equal amount of samples!');

%% Define suitable model searching parameters:
%Read from data:
number_of_inputs    = size(DU,1); 
number_of_outputs   = size(DY,1);
number_of_max       = max(number_of_outputs,number_of_inputs);

model_min         = 2;  %The estimated model dimension min...
model_max         = 20; %and max.
observation_min   = model_max+1; %"Strictly greater than model dimensions"
observation_step  = 2;           %Freely chosen step size.
%Data Hankel block matrices are used according to Chapter 6.5 of Katayama's
%book. No "future data" is actually needed, but "number of observations
%used to model identification" should still be restricted a bit so we can
%see the performance also on "untrained" data. But the main constraint is
%that the final data matrices (U and Y) should have at least double the
%columns compared to rows.
observation_max     = floor(size(DY,2) / (2*number_of_max+1)); %To avoid indexing problems with finite data

%% Variable initializations:
%Initialize model saving variables:
tablerow            = 1;
compare_eSum        = 1e10;
data_to_use_final   = 0;

%% Computation
%Search for stable system matrices by changing data size and model
%dimensions. Save the obtained information about stability to table. Save
%the best estimate parameters for usage in next steps:
combination_total   = size(model_min:model_max,2)*size(observation_min:observation_step:observation_max,2);
stabilityData       = table('Size',[combination_total 5],'VariableTypes',{'double','double','double','double','double'},'VariableNames',{'System dimension','Observation amount','Unstable modes','Largest abs(eigenvalue)','min(svd(obsv(Ae,Ce)))'});
for model_order_try = model_min:model_max
    for number_of_observations_try = observation_min:observation_step:observation_max
        clc
        disp([ 'Evaluating model/data combination ' num2str(tablerow) '/' num2str(combination_total) ] )

        data_to_use = 2*number_of_max*number_of_observations_try;
        stabilityData(tablerow,1) = num2cell(model_order_try);
        stabilityData(tablerow,2) = num2cell(number_of_observations_try);
        
        U_Data = Create_Data_Hankel_Matrix( DU, number_of_observations_try, data_to_use );
        Y_Data = Create_Data_Hankel_Matrix( DY, number_of_observations_try, data_to_use );
        [A_try,C_try,L21_try,U2_try,L11_try] = Create_System_Matrix(U_Data,Y_Data,number_of_outputs,model_order_try);
        
        stabilityData(tablerow,3) = num2cell(size(find(abs(eig(A_try))>1),1));
        stabilityData(tablerow,4) = num2cell(max(abs(eig(A_try))));
        
        if stabilityData{tablerow,3} == 0
            disp('Stable system estimate found!')
            if min(abs(diag(L11_try))) > 1e-6
                [obsv_value,D_try,B_try] = testObservability(A_try,C_try,number_of_observations_try,U2_try,L21_try,L11_try,number_of_inputs,number_of_outputs);
                stabilityData(tablerow,5) = num2cell(obsv_value);
                %Select the model by best LS fit to data:
                updateRequest = false;
                modelCell = {A_try, B_try, C_try, D_try};
                simuData = identitySimulator1by1(modelCell,DU,DY);
                if sum(sum((DY-simuData{1}).^2)) < compare_eSum
                    disp(['The best fit! ', num2str(sum(sum((DY-simuData{1}).^2))) ]);pause(1)
                    updateRequest = true;
                    compare_eSum = sum(sum((DY-simuData{1}).^2));
                end
                if updateRequest
                    A = A_try;
                    B = B_try;
                    C = C_try;
                    D = D_try;
                    model_order = model_order_try;
                    data_to_use_final = data_to_use;
                end
            else
                obsv_value = 0;
                stabilityData(tablerow,5) = num2cell(obsv_value);
                disp('WARNING! Input data matrix is singular!');pause(1)
            end
        end
        tablerow = tablerow + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  FUNCTION DEFINITIONS:  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function U_data = Create_Data_Hankel_Matrix( Input_Sequence, Number_of_Rows, Number_of_Columns  )
%% Done by Istvan Selek
% Following the algorithm by Tohru Katayama's book.
input_dimension = size( Input_Sequence, 1 );

U_data = zeros( input_dimension*Number_of_Rows, Number_of_Columns );

for i=1:Number_of_Columns
   
    U_data(:,i) = reshape( Input_Sequence(:, i:(i+Number_of_Rows-1) ), [], 1 );
    
end

end

function [A,C,L21,U2,L11] = Create_System_Matrix(U_Data,Y_Data,number_of_outputs,model_order)
%% Done by Istvan Selek
% Following the algorithm by Tohru Katayama's book.
W = [ U_Data ; Y_Data ]; % Data matrix
[ ~,R ] = qr( (W')  );
L = (R');

n_UData = size( U_Data , 1 );

L11   = L( 1:n_UData       , 1:n_UData );
L21   = L( (n_UData+1):end , 1:n_UData );
L22   = L( (n_UData+1):end , (n_UData+1):end ); 

[Us,S,~] = svd(L22);

U1 = Us(: , 1:model_order );
U2 = Us(: , (model_order+1):end );
S1 = S(1:model_order, 1:model_order );

O = U1*sqrt(S1); 


C = O(1:number_of_outputs,:);
if min(svd(O( 1:model_order,: ))) > 1e-6
    A = ( O( 1:model_order,: ) )\(  O( (number_of_outputs+1):( number_of_outputs + model_order ) , : ) );
else
    A = 5;
end


end

function simuData = identitySimulator1by1(systemCell,DU,DY)
%% Comparison of identified LTI system model to the original IO-data

% Initialize:
Asim = systemCell{1};
Bsim = systemCell{2};
Csim = systemCell{3};
Dsim = systemCell{4};

simulation_horizon = size(DY,2);

X_hat      = zeros( size(Asim,1) , simulation_horizon + 1 );
Y_hat      = zeros( size(Csim,1) , simulation_horizon     );
Y_hat2     = zeros( size(Csim,1) , simulation_horizon     );

% Simulation:
for k=1:simulation_horizon
   
    X_hat(:,k+1) = Asim*X_hat(:,k) + Bsim*DU(:,k);
    Y_hat(:,k)   = Csim*X_hat(:,k); %Without identified static gain.
    Y_hat2(:,k)   = Csim*X_hat(:,k) + Dsim*DU(:,k); %With static gain.
    
end

simuData = {Y_hat,Y_hat2};
end

function [obsv_value,D,B] = testObservability(A,C,number_of_observations,U2,L21,L11,number_of_inputs,number_of_outputs)
%% By Istvan Selek
O = Calculate_Observability_Matrix( A, C, number_of_observations );
M = ((U2')*L21)/L11;

[ML , MR] = Calculate_BD_Multiplier_Matrix( U2, O, M, number_of_inputs, number_of_outputs );

X = pinv(ML)*MR;

D = X( 1:number_of_outputs ,: );
B = X( (number_of_outputs+1):end ,: );

if min(svd(obsv( A , C ))) > 1e-6
    obsv_value = min(svd(obsv( A , C )));
    disp('   --> System also observable!');pause(0.1)
else
    obsv_value = 0;
    disp('Unfortunately the extended system not observable...');pause(0.1)
end
end

function [ML , MR] = Calculate_BD_Multiplier_Matrix( U2, O, M,  Input_Dimension, Ouput_Dimension )
%% By Istvan Selek

U2T = U2';

n_of_observations = size(U2T,2)/Ouput_Dimension;

Ml_D_block = cell( n_of_observations ,1 );
Ml_B_block = cell( n_of_observations ,1 );
Mr_Block   = cell( n_of_observations ,1 );

for i=1:(n_of_observations-1)
    
    Ml_D_block{i} =   U2T( :, ( Ouput_Dimension*(i-1) + 1  ):( Ouput_Dimension*i )  );
    Ml_B_block{i} = ( U2T( :, ( Ouput_Dimension*i+1 ):end ) )*O( 1:(end - Ouput_Dimension*(i-1) ) , :  );
   
    Mr_Block{i}   = M( : ,  ( Input_Dimension*(i-1) + 1  ):( Input_Dimension*i )  );
    
end

Ml_D_block{n_of_observations} = U2T( :, ( Ouput_Dimension*(n_of_observations-1) + 1 ):( Ouput_Dimension*n_of_observations )  );
Ml_B_block{n_of_observations} = zeros( size( Ml_B_block{ n_of_observations-1 })  );

Mr_Block{ n_of_observations } = M( : ,  ( Input_Dimension*(n_of_observations-1) + 1  ):end  );

ML = [ cell2mat( Ml_D_block ) cell2mat(Ml_B_block)  ];
MR =  cell2mat( Mr_Block );

end

function O = Calculate_Observability_Matrix( A, C, Number_of_observations )
%% By Istvan Selek
O_block     = cell( Number_of_observations-1 ,1 );
O_block{1}  = C;

for i=2:(Number_of_observations-1)
    
    O_block{i} = O_block{i-1}*A;
    
end

O = cell2mat( O_block ) ; 
end


end