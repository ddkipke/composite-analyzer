%To run: put excel file "MaterialDat.xlsx" in the same folder as this script

clear
clc

numCer=8; %number of ceramic values stored
numMet=15; %"          metal "
numPoly=13; %number of polymer  values stored


filename='MaterialDat.xlsx';
ValueTable=xlsread(filename);


disp('You can input the properties of your material or choose the materials from a list')
InputType=input('Would you like to input your own properties? Y/N [Y] ','s');

if isempty(InputType) %default is yes
    InputType = 'Y';
end

%Routine to allow user to choose from stored values
validInput=1; %check to see if user put in acceptable value
while validInput %if user types invalid input, this loop will restart
    validInput=0;
	if InputType == 'N' 
		clc
		fprintf('First, choose a matrix type\nThe choices are:\n1 Ceramic\n2 Metal\n3 Polymer\n') %display choices
		MatrixClass=input('Type in the number of your choice: ','s'); %input choice
		validInput=1;
		while(validInput)
			validInput=0;
			if MatrixClass == '1'
				clc
				fprintf('The ceramic matrix choices are:\n1 Alumina\n2 Boron Carbide\n3 Cordierite\n4 Graphite\n5 Mullite\n6 Silicon Carbide\n7 Steatite\n8 Zirconia\n');
				MatrixType=input('Type in the number of your choice: ','s');
				validInput=1;
				while(validInput)
					validInput=0;
					if MatrixType == '1'
						MatrixName = 'Alumina';
					elseif MatrixType == '2' 
						MatrixName = 'Boron Carbide';
					elseif MatrixType == '3' 
						MatrixName = 'Cordierite';
					elseif MatrixType == '4' 
						MatrixName = 'Graphite';
					elseif MatrixType == '5' 
						MatrixName = 'Mullite';
					elseif MatrixType == '6'
						MatrixName = 'Silicon Carbide';
					elseif MatrixType == '7' 
						MatrixName = 'Steatite';
					elseif MatrixType == '8' 
						MatrixName = 'Zirconia';
					else
						clc
						validInput=1;
						disp('Invalid input, try again')
						fprintf('The ceramic matrix choices are:\n1 Alumina\n2 Boron Carbide\n3 Cordierite\n4 Graphite\n5 Mullite\n6 Silicon Carbide\n7 Steatite\n8 Zirconia\n');
						MatrixType=input('Type in the number of your choice: ','s');
					end
				end
				Em=ValueTable(str2double(MatrixType),1); %modulus of matrix
				Sm=ValueTable(str2double(MatrixType),2); %strength of matrix
				rhoM=ValueTable(str2double(MatrixType),3); %density of matrix
				x=['The strength of ',MatrixName,' is ',num2str(Sm),' MPa and the Modulus is ',num2str(Em),' GPa'];
				disp(x)
			elseif MatrixClass == '2'
				clc
				fprintf('The metal matrix choices are:\n1 Aluminum\n2 Beryllium\n3 Cobalt\n4 Copper\n5 Gold\n6 Iron\n7 Lead\n8 Magnesium\n9 Nickel\n10 Palladium\n11 Silver\n12 1080 Steel\n13 Tin\n14 Titanium\n15 Tungsten\n');
				MatrixType=input('Type in the number of your choice: ','s');
				validInput=1;
				while(validInput)
					validInput=0;
					if MatrixType == '1'
						MatrixName = 'Aluminum';
					elseif MatrixType == '2' 
						MatrixName = 'Beryllium';
					elseif MatrixType == '3' 
						MatrixName = 'Cobalt';
					elseif MatrixType == '4' 
						MatrixName = 'Copper';
					elseif MatrixType == '5' 
						MatrixName = 'Gold';
					elseif MatrixType == '6'
						MatrixName = 'Iron';
					elseif MatrixType == '7' 
						MatrixName = 'Lead';
					elseif MatrixType == '8' 
						MatrixName = 'Magnesium';
					elseif MatrixType == '9' 
						MatrixName = 'Nickel';
					elseif MatrixType == '10'
						MatrixName = 'Palladium';
					elseif MatrixType == '11' 
						MatrixName = 'Silver';
					elseif MatrixType == '12' 
						MatrixName = '1080 Steel';
					elseif MatrixType == '13' 
						MatrixName = 'Tin';
					elseif MatrixType == '14' 
						MatrixName = 'Titanium';
					elseif MatrixType == '15'
						MatrixName = 'Tungsten';
					else
						clc
						validInput=1;
						disp('Invalid input, try again')
						fprintf('The metal matrix choices are:\n1 Aluminum\n2 Beryllium\n3 Cobalt\n4 Copper\n5 Gold\n6 Iron\n7 Lead\n8 Magnesium\n9 Nickel\n10 Palladium\n11 Silver\n12 1080 Steel\n13 Tin\n14 Titanium\n15 Tungsten\n');
						MatrixType=input('Type in the number of your choice: ','s');
					end
				end
				Em=ValueTable(str2double(MatrixType)+numCer,1); %modulus of matrix
				Sm=ValueTable(str2double(MatrixType)+numCer,2); %strength of matrix
				rhoM=ValueTable(str2double(MatrixType)+numCer,3); %density of matrix
				x=['The strength of ',MatrixName,' is ',num2str(Sm),' MPa and the Modulus is ',num2str(Em),' GPa'];
				disp(x)
			elseif MatrixClass == '3'
				clc
				fprintf('The polymer matrix choices are:\n1 Acrylic\n2 Cellulose Acetate\n3 Cellulose Nitrate\n4 Epoxy\n5 Nylon\n6 Polyamide\n7 Polyester\n8 Polyethylene\n9 Polypropylene\n10 Polystyrene\n11 PTFE\n12 PVC\n13 Vinyl Ester\n');
				MatrixType=input('Type in the number of your choice: ','s');
				validInput=1;
				while(validInput)
				validInput=0;
					if MatrixType == '1'
						MatrixName = 'Acrylic';
					elseif MatrixType == '2' 
						MatrixName = 'Cellulose Acetate';
					elseif MatrixType == '3' 
						MatrixName = 'Cellulose Nitrate';
					elseif MatrixType == '4' 
						MatrixName = 'Epoxy';
					elseif MatrixType == '5' 
						MatrixName = 'Nylon';
					elseif MatrixType == '6'
						MatrixName = 'Polyamide';
					elseif MatrixType == '7' 
						MatrixName = 'Polyester';
					elseif MatrixType == '8' 
						MatrixName = 'Polyethylene';
					elseif MatrixType == '9' 
						MatrixName = 'Polypropylene';
					elseif MatrixType == '10'
						MatrixName = 'Polystyrene';
					elseif MatrixType == '11' 
						MatrixName = 'PTFE';
					elseif MatrixType == '12' 
						MatrixName = 'PVC';
					elseif MatrixType == '13' 
						MatrixName = 'Vinyl Ester';
					else
						clc
						validInput=1;
						disp('Invalid input, try again')
						fprintf('The polymer matrix choices are:\n1 Acrylic\n2 Cellulose Acetate\n3 Cellulose Nitrate\n4 Epoxy\n5 Nylon\n6 Polyamide\n7 Polyester\n8 Polyethylene\n9 Polypropylene\n10 Polystyrene\n11 PTFE\n12 PVC\n13 Vinyl Ester\n');
						MatrixType=input('Type in the number of your choice: ','s');
					end
				end
				Em=ValueTable(str2double(MatrixType)+numCer+numMet,1); %modulus of matrix
				Sm=ValueTable(str2double(MatrixType)+numCer+numMet,2); %strength of matrix
				rhoM=ValueTable(str2double(MatrixType)+numCer+numMet,3); %density of matrix
				x=['The strength of ',MatrixName,' is ',num2str(Sm),' MPa and the Modulus is ',num2str(Em),' GPa'];
				disp(x)

			else
				clc
				validInput=1;
				disp('Invalid input, try again');
				fprintf('First, choose a matrix type\nThe choices are:\n1 Ceramic\n2 Metal\n3 Polymer\n')
				MatrixClass=input('Type in the number of your choice: ','s');
			end
		end
		
		fprintf('\n\nNext, choose a fiber type\nThe choices are:\n1 Ceramic\n2 Metal\n3 Polymer\n')
		FiberClass=input('Type in the number of your choice: ','s');
		validInput=1;
		while(validInput)
			validInput=0;
			if FiberClass == '1'
				clc
				fprintf('The ceramic fiber choices are:\n1 Alumina\n2 Boron Carbide\n3 Cordierite\n4 Graphite\n5 Mullite\n6 Silicon Carbide\n7 Steatite\n8 Zirconia\n');
				FiberType=input('Type in the number of your choice: ','s');
				validInput=1;
				while(validInput)
					validInput=0;
					if FiberType == '1'
						FiberName = 'Alumina';
					elseif FiberType == '2' 
						FiberName = 'Boron Carbide';
					elseif FiberType == '3'
						FiberName = 'Cordierite';
					elseif FiberType == '4' 
						FiberName = 'Graphite';
					elseif FiberType == '5' 
						FiberName = 'Mullite';
					elseif FiberType == '6'
						FiberName = 'Silicon Carbide';
					elseif FiberType == '7' 
						FiberName = 'Steatite';
					elseif FiberType == '8' 
						FiberName = 'Zirconia';
					else
						clc
						validInput=1;
						disp('Invalid input, try again')
						fprintf('The ceramic fiber choices are:\n1 Alumina\n2 Boron Carbide\n3 Cordierite\n4 Graphite\n5 Mullite\n6 Silicon Carbide\n7 Steatite\n8 Zirconia\n');
						FiberType=input('Type in the number of your choice: ','s');
					end
				end
				Ef=ValueTable(str2double(FiberType),1); %modulus of matrix
				Sf=ValueTable(str2double(FiberType),2); %strength of matrix
				rhoF=ValueTable(str2double(FiberType),3); %density of matrix
				x=['The strength of ',FiberName,' is ',num2str(Sf),' MPa and the Modulus is ',num2str(Ef),' GPa'];
				disp(x)
			elseif FiberClass == '2'
				clc
				fprintf('The metal fiber choices are:\n1 Aluminum\n2 Beryllium\n3 Cobalt\n4 Copper\n5 Gold\n6 Iron\n7 Lead\n8 Magnesium\n9 Nickel\n10 Palladium\n11 Silver\n12 1080 Steel\n13 Tin\n14 Titanium\n15 Tungsten\n');
				FiberType=input('Type in the number of your choice: ','s');
				validInput=1;
				while(validInput)
					validInput=0;
					if FiberType == '1'
						FiberName = 'Aluminum';
					elseif FiberType == '2' 
						FiberName = 'Beryllium';
					elseif FiberType == '3' 
						FiberName = 'Cobalt';
					elseif FiberType == '4' 
						FiberName = 'Copper';
					elseif FiberType == '5' 
						FiberName = 'Gold';
					elseif FiberType == '6'
						FiberName = 'Iron';
					elseif FiberType == '7' 
						FiberName = 'Lead';
					elseif FiberType == '8' 
						FiberName = 'Magnesium';
					elseif FiberType == '9' 
						FiberName = 'Nickel';
					elseif FiberType == '10'
						FiberName = 'Palladium';
					elseif FiberType == '11' 
						FiberName = 'Silver';
					elseif FiberType == '12' 
						FiberName = '1080 Steel';
					elseif FiberType == '13' 
						FiberName = 'Tin';
					elseif FiberType == '14' 
						FiberName = 'Titanium';
					elseif FiberType == '15'
						FiberName = 'Tungsten';
					else
						clc
						validInput=1;
						disp('Invalid input, try again')
						fprintf('The metal fiber choices are:\n1 Aluminum\n2 Beryllium\n3 Cobalt\n4 Copper\n5 Gold\n6 Iron\n7 Lead\n8 Magnesium\n9 Nickel\n10 Palladium\n11 Silver\n12 1080 Steel\n13 Tin\n14 Titanium\n15 Tungsten\n');
						FiberType=input('Type in the number of your choice: ','s');
					end
				end
				Ef=ValueTable(str2double(FiberType)+numCer,1); %modulus of matrix
				Sf=ValueTable(str2double(FiberType)+numCer,2); %strength of matrix
				rhoF=ValueTable(str2double(FiberType)+numCer,3); %density of matrix
				x=['The strength of ',FiberName,' is ',num2str(Sf),' MPa and the Modulus is ',num2str(Ef),' GPa'];
				disp(x)
			elseif FiberClass == '3'
				clc
				fprintf('The polymer fiber choices are:\n1 Acrylic\n2 Cellulose Acetate\n3 Cellulose Nitrate\n4 Epoxy\n5 Nylon\n6 Polyamide\n7 Polyester\n8 Polyethylene\n9 Polypropylene\n10 Polystyrene\n11 PTFE\n12 PVC\n13 Vinyl Ester\n');
				FiberType=input('Type in the number of your choice: ','s');
				validInput=1;
				while(validInput)
					validInput=0;
					if FiberType == '1'
						FiberName = 'Acrylic';
					elseif FiberType == '2' 
						FiberName = 'Cellulose Acetate';
					elseif FiberType == '3' 
						FiberName = 'Cellulose Nitrate';
					elseif FiberType == '4' 
						FiberName = 'Epoxy';
					elseif FiberType == '5' 
						FiberName = 'Nylon';
					elseif FiberType == '6'
						FiberName = 'Polyamide';
					elseif FiberType == '7' 
						FiberName = 'Polyester';
					elseif FiberType == '8' 
						FiberName = 'Polyethylene';
					elseif FiberType == '9' 
						FiberName = 'Polypropylene';
					elseif FiberType == '10'
						FiberName = 'Polystyrene';
					elseif FiberType == '11' 
						FiberName = 'PTFE';
					elseif FiberType == '12' 
						FiberName = 'PVC';
					elseif FiberType == '13' 
						FiberName = 'Vinyl Ester';
					else
						clc
						validInput=1;
						disp('Invalid input, try again')
						fprintf('The polymer fiber choices are:\n1 Acrylic\n2 Cellulose Acetate\n3 Cellulose Nitrate\n4 Epoxy\n5 Nylon\n6 Polyamide\n7 Polyester\n8 Polyethylene\n9 Polypropylene\n10 Polystyrene\n11 PTFE\n12 PVC\n13 Vinyl Ester\n');
						FiberType=input('Type in the number of your choice: ','s');
					end
				end
				Ef=ValueTable(str2double(FiberType)+numCer+numMet,1); %modulus of matrix
				Sf=ValueTable(str2double(FiberType)+numCer+numMet,2); %strength of matrix
				rhoF=ValueTable(str2double(FiberType)+numCer+numMet,3); %density of matrix
				x=['The strength of ',FiberName,' is ',num2str(Sf),' MPa and the Modulus is ',num2str(Ef),' GPa'];
				disp(x)

			else
				clc
				validInput=1;
				disp('Invalid input, try again');
				fprintf('Next, choose a fiber type\nThe choices are:\n1 Ceramic\n2 Metal\n3 Polymer\n')
				FiberClass=input('Type in the number of your choice: ','s');
			end
    
		end
	  
	%Routine to allow user to input own values    
	elseif InputType == 'Y'
		Ef=input('Modulus of elasticity of fiber? (GPa)');
		Em=input('Modulus of elasticity of matrix? (GPa)');
	   
		Sf=input('Tensile strength of fiber? (MPa)');
		Sm=input('Tensile strength of matrix? (MPa)');
		
		rhoF=input('Density of the fiber? (g/cm^3)');
		rhoM=input('Density of the matrix? (g/cm^3)');
	  
	 
	else
		clc
		validInput=1;
		disp('Invalid input, try again');
		disp('You can input the properties of your material or choose the materials from a list')
		InputType=input('Would you like to input your own properties? Y/N [Y] ','s');
	end
end

Vf=input('\n\nVolume fraction of the fiber? (0-1)');
clc
Vm=1-Vf;


UB=0; %upper bound of E for long fibers
LB=0; %lower bound for E long fibers
Es=0; %E for short fibers


fprintf('To optimize longitudinal properties, fibers will be aligned.\n')
fprintf('Does your composite have continuous fibers? ');
fiblength=input('Y/N [Y] ','s');

if isempty(fiblength) %default is yes
    fiblength = 'Y';
end

if fiblength == 'Y'
    lengthFactor=1;%equations for continuous vs critical length fibers differ only by factor of 2 for fiber terms
elseif fiblength == 'N'
    lengthFactor=0.5;
    disp('You chose non-continuous fibers. To optimize properties, assume the fibers are of critical length')
end

clc
UB=lengthFactor*Vf*Ef+Vm*Em; %calculate upper bound of modulus
LB=Ef*Em/(lengthFactor*Vf*Em+(1-Vf)*Ef); %calculate lower bound
%Check if matrix or fiber fails first
FailStrainFib=Sf/Ef/1000; %divide by 10^3 for unitless strain answer
FailStrainMat=Sm/Em/1000;
    
% Make vectors to graph stress-strain later
StrainMat=0:0.0001:FailStrainMat;
StressMat=StrainMat*Em;
StrainFib=0:0.0001:FailStrainFib;
StressFib=StrainFib*Ef;

%calculate strength
if FailStrainFib > FailStrainMat %Check which part will fail first
    Sc=(lengthFactor*Ef*FailStrainMat*Vf+Em*FailStrainMat*Vm)*1000; %divide by 10^3 for answer in MPa rather than GPa
    StrainCom=0:0.001:FailStrainMat;
else
    Sc=(lengthFactor*Ef*FailStrainFib*Vf+Em*FailStrainFib*Vm)*1000;
    StrainCom=0:0.001:FailStrainFib;
end
StressCom=StrainCom*UB; %Hooke's law

%calculate specific strength
rho=rhoF*Vf+rhoM*Vm;
specstr=Sc/rho; %MPa/(g/cm^3) gives specific strenght in kN*m/kg
    
%Tell user name of final composite
if InputType == 'N'
     w=['You chose a ',MatrixName,' matrix with ',FiberName,' fibers.'];
     fprintf('\n\n')
     disp(w)
end
    
%Output results
x=['The longitudinal modulus of elasticity of your composite is ', num2str(UB), 'GPa'];
disp(x)
y=['The transverse modulus is ', num2str(LB), 'GPa'];
disp(y)
z=['The longitudinal tensile strength is ',num2str(Sc), 'MPa'];
disp(z)
zz=['The specific strength (strength to density ratio) is ',num2str(specstr),'kN*m/kg'];
disp(zz)
%graphing routine for Modulus vs volume fraction
Vf=0:0.01:1;
UB=Vf.*Ef+(1-Vf).*Em;
LB=Ef.*Em./(Vf.*Em+(1-Vf).*Ef);
plot(Vf,UB,'-r',Vf,LB,'--b')
xlabel('Volume fraction of fiber')
ylabel('E (GPa)')
legend('Longitudinal','Transverse')
    
%graphing routine for stress-strain
figure(2)
plot(StrainMat,StressMat,'-',StrainFib,StressFib,'--',StrainCom,StressCom,'-.')
xlabel('Strain (mm/mm)')
ylabel('Stress (MPa)')
if InputType == 'N'
    legend(MatrixName,FiberName,'Composite')
else
    legend('Matrix','Fiber','Composite')
end
    
  
