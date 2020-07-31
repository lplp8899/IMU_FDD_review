%-----------------------------------------------------------------------
%   
This is the code for the following paper:
%

%   Q. He, W. Zhang, P. Lu and J. Liu
%   
Performance comparison of representative model-based fault 
reconstruction algorithms for aircraft sensor fault detection 
and diagnosis
%  
%
%   if you find this code helpful, please cite the above paper.
% 

%   If you have more questions, please contact the correponding author:

%   Peng Lu, The Hong Kong Polytechnic University   

%   Email:  peng.lu@polyu.edu.hk
%   
Released in 2019

%----------------------------------------------------------------------


This package is the source code for the IMU fault diagnosis of aircraft.

Just run the 4 main function, it will display the results. 


Abbreviation:
	
1.ATSEKF: adaptive tow stage extended kalman filter
	
2.IOTSEKF: iterated optimal two stage extended kalman filter
	
3.NDO: nonlinear disturbance observer
	
4.SMO: sliding mode observer



%
Files describtion:


SCRIPT:
	
main_ATSEKF.m :main function of the ATSEKF method.
	
main_IOTSEKF.m :main function of the IOTSEKF method.
	
main_NDO.m :main function of the NDO method.
	
main_SMO.m :main function of the SMO method.
	

FUNCTION:
	E_lin.m
	F_Ja.m
	F_lin.m
	fx.m
	G_Ja.m
	G_lin.m
	H_Ja.m
	Matrix_get.m
	model_fx.m
	model_output.m
	sampling_g.m


DATA:
	new_sim_Ma_06_H_1000_notur.mat
