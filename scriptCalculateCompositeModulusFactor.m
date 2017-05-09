M=importdata('./finiteEltComsolData/Energy_vs_Modulus_CompositeTissue.txt');
E=M.data(:,1)/1000; %Matrix modulus in kPa

%next line plots raw data of elastic energy vs matrix modulus
%figure(20);plot(E,M.data(:,2),'b-',E,M.data(:,3),'k-',E,M.data(:,4),'r-',E,M.data(:,5),'bo',E,M.data(:,6),'ko',E,M.data(:,7),'ro');
%next line plots the energy divided by strain^2 vs matrix modulus: expect collapse in limit of linear elasticity!
%figure(21);plot(E,M.data(:,2)*40^2,'b-',E,M.data(:,3)*20^2,'k-',E,M.data(:,4)*(40/3)^2,'r-',E,M.data(:,5)*40^2,'bo',E,M.data(:,6)*20^2,'ko',E,M.data(:,7)*(40/3)^2,'ro');
%legend('1um displacement, composite','2um displacement, composite','3um displacement, composite','1um displacement, homogeneous','2um displacement, homogeneous','3um displacement, homogeneous','Location','NorthWest');
%title('Normalized strain energy vs modulus');
%figure(20);title('Strain energy vs modulus');legend('1um displacement, composite','2um displacement, composite','3um displacement, composite','1um displacement, homogeneous','2um displacement, homogeneous','3um displacement, homogeneous','Location','NorthWest');

%next line plots the factor relating matrix modulus to composite modulus
%figure(22);plot(E,M.data(:,2)./M.data(:,5),'b-',E,M.data(:,3)./M.data(:,6),'k-',E,M.data(:,4)./M.data(:,7),'r-');title('Effective Modulus Normalization');
%xlabel('E_m (kPa)');ylabel('Norm Factor');
%figure(23);loglog(E,M.data(:,2)./M.data(:,5),'b-');figure(24);semilogx(E,M.data(:,2)./M.data(:,5),'b-');
fo_power2=fit(E,M.data(:,2)./M.data(:,5),'power2');
[Em,str]=modulusRenormPow2Fit(fo_power2);
%figure(25);plot(Em,str,'b-');

normFactor=mean([M.data(:,2)./M.data(:,5),M.data(:,3)./M.data(:,6),M.data(:,4)./M.data(:,7)],2);
