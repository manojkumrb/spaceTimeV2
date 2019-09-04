function  [allCorrect,pYttSave]=findInOutTol(devt,ytt,varYtt,tol,cutoff)
% points set to zero by cutoff





pYtt	=1-cdf('Normal',tol,ytt,sqrt(varYtt))+...
	cdf('Normal',-tol,ytt,sqrt(varYtt));

pYttSave	= pYtt;

pYtt(pYtt<cutoff)		=0;
pYtt(pYtt>=cutoff)		=1;

% Keep the order, else the logic changes
devTol					=zeros(size(devt));
devTol(devt<-tol)		=1;
devTol(devt>=tol)		=1;



allCorrect				=~xor(pYtt,devTol);


end