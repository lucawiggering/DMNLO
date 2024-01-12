* model.h

#include "../util/DMNLO_Model.h"
	
	double precision Pcm
      integer n1,n2,n3,n4
      common /global/ Pcm,n1,n2,n3,n4
	
* Additional code Indizes and masses
	double precision m1, m2, m3, m4
	double precision m1red, m2red, m3red, m4red
	integer chimeigen, sfmeigen, sftype, sfgen, ftype, fgen, xtype
	integer pstype, psgen, psmeigen, pttype, ptgen, ptmeigen, putype, pugen, pumeigen
	double precision testvar, testvar1, testvar2, testvar3, testvar4
	double precision treeperc(7)

	common /add_para/ m1, m2, m3, m4
	common /add_para/ m1red, m2red, m3red, m4red
	common /add_para/ chimeigen, sfmeigen, sftype, sfgen, ftype, fgen, xtype
	common /add_para/ pstype, psgen, psmeigen, pttype, ptgen, ptmeigen, putype, pugen, pumeigen
	common /add_para/ testvar, testvar1, testvar2, testvar3, testvar4
      common /add_para/ treeperc


