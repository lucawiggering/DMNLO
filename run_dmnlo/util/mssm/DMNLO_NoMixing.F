* NoMixing.f
* this file is part of FormCalc
* last modified 15 Jun 04 th

	subroutine NoMixing(m, msq, U, M11, M22, M12, *)
	implicit none
	
	double precision m(*), msq(*)
	double complex U(2,*)
	double precision M11, M22
	double complex M12
	
	
	
	msq(1) = M11
	m(1) = sqrt(M11)
	msq(2) = M22
	m(2) = sqrt(M22)
	U(1,1) = 1
	U(2,2) = 1
	U(2,1) = 0
	U(1,2) = 0
	
	end
