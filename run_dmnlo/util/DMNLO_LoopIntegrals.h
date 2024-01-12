	double precision Lpi, LEg, Ln4pi

	parameter (Lpi   = 3.1415926535897932384626433832795029D0)
	parameter (LEg   = 0.5772156649015328606065120900824024D0)
	parameter (Ln4pi = 2.5310242469692907929778915943677900D0)

	complex*16 xa0
	complex*16 xb0,xb1,xb00,xb11
	complex*16 xbp0,xbp1,xbp00,xbp11
	complex*16 xc0,xc1,xc2,xc00,xc11,xc12,xc22
	complex*16 xc001,xc002,xc111,xc112,xc122,xc222
	complex*16 xd0,xd1,xd2,xd3,xd00,xd11,xd12,xd13,xd22,xd23,xd33
	complex*16 xd001,xd002,xd003,xd111,xd112,xd113,xd122,xd123,xd133,xd222,xd223,xd233,xd333

	double precision UVdiv, IRdiv, Div, xicoll, xiIR, muSc, EpsPole

	common /Aredcoef/ xa0
	common /Bredcoef/ xb0,xb1,xb00,xb11
	common /BPredcoef/ xbp0,xbp1,xbp00,xbp11
	common /Credcoef/ xc0,xc1,xc2,xc00,xc11,xc12,xc22,xc001,xc002,xc111,xc112,xc122,xc222
	common /Dredcoef/ xd0,xd1,xd2,xd3,xd00,xd11,xd12,xd13,xd22,xd23,xd33
	common /Dredcoef/ xd001,xd002,xd003,xd111,xd112,xd113,xd122,xd123,xd133,xd222,xd223,xd233,xd333
	common /renScales/ UVdiv, IRdiv, Div, xicoll, xiIR, muSc, EpsPole
