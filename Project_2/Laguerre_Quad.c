/* Laguerre_Quad.f -- translated by f2c (version 20100827)
                      and then modified by hand (20181106).
*/

#include <stdio.h>
#include <math.h>

/* Main program */ int main(void)
{
    /* System generated locals */
    long int i__1;
    double d__1, d__2, d__3;

    /* Local variables */
    static long int i__, m, n;
    extern /* Subroutine */ int laguerreq_(long int *, double *, 
	    double *);
    static double wa[24], xa[24], ro, r3s, sum;

    for (n = 4; n <= 24; n += 4) {
	laguerreq_(&n, xa, wa);
	printf("\n\n\n");
	printf("Quadrature à %ld points :\n",n);
	i__1 = n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    printf("           %25.17E %25.17E \n",xa[i__ - 1],wa[i__ - 1]);
	}
/*     Vérification \int_0^\infty x^p.e^{-x}dx=p! */
	sum = 0.;
	r3s = 0.;
	i__1 = n;
	for (m = 1; m <= i__1; ++m) {
/*           sum=sum+wa(m)*dexp(-Xa(m)) */
/* Computing 2nd power */
	    d__1 = xa[m - 1];
	    sum += wa[m - 1] * (d__1 * d__1 - xa[m - 1] * 2 + 3);
	    ro = xa[m - 1] * 3 / 2;
/* Computing 2nd power */
	    d__2 = ro;
/* Computing 2nd power */
	    d__1 = 2 / sqrt(3.) / 81 * (27 - 18 * ro + 2 * (d__2 * d__2));
/* Computing 3rd power */
	    d__3 = ro;
	    r3s += wa[m - 1] * (d__1 * d__1) * (d__3 * (d__3 * d__3));
	}
	printf("Vérification :%25.17E, vs 3\n", sum);
	d__1 = r3s * 3 / 2;
	printf("<r>_3s :%25.17E, vs. 13.5\n", d__1);
    }
    return 0;
}

/* Subroutine */ int laguerreq_(long int *nx, double *x_q__, double *
	w_q__)
{
    /* Initialized data */

    static struct {
	double e_1[144];
	} equiv_15 = { .3225476896193923, .170279632305101, .1157221173580207,
		 .08764941047892784, .07053988969198876, .05901985218150798, 
		1.745761101158347, .9037017767993799, .6117574845151307, 
		.4626963289150808, .3721268180016115, .3112391461984837, 
		4.536620296921128, 2.251086629866131, 1.512610269776419, 
		1.141057774831227, .9165821024832736, .7660969055459366, 
		9.395070912301133, 4.266700170287659, 2.833751337743507, 
		2.129283645098381, 1.707306531028344, 1.425597590803613, 0., 
		7.045905402393466, 4.599227639418348, 3.437086633893207, 
		2.749199255309432, 2.29256205863219, 0., 10.75851601018099, 
		6.844525453115177, 5.078018614549768, 4.048925313850887, 
		3.370774264208998, 0., 15.740678641278, 9.621316842456866, 
		7.070338535048234, 5.615174970861617, 4.665083703467171, 0., 
		22.86313173688926, 13.00605499330635, 9.438314336391938, 
		7.459017453671063, 6.181535118736765, 0., 0., 
		17.11685518746226, 12.21422336886616, 9.594392869581097, 
		7.927539247172152, 0., 0., 22.15109037939701, 
		15.44152736878162, 12.03880254696432, 9.912098015077706, 0., 
		0., 28.487967250984, 19.18015685675314, 14.81429344263074, 
		12.14610271172977, 0., 0., 37.09912104446692, 
		23.51590569399191, 17.94889552051938, 14.64273228959667, 0., 
		0., 0., 28.57872974288214, 21.47878824028501, 
		17.41799264650898, 0., 0., 0., 34.58339870228663, 
		25.45170279318691, 20.49146008261642, 0., 0., 0., 
		41.94045264768833, 29.93255463170061, 23.88732984816973, 0., 
		0., 0., 51.70116033954332, 35.013434240479, 27.63593717433271,
		 0., 0., 0., 0., 40.83305705672857, 31.77604135237472, 0., 0.,
		 0., 0., 47.6199940473465, 36.35840580165162, 0., 0., 0., 0., 
		55.8107957500639, 41.45172048487077, 0., 0., 0., 0., 
		66.52441652561575, 47.15310644515632, 0., 0., 0., 0., 0., 
		53.60857454469506, 0., 0., 0., 0., 0., 61.05853144721876, 0., 
		0., 0., 0., 0., 69.96224003510503, 0., 0., 0., 0., 0., 
		81.49827923394888 };

    static struct {
	double e_1[144];
	} equiv_7 = { .6031541043416336, .3691885893416375, .2647313710554432,
		 .206151714957801, .1687468018511139, .1428119733347819, 
		.3574186924377997, .418786780814343, .377759275873138, 
		.3310578549508842, .2912543620060683, .2587741075174239, 
		.03888790851500538, .1757949866371718, .2440820113198776, 
		.2657957776442141, .2666861028670013, .2588067072728698, 
		5.392947055613274e-4, .03334349226121565, .09044922221168093, 
		.1362969342963775, .1660024532695068, .183322688977778, 0., 
		.002794536235225673, .0201023811546341, .04732892869412522, 
		.07482606466879237, .09816627262991889, 0., 
		9.076508773358213e-5, .002663973541865316, .01129990008033945,
		 .02496441730928322, .04073247815140865, 0., 
		8.485746716272531e-7, 2.032315926629994e-4, 
		.001849070943526311, .006202550844572237, .01322601940512016, 
		0., 1.04800117487151e-9, 8.365055856819799e-6, 
		2.042719153082785e-4, .001144962386476908, 
		.003369349058478304, 0., 0., 1.66849387654091e-7, 
		1.48445868739813e-5, 1.55741773027812e-4, 
		6.721625640935479e-4, 0., 0., 1.342391030515004e-9, 
		6.8283193308712e-7, 1.540144086522492e-5, 
		1.044612146592752e-4, 0., 0., 3.061601635035021e-12, 
		1.881024841079673e-8, 1.086486366517982e-6, 
		1.254472197799333e-5, 0., 0., 8.148077467426242e-16, 
		2.862350242973882e-10, 5.330120909556715e-8, 
		1.15131581273728e-6, 0., 0., 0., 2.127079033224103e-12, 
		1.757981179050582e-9, 7.96081295913363e-8, 0., 0., 0., 
		6.297967002517868e-15, 3.725502402512321e-11, 
		4.07285898755e-9, 0., 0., 0., 5.050473700035513e-18, 
		4.767529251578191e-13, 1.507008226292585e-10, 0., 0., 0., 
		4.161462370332855e-22, 3.372844243362438e-15, 
		3.917736515058451e-12, 0., 0., 0., 0., 1.155014339500399e-17, 
		6.894181052958086e-14, 0., 0., 0., 0., 1.539522140582344e-20, 
		7.819800382459448e-16, 0., 0., 0., 0., 5.286442725569158e-24, 
		5.350188813010038e-18, 0., 0., 0., 0., 1.656456612499023e-28, 
		2.010517464555503e-20, 0., 0., 0., 0., 0., 
		3.605765864552959e-23, 0., 0., 0., 0., 0., 
		2.451818845878403e-26, 0., 0., 0., 0., 0., 
		4.088301593680658e-30, 0., 0., 0., 0., 0., 
		5.575345788328357e-35 };


    /* System generated locals */
    long int i__1;

    /* Local variables */
#define a ((double *)&equiv_7)
    static long int j, n;
#define x ((double *)&equiv_15)
#define a1 ((double *)&equiv_7)
#define a2 ((double *)&equiv_7 + 19)
#define a3 ((double *)&equiv_7 + 38)
#define a4 ((double *)&equiv_7 + 57)
#define a5 ((double *)&equiv_7 + 76)
#define a6 ((double *)&equiv_7 + 95)
#define a7 ((double *)&equiv_7 + 114)
#define a8 ((double *)&equiv_7 + 133)
#define x1 ((double *)&equiv_15)
#define x2 ((double *)&equiv_15 + 19)
#define x3 ((double *)&equiv_15 + 38)
#define x4 ((double *)&equiv_15 + 57)
#define x5 ((double *)&equiv_15 + 76)
#define x6 ((double *)&equiv_15 + 95)
#define x7 ((double *)&equiv_15 + 114)
#define x8 ((double *)&equiv_15 + 133)


/*     Gauss-Laguerre quadrature for nX points(nX/4=1,6) */
/*     associated to the weight function W(x)=exp(-x) */
/*     Local variables */
    /* Parameter adjustments */
    --w_q__;
    --x_q__;

    /* Function Body */

    n = *nx / 4;
    *nx = n << 2;
    if (n < 1 || n > 6) {
	printf("N/4 out of range : %ld\n",n);		
	return -1;
    }
    i__1 = *nx;
    for (j = 1; j <= i__1; ++j) {
	x_q__[j] = x[n + j * 6 - 7];
	w_q__[j] = a[n + j * 6 - 7];
    }
    return 0;
} /* laguerreq_ */

#undef x8
#undef x7
#undef x6
#undef x5
#undef x4
#undef x3
#undef x2
#undef x1
#undef a8
#undef a7
#undef a6
#undef a5
#undef a4
#undef a3
#undef a2
#undef a1
#undef x
#undef a


