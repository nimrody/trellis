#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define Ninputs (1 << Nin)
#define Noutputs (1 << Nout)

#define INF (1.0e10)
#define sqr(x) ((x)*(x))

#define TRUNCLEN  20
#define MAXSTATES 64
#define MAXINPUTS  6


typedef unsigned char BYTE;				/* data - 8 bits per byte */
typedef unsigned long WORD;
typedef struct {float i,q;} SYMBOL;		/* 2D channel symbol */


/*
   This version fits a constant constellation: 8PSK and 16QAM are given
   as examples (with 4 and 8 state encoders respectively)
*/

/* trellis encoder

This encoder is designed to encoder Ungerboeck type trellis codes.
All such codes can be realized as a systematic encoder with feedback.

Therefore this implementation uses one linear feedback shift register
(example is for a rate 2/3 8-state encoder):

MSB x3[n] ------------------------------> z3[n] uncoded bit

	x2[n] --------+---------------------> z2[n]
			      | poly=02
LSB x1[n] --------|---+-----------------> z1[n]
		 	      |   |
			      V   V poly=01
			.-----------.
			|MSB|   |LSB|---+-----------> z0[n]
			`-----------'    |
			  ^              |
  poly=04	  `--------------'



or a rate 2/3 4-state encoder (say for 8PSK):

	x2[n] ------------------------------> z2[n]

LSB x1[n] ------------+-----------------> z1[n]
		 	          |
			          V poly=01
			    .-------.
			    |   |LSB|---+-----------> z0[n]
			    `-------'    |
			      ^          |
		 poly=02  `----------'




*/


unsigned poly8state2input[] = {
	04,
	01,
	02
};

unsigned poly8state3input[] = {
	04,
	01,
	02,

	00

};


unsigned poly8state4input[] = {
	04,
	01,
	02,

	00,
	00

};


unsigned poly4state2input[] = {
	02,
	01,
	00
};


int states_of_encoder(unsigned poly[])
{
	if (poly == poly4state2input)
		return 4;
	else if (poly == poly8state3input  || poly == poly8state4input
		||  poly == poly8state2input)
		return 8;
	else
		abort();
}


/* -- Constellation Mapping -- */

#define P8 2.73


double QAM8[8][2] = {
	{ P8,  0},          /*000 */
	{  1,  1},          /*001 */
	{  0, P8},          /*010 */
	{ -1,  1},          /*011 */
	{-P8,  0},          /*100 */
	{ -1, -1},          /*101 */
	{  0,-P8},          /*110 */
	{  1, -1}           /*111 */
};

double GRAY_QAM8[8][2] = {
	{-P8,  0},          /*000 */
	{ -1,  1},          /*001 */
	{  1,  1},          /*010 */
	{  0, P8},          /*011 */
	{ -1, -1},          /*100 */
	{  0,-P8},          /*101 */
	{ P8,  0},          /*110 */
	{  1, -1}           /*111 */
};

double QAM16[16][2] = {
	{1,1},
	{1,2},
	{0,0},
	{0,3},
	{1,3},
	{1,0},
	{0,2},
	{0,1},
	{3,3},
	{3,0},
	{2,2},
	{2,1},
	{3,1},
	{3,2},
	{2,0},
	{2,3}
};

double GRAY_QAM16[16][2] = {
	{3,0},          /*0 000 */
	{3,3},          /*0 001 */
	{3,1},          /*0 010 */
	{3,2},          /*0 011 */
	{0,0},          /*0 100 */
	{0,3},          /*0 101 */
	{0,1},          /*0 110 */
	{0,2},          /*0 111 */
	{2,0},          /*1 000 */
	{2,3},          /*1 001 */
	{2,1},          /*1 010 */
	{2,2},          /*1 011 */
	{1,0},          /*1 100 */
	{1,3},          /*1 101 */
	{1,1},          /*1 110 */
	{1,2},          /*1 111 */
};

double QAM32[32][2] = {
	{1,0},          /*00 000 */
	{0,2},          /*00 001 */
	{3,0},          /*00 010 */
	{5,1},          /*00 011 */
	{4,5},          /*00 100 */
	{5,3},          /*00 101 */
	{2,5},          /*00 110 */
	{0,4},          /*00 111 */
	{3,2},          /*01 000 */
	{2,0},          /*01 001 */
	{5,2},          /*01 010 */
	{3,3},          /*01 011 */
	{2,3},          /*01 100 */
	{3,5},          /*01 101 */
	{0,3},          /*01 110 */
	{2,2},          /*01 111 */


	{1,4},          /*10 000 */
	{2,4},          /*10 001 */
	{1,2},          /*10 010 */
	{1,1},          /*10 011 */
	{4,1},          /*10 100 */
	{3,1},          /*10 101 */
	{4,3},          /*10 110 */
	{4,4},          /*10 111 */
	{5,4},          /*11 000 */
	{4,2},          /*11 001 */
	{3,4},          /*11 010 */
	{1,5},          /*11 011 */
	{0,1},          /*11 100 */
	{1,3},          /*11 101 */
	{2,1},          /*11 110 */
	{4,0},          /*11 111 */
};



#define P2 0.707
double PSK8[8][2] = {
	{  1,  0},
	{ P2, P2},
	{  0,  1},
	{-P2, P2},
	{ -1,  0},
	{-P2,-P2},
	{  0, -1},
	{ P2,-P2}
};


double GRAY_PSK8[8][2] = {
	{  1,  0},      /* 000 */
	{ P2, P2},      /* 001 */
	{-P2, P2},      /* 010 */
	{  0,  1},		/* 011 */
	{ P2,-P2},		/* 100 */
	{  0, -1},		/* 101 */
	{ -1,  0},		/* 110 */
	{-P2,-P2}		/* 111 */
};

void init_constellations(void)
{
	int i;
	static alreadydone=0;

	if (alreadydone)
		return;


	/* Init TCM 16QAM & GRAY 16QAM: */
	for (i=0; i<16; i++) {
		QAM16[i][0]  = 2*QAM16[i][0] - 3.0;
		QAM16[i][1]  = 2*QAM16[i][1] - 3.0;

		GRAY_QAM16[i][0]  = 2*GRAY_QAM16[i][0] - 3.0;
		GRAY_QAM16[i][1]  = 2*GRAY_QAM16[i][1] - 3.0;
	}

	/* Init TCM 32QAM: */
	for (i=0; i<32; i++) {
		QAM32[i][0]  = 2*QAM32[i][0] - 5.0;
		QAM32[i][1]  = 2*QAM32[i][1] - 5.0;
	}

	/* 8PSK, 8QAM needs no modifications */


	alreadydone=1;
}

double constellation_power(double cons[][2], int Nbits)
{
	int i, n;
	double p=0;

	n = 1<<Nbits;
	for (i=0; i < n ; i++)
		p += sqr(cons[i][0]) + sqr(cons[i][1]);

	return p/n;
}

SYMBOL mapper(double cons[][2], int Nbits, BYTE bits[])
{
	SYMBOL s;
	BYTE j;
	int i;

	/* note that uncoded bits are HIGH order bits */
	j=0;
	for(i=0; i<Nbits; i++)
		j |= bits[i]<<i;

	s.i = cons[j][0];
	s.q = cons[j][1];

	return s;
}


/* -- Encoder --*/

/*
	encode() - TCM encoder

		cons[][2]  output symbols constellation
		poly[]     vector of polynomials of a systematic recursive encoder
		Nin        number of bits as input to encoder (Nout:=Nin+1)
		data       vector of binary data (0,1 per element - MSB first)
		Ndata      size of data vector

		symbols    returned vector of I&Q symbols generated by encoder+mapper
		returns    Nsymbols - length of generated symbols[] vector.
*/


int encode(double cons[][2], unsigned poly[], int Nin, BYTE *data, int Ndata,
		   SYMBOL *symbols)
{
	int Nout=Nin+1;
	unsigned enc_state=0;
	BYTE output[Nout];
	int j, nsymbols=0;

	while(Ndata > 0) {				/* for all bytes in input stream */

			/* -- processing per output/input symbol -- */

			/* first all input bits and rightmost register bit
			   are moved to the output: */

			output[0] = enc_state & 1;

			for(j=0; j<Nin; j++) {
				output[Nout-j - 1] = *data & 1;
				data++;
				Ndata--;
			}

			*symbols = mapper(cons, Nout, output);
			symbols++;
			nsymbols++;

			enc_state >>= 1;

			/* next we let output[] affect the register: */

			for (j=0; j<Nout; j++)
				enc_state ^=  output[j] ? poly[j] : 0;
	}
	return nsymbols;
}


/*
	modulate() - general modulator (no encoding. should be used w/gray mapper.

		cons[][2]  output symbols constellation
		Nin        number of bits per symbol
		data       vector of binary data (0,1 per element - MSB first)
		Ndata      size of data vector

		symbols    returned vector of I&Q symbols generated by encoder+mapper
		returns    Nsymbols - length of generated symbols[] vector.
*/


int modulate(double cons[][2], int Nin, BYTE *data, int Ndata,
			 SYMBOL *symbols)
{
	unsigned enc_state=0;
	BYTE output[Nin];
	int j, nsymbols=0;

	while(Ndata > 0) {				/* for all bytes in input stream */

			for(j=0; j<Nin; j++) {
				output[j] = *data & 1;
				data++;
				Ndata--;
				if (Ndata<1) {
					return nsymbols;
				}

			}

			*symbols = mapper(cons, Nin, output);
			symbols++;
			nsymbols++;

	}
	return nsymbols;
}




/* -- Decoding functions -- */

/* state transition table - used by decoder */
int src2dst[MAXSTATES][1 << MAXINPUTS];


/*
	generate_state_table - generates state table for decoder

		cons[][2]  output symbols constellation
		poly[]     vector of polynomials of a systematic recursive encoder
		Nin        number of bits as input to encoder (Nout:=Nin+1)
		printit	   1 if the table is to printed out (for manual verification)
*/


void generate_state_table(double cons[][2], unsigned poly[], int Nin, int printit)
{
	int Nout=Nin+1;
	int STATES = states_of_encoder(poly);
	unsigned enc_state, state;
	BYTE output[Nout];
	int i,input;

	if (printit) {
		printf("\n enc_state input --> nextstate  output\n");
		printf(  " --------- -----     ---------  ------\n\n");
	}

	for (enc_state=0; enc_state < STATES; enc_state++) {
		/* enc_state is the original state */


		/* for all possible inputs, we find the next state */
		for(input=0; input< Ninputs; input++) {
			for (i=0; i < Nin; i++)
				output[i+1] = (input >> i)&1;

			state = enc_state;
			output[0] = state & 1;

			state >>= 1;

			/* next we let output[] affect the register: */
			for (i=0; i<Nout; i++)
				state ^=  output[i] ? poly[i] : 0;

			/* and fill the table: */
			src2dst[enc_state][input] = state;
			if (printit) {
				SYMBOL s = mapper(cons, Nout, output);
				printf(" %10X %5X     %8X (%5.1f %5.1f) \n",
				  enc_state, input, state, s.i, s.q);
			}

		} /* for input */
	} /* for enc_state */
}



double compute_metric(double cons[][2], int Nout, WORD src_state, WORD input, SYMBOL *s)
{

	BYTE output[Nout];
	SYMBOL s1;
	double mi,mq;
	int j;

	/* in a systematic encoder the output is the input + one bit
	from the encoder state */

	output[0] = src_state & 1;
	for(j=1; j<Nout; j++)
		output[j] = (input >> (j-1))&1;

	s1 = mapper(cons, Nout, output);

	/* compute euclidean distance */

	mi = s->i - s1.i;
	mq = s->q - s1.q;

	return sqr(mi) + sqr(mq);
}



/*
	decode - TCM decoder

		cons[][2]  output symbols constellation
		poly[]     vector of polynomials of a systematic recursive encoder
		Nin        number of bits as input to encoder (Nout:=Nin+1)
		symbols    received symbols (I&Q)
		Nsymbols   length of symbols[] vector

		data       returned bits (one per vector element. MSB first)
		returns Ndata - length of data[] vector

*/

int decode(double cons[][2], unsigned poly[], int Nin, SYMBOL *symbols, int Nsymbols, BYTE *data)
{
	int STATES=states_of_encoder(poly);
	int Nout=Nin+1;
	double metric[MAXSTATES], metric1[MAXSTATES];
	WORD  src_state, dst_state, input;

	int i,j;
	unsigned symbolNo;
	int decodedbits=0;


	/* paths array stores survivors (one per current state). Nin
	bits are stored per decoded symbol. Constraint lengths of no
	more than 32 symbols are allowed when WORD is int32 */

	WORD paths[MAXSTATES][MAXINPUTS];
	WORD paths1[MAXSTATES][MAXINPUTS];

	/* clear path table */
	for(i=0; i<STATES; i++)
		for(j=0; j<Nin; j++) {
			paths[i][j] = 0L;
			paths1[i][j] = 0L;
		}


	/* initialize metric table so that the all-zero state is
	   chosen as the starting state (which it is in the encoder!)
	*/


	metric[0]=-INF;
	for(i=1; i<STATES; i++)
		metric[i] = 0;

	symbolNo = 0;

	while (symbolNo < Nsymbols) {
		/* first let's see if we have some bits to output (does not happen
		till we move past TRUNCLEN symbols) */

		if (symbolNo >= TRUNCLEN)  {
			double m1 = +INF; /* we look for the lowest metric */
			int    s1 = 0;

			/* find lowest metric path */
			for(i=0; i<STATES; i++)
				if (metric[i] < m1) {
					m1 = metric[i];
					s1 = i;
				}

			/* we output the oldest symbol(=input bits) of path s1 */
			for(i=Nin-1; i>=0; i--) {
				*data++ = (paths[s1][i] >> (TRUNCLEN-1) )&1;
				decodedbits++;
			}

		} /* if > TRUNCLEN*/

		/* for each destination state we iterate over all possible
		originating states and inputs [including parallel transitions!] .*/

		for (dst_state=0; dst_state < STATES; dst_state++) {

			double met1 = +INF; /* we look for the lowest metric */
			double m;
			int    in1 = 0;
			int    state1 = 0;

			for(src_state=0; src_state < STATES; src_state++) {
				for (input=0; input < Ninputs; input++) {
					if (src2dst[src_state][input] == dst_state) {

					 	m = compute_metric(cons, Nout, src_state, input, symbols);
						m += metric[src_state];
						if (m < met1) {
							met1 = m;
							in1  = input;
							state1 = src_state;
						}
					}
				} /* input */

			} /* src_state */

			/* fill in the new metric table (don't touch the current one
			as it is still used for further computations (dst_state might
			be a src_state for some other destination state) */

			metric1[dst_state] = met1;

			/* and ofcourse, dst_state takes the path 'state1' + chosen
			input bits 'in1', not forgetting to shift the registers to the
			left to make room for the new input symbol */

			for(i=0; i<Nin; i++) {
				paths1[dst_state][i] = paths[state1][i] << 1;
				paths1[dst_state][i] |= (in1 >> i)&1;
			}


		} /* dst_state */

		memcpy(paths, paths1, sizeof(paths));
		memcpy(metric, metric1, sizeof(metric));


		symbolNo++;	/* next received symbol */
		symbols++;
	}

	return decodedbits;
}



/*
	demodulate - simple minimum distance demodulation

		cons[][2]  output symbols constellation
		Nin        bits per symbol
		symbols    received symbols (I&Q)
		Nsymbols   length of symbols[] vector

		data       returned bits (one per vector element. MSB first)
		returns Ndata - length of data[] vector

*/

int demodulate(double cons[][2], int Nin, SYMBOL *symbols, int Nsymbols, BYTE *data)
{
	WORD  input;

	int i,j;
	unsigned symbolNo;
	int decodedbits=0;


	symbolNo = 0;

	while (symbolNo < Nsymbols) {
		double sq,si;
		double m, m1 = +INF;
		int    s1;

		/* we go over all constellation symbols to find the one
		closet to the received symbol */
		for(i=0; i<Ninputs; i++) {
			si = cons[i][0];
			sq = cons[i][1];
			m = sqr(si - symbols->i) + sqr(sq - symbols->q);
			if (m < m1) {
				s1 = i;
				m1 = m;
			}
		}

		for(i=0; i<Nin; i++) {
			*data++ = (s1 >> i) &1;
			decodedbits++;
		}

		symbolNo++;	/* next received symbol */
		symbols++;
	}

	return decodedbits;
}


/* SIMULATION --------------------------------------------------*/

/* channel */


#define	MAX_RANDOM	0x7fffffff

/* Generate gaussian random double with specified mean and std_dev */
double normal_rand(double mean, double std_dev)
{
	double fac,rsq,v1,v2;
	static double gset;
	static int iset;

	if(iset){
		/* Already got one */
		iset = 0;
		return mean + std_dev*gset;
	}
	/* Generate two evenly distributed numbers between -1 and +1
	 * that are inside the unit circle
	 */
	do {
		v1 = 2.0 * (double)random() / MAX_RANDOM - 1;
		v2 = 2.0 * (double)random() / MAX_RANDOM - 1;
		rsq = v1*v1 + v2*v2;
	} while(rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	gset = v1*fac;
	iset++;
	return mean + std_dev*v2*fac;
}


/* adds noise with power 2*stddev^2 (power divided equally
between inphase and quadrature components */

void addnoise(SYMBOL *s, int Nsymbols, double stddev)
{

	while (Nsymbols-- >0) {
		s->i = normal_rand( s->i, stddev);
		s->q = normal_rand( s->q, stddev);
		s++;
	}
}



#define TESTLEN  500000UL
#define TESTITER 5
#define TERMINATION 10



BYTE input[TESTLEN];
SYMBOL output[TESTLEN];
BYTE decoded[TESTLEN+10];

void print_bits(BYTE *d, int n)
{
	while(n-- > 0)
		printf("%1d", *d++);
	putchar('\n');

}



void print_symbols(SYMBOL *d, int n)
{
	while(n-- > 0) {
		printf("(%3.0f,%3.0f)", d->i, d->q);
		d++;
	}
	putchar('\n');

}


void testnoise()
{
	int i,j;
	double p,m;


	for (j=0; j<TESTITER;j++) {
		for (i=0; i<TESTLEN; i++)
			output[i].q = output[i].i = 0.0;

		addnoise(output, TESTLEN, 1.0);

		p=0;
		for (i=0; i<TESTLEN; i++) {
			p += sqr(output[i].i) + sqr(output[i].q);
			m += output[i].i + output[i].q;
		}

		p/= TESTLEN;
		m/= TESTLEN;

		printf("%d> Noise mean / power: %10f %10f\n", j, m,p);
	}
}


void test_coded()
{
	int s,d, i,j, err, totalerr=0;
	double S,Eb,N0, EbNo_db, EbNo;


	int Nin=4;
	int Nout=Nin+1;


	init_constellations();
	generate_state_table(QAM32, poly8state4input, Nin, 0);

	S = constellation_power(QAM32,Nout);
	printf("constellation power [db] = %f\n\n", 10*log10(S));

	EbNo_db = 8.5;
	EbNo = pow(10,EbNo_db/10.0);


	Eb = S / Nin;
	N0 = Eb / EbNo ;

	printf("Nin=%d ; Nout=%d\n",Nin, Nout);
	for (j=0; j<TESTITER;j++) {


		#if 1
 		for (i=0; i<TESTLEN-TERMINATION; i++)
			input[i]=rand()&1;

		for (; i<TESTLEN; i++)
			input[i]=0;
		#else
		BYTE input[]={1,1, 0,0, 1,1, 0,0,   1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		print_bits(input, 30);
		#endif

		s = encode(QAM32, poly8state4input, Nin, input, sizeof(input), output);
	//	print_symbols(output,s);
		addnoise(output,s, sqrt(N0/2) );
		d = decode(QAM32, poly8state4input, Nin, output, s, decoded);
	//	print_bits(decoded, d);

		err = 0;
		for (i=0; i<d; i++)
			if(input[i] != decoded[i])
				err++;

		printf("%5d> Eb/No[dB]=%5.2f  BER=%f (total symbols %d, bits %d, errs %d)\n",
		       j,
		       10.0*log10(Eb/N0),
		       err/(double)d,
		       s,d,err);

		totalerr += err;
	}

	printf("\n\n%5d iterations Eb/No[dB]=%5.2f  average BER=%f\n",
			   j,
		       10.0*log10(Eb/N0),
		       totalerr/((double)d*j));

}


int test_coded_QAM32(double EbNo_db)
{
	int s,d, i, err;
	double S,Eb,N0, EbNo;


	int Nin=4;
	int Nout=Nin+1;


	init_constellations();
	generate_state_table(QAM32, poly8state4input, Nin, 0);

	S = constellation_power(QAM32,Nout);

	EbNo = pow(10,EbNo_db/10.0);
	Eb = S / Nin;
	N0 = Eb / EbNo ;

	for (i=0; i<TESTLEN-TERMINATION; i++)
		input[i]=rand()&1;
	for (; i<TESTLEN; i++)
		input[i]=0;

	s = encode(QAM32, poly8state4input, Nin, input, sizeof(input), output);
	addnoise(output,s, sqrt(N0/2) );
	d = decode(QAM32, poly8state4input, Nin, output, s, decoded);

	err = 0;
	for (i=0; i<d; i++)
		if(input[i] != decoded[i])
			err++;

	return err;
}


int test_coded_QAM16(double EbNo_db)
{
	int s,d, i, err;
	double S,Eb,N0, EbNo;


	int Nin=3;
	int Nout=Nin+1;


	init_constellations();
	generate_state_table(QAM16, poly8state4input, Nin, 0);

	S = constellation_power(QAM16,Nout);

	EbNo = pow(10,EbNo_db/10.0);
	Eb = S / Nin;
	N0 = Eb / EbNo ;

	for (i=0; i<TESTLEN-TERMINATION; i++)
		input[i]=rand()&1;
	for (; i<TESTLEN; i++)
		input[i]=0;

	s = encode(QAM16, poly8state3input, Nin, input, sizeof(input), output);
	addnoise(output,s, sqrt(N0/2) );
	d = decode(QAM16, poly8state3input, Nin, output, s, decoded);

	err = 0;
	for (i=0; i<d; i++)
		if(input[i] != decoded[i])
			err++;

	return err;
}


int test_coded_PSK8(double EbNo_db)
{
	int s,d, i, err;
	double S,Eb,N0, EbNo;


	int Nin=2;
	int Nout=Nin+1;


	init_constellations();
	generate_state_table(PSK8, poly8state2input, Nin, 0);

	S = constellation_power(PSK8,Nout);

	EbNo = pow(10,EbNo_db/10.0);
	Eb = S / Nin;
	N0 = Eb / EbNo ;

	for (i=0; i<TESTLEN-TERMINATION; i++)
		input[i]=rand()&1;
	for (; i<TESTLEN; i++)
		input[i]=0;

	s = encode(PSK8, poly8state2input, Nin, input, sizeof(input), output);
	addnoise(output,s, sqrt(N0/2) );
	d = decode(PSK8, poly8state2input, Nin, output, s, decoded);

	err = 0;
	for (i=0; i<d; i++)
		if(input[i] != decoded[i])
			err++;

	return err;
}



int test_coded_QAM8(double EbNo_db)
{
	int s,d, i, err;
	double S,Eb,N0, EbNo;


	int Nin=2;
	int Nout=Nin+1;


	init_constellations();
	generate_state_table(QAM8, poly8state2input, Nin, 0);

	S = constellation_power(QAM8,Nout);

	EbNo = pow(10,EbNo_db/10.0);
	Eb = S / Nin;
	N0 = Eb / EbNo ;

	for (i=0; i<TESTLEN-TERMINATION; i++)
		input[i]=rand()&1;
	for (; i<TESTLEN; i++)
		input[i]=0;

	s = encode(QAM8, poly8state2input, Nin, input, sizeof(input), output);
	addnoise(output,s, sqrt(N0/2) );
	d = decode(QAM8, poly8state2input, Nin, output, s, decoded);

	err = 0;
	for (i=0; i<d; i++)
		if(input[i] != decoded[i])
			err++;

	return err;
}

void test_coded2(void)
{
	double ebno;
	int err, i;

	for(ebno=3.0; ebno<14; ebno+=0.5) {
		err=0;
		i=0;
		while(err<25) {
			err += test_coded_QAM8(ebno);
			i++;
		}
		printf("%e %e\n",  ebno, err/((double)i*(TESTLEN-TERMINATION)));
		fflush(stdout);
	}
}


void test_uncoded()
{
	int s,d, i,j, err, totalerr=0;
	double S,Eb,N0, EbNo_db, EbNo;


	int Nin=3; // 8QAM

	int Nout=Nin+1;


	init_constellations();

	S = constellation_power(GRAY_QAM8,Nin);
	printf("constellation power [db] = %f\n\n", 10*log10(S));

	EbNo_db = 13.0;
	EbNo = pow(10,EbNo_db/10.0);


	Eb = S / Nin;
	N0 = Eb / EbNo ;

	printf("Nin=%d ; Nout=%d\n",Nin, Nout);
	for (j=0; j<TESTITER;j++) {

 		for (i=0; i<TESTLEN-TERMINATION; i++)
			input[i]=rand()&1;

		s = modulate(GRAY_QAM8, Nin, input, sizeof(input), output);
		//print_symbols(output, s);
		addnoise(output,s, sqrt(N0/2) );
		d = demodulate(GRAY_QAM8, Nin, output, s, decoded);
		//print_bits(decoded, d);

		err = 0;
		for (i=0; i<d; i++)
			if(input[i] != decoded[i])
				err++;

		printf("%5d> Eb/No[dB]=%5.2f  BER=%f (total symbols %d, bits %d, errs %d)\n",
		       j,
		       10.0*log10(Eb/N0),
		       err/(double)d,
		       s,d,err);

		totalerr += err;
	}

	printf("\n\n%5d iterations Eb/No[dB]=%5.2f  average BER=%f\n",
			   j,
		       10.0*log10(Eb/N0),
		       totalerr/((double)d*j));
}



main()
{

#if 1
//	testnoise();
// 	test_coded();

//	test_uncoded();

	test_coded2();

#else
	BYTE input[]={1,1, 0,0, 1,1, 0,0,   1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	SYMBOL output[100];
	BYTE decoded[100];

	int i;
	int s, d;

	generate_state_table(0);
	s  = encode(output, input, sizeof(input));

	memset(decoded, 0, sizeof(decoded));
	d = decode(output, decoded,s);

	printf("Nin=%d ; Nout=%d\n",Nin, Nout);
	printf("s=%d\n\n", s);
	for (i=0; i<s; i++){
		printf("(%3.1f,%3.1f) ", output[i].i, output[i].q);
	}


	printf("\n\n %d bits decoded ----------------------------\n", d);


	for (i=0; i<d; i++){
		printf("%d ", decoded[i]);
	}
#endif


}

