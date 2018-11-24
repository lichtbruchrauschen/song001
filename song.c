#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define HP 1.5707963267948966
#define OP 3.141592653589793
#define DP 6.283185307179586
#define SR 96000.0
#define HR 48000.0
#define HN 0.9438743125896005
#define LC 32.703195672147295

typedef struct {
	int len;
	double buf[256];
} IREV, *PREV;

typedef struct {
	double phase[8];
	double inc[8];
	double vol;
	double hi;
	double buf[128];
} IHOSC, *PHOSC;

typedef struct {
	int i;
	double phase;
	double f;
	double fdec;
	double fmul;
	double vol;
	double vdec;
	int len;
	char on;
} IBD, *PBD;

typedef struct {
	double phase[5];
	double inc[5];
	double vol;
	double lo;
	double xlo;
	double qphase;
	double qinc;
} ISAW, *PSAW;

typedef struct {
	double lo;
	double hi;
	double max;
	double mul;
	double qlo;
	double qhi;
} IPASS, *PPASS;

typedef struct {
	IPASS low;
	IPASS mid;
	IPASS hih;
} IFX, *PFX;

void init_bpass(PPASS pass, double von, double bis) {
	pass->lo = 0.0;
	pass->hi = 0.0;
	pass->max = 0.5;
	pass->mul = 1.0;
	pass->qhi = von / HR;
	pass->qlo = bis / HR;
}

void init_fx(PFX fx) {
	init_bpass(&fx->low, 30.0, 300.0);
	init_bpass(&fx->mid, 400.0, 4000.0);
	init_bpass(&fx->hih, 5000.0, 20000.0);
}

double bpass(PPASS pass, double samin) {
	double out;

	pass->lo -= (pass->lo - samin) * pass->qlo;
	out = pass->lo;
	pass->hi += (out - pass->hi) * pass->qhi;
	out -= pass->hi;

	double v;
	double m;

	v = out;
	if (v < 0.0) v = -v;
	if (v > pass->max) pass->max =v;
	else {
		if (pass->max > 0.3) pass->max *= 0.9999875;
	}
	m = 0.7 / pass->max;
	pass->mul *= 15.0;
	pass->mul += m;
	pass->mul /= 16.0;

	out *= pass->mul;
	return out;
}

double run_fx(PFX fx, double samin) {
	double out;

	out = bpass(&fx->low, samin);
	out += bpass(&fx->mid, samin);
	out += bpass(&fx->hih, samin);

	if (out < -1.0) out = -1.0;
	if (out > 1.0) out = 1.0;

	return out;
}

void start_bdrum(PBD bd) {
	bd->i = 0;
	bd->phase = 0.0;
	bd->f = 500.0;
	bd->vol = 0.8;
	bd->len = (int)(0.35 * SR);
	bd->vdec = 0.8 / (double)bd->len;
	double wlen;
	wlen = SR / 1000.0;
	bd->fmul = 0.2;
	bd->fdec = (bd->f * bd->fmul) / wlen;
	bd->on = 1;
}

double bdrum(PBD bd, double samin) {
	if (bd->on == 0)
		return samin;

	double out;
	out = samin;
	double sil;
	sil = 1.0 - bd->vol;
	out *= sil;
	out += (sin(bd->phase) * bd->vol);

	if (bd->i > bd->len / 2) bd->vol -= bd->vdec * 2.0;
	if (bd->vol < 0.0) bd->vol = 0.0;

	double wlen;
	double inc;
	double intense;
	intense = 16.0 / (double)(bd->i + 1);
	intense *= 3.0;
	intense += 1.0;
	intense /= 4.0;
	inc = (DP / SR) * bd->f * (1.0 + (intense * cos(bd->phase)));
	if (bd->f > 50.0)
		bd->f -= bd->fdec;
	else
		bd->f *= 0.999975;
	bd->phase += inc;
	if (bd->phase > DP) {
		bd->phase -= DP;
		if (bd->fmul > 0.1) bd->fmul -= 0.001;
		wlen = SR / bd->f;
		bd->fdec = (bd->f * bd->fmul) / wlen;
	}

	bd->i++;
	if (bd->i >= bd->len) bd->on = 0;

	return out;
}

double bassor(PSAW insaw, PSAW dstsaw) {
	double out;
	double v;
	int i;

	out = 0.0;
	for (i = 0; i < 5; i++) {
		out += cos(insaw->phase[i]);
		insaw->phase[i] += insaw->inc[i];
		if (insaw->phase[i] > OP)
			insaw->phase[i] -= OP;
		insaw->inc[i] *= 31.0;
		insaw->inc[i] += dstsaw->inc[i];
		insaw->inc[i] /= 32.0;
	}
	v = out; if (v < 0.0) v = -v;
	out *= v * insaw->vol;

	insaw->vol *= 31.0;
	insaw->vol += dstsaw->vol;
	insaw->vol /= 32.0;

	insaw->qphase += insaw->qinc;
	if (insaw->qphase > DP) insaw->qphase -= DP;

	insaw->lo -= (insaw->lo - out) * (0.07 + (0.065 * sin(insaw->qphase)));
	out = insaw->lo;
	insaw->xlo -= (insaw->xlo - out) * (0.07 + (0.065 * sin(insaw->qphase)));
	out = insaw->xlo;

	return out;
}

double hihat(PHOSC osc) {
	double out;
	int i;
	double f;
	int n;
	int tune;

	out = 0.0;
	for (i = 0; i < 4; i++) {
		out += osc->vol * cos(osc->phase[i]);
		osc->phase[i] += osc->inc[i];
		if (osc->phase[i] > OP || osc->inc[i] < 0.000001) {
			osc->phase[i] -= OP;
			if (rand()%10 == 0) {
			tune = (rand()%4) * 3;
			f = 300.0; for (n = 0; n < tune; n++) f /= HN;
			osc->inc[i] = f;
			}
		}
	}
	osc->vol *= 0.9995;

	f = osc->buf[0] * 0.9;
	memcpy(osc->buf, &osc->buf[1], sizeof(double) * 127);
	osc->buf[127] = f + (out * 0.4);
	out += osc->buf[0];

	osc->hi += (out - osc->hi) * 0.15;
	out -= osc->hi;

	return out;
}

void set_hat(PHOSC osc) {
	osc->vol = 0.25;
}

void set_saw(PSAW saw, int tune, double volume) {
	double f;
	int i;
	double val;

	f = LC; for (i = 0; i < tune; i++) f /= HN;

	val = 1.0;
	for (i = 0; i < 5; i++) {
		saw->inc[i] = (DP / SR) * f * val;
		val += 0.505;
		saw->vol = volume / ((double)(i+i+1));
		if (i > 1) saw->vol *= 0.3;
	}
}

int main() {
	int i;
	int x;

	double ln32;
	double bpm;
	int songlen;
	int num;
	double lfo;

	num = 32 * 64;
	bpm = 125.0;
	ln32 = (7.5 / bpm) * SR;
	songlen = ln32 * (double)num;
	lfo = (bpm / 60.0);

	songlen += SR * 8.0;
	double *mono;
	mono = malloc(sizeof(double) * songlen);
	songlen -= SR * 7.0;

	int c;
	c = 0;
	double pos;
	pos = ln32;
	int next;
	next = (int)ln32;

	double e;
	double sam;

	double *echo;
	int elen;
	elen = (int)(ln32 * 6.0);
	elen++;
	echo = (double*) malloc(sizeof(double) * elen);
	memset(echo, 0, sizeof(double) * elen);
	elen--;

	IREV rev[8];
	for (i = 0; i < 8; i++) {
		rev[i].len = 96 + (rand()%192);
		memset(&rev[i].buf, 0, sizeof(double) * 256);
	}

	double buf[256];
	memset(buf, 0, sizeof(double) * 256);
	double lft;
	double rgt;
	int xsam;

	IFX lfx;
	IFX rfx;
	init_fx(&lfx);
	init_fx(&rfx);

	IBD bd;
	memset(&bd, 0, sizeof(IBD));

	FILE *outf;
	outf = fopen("tmp.raw", "wb");

	char txt[80];
	int tln;
	int oln;
	oln = 0;

	char ontick;
	ontick = 1;

	ISAW isbass;
	ISAW dstbass;
	set_saw(&isbass, 4, 0.1);
	isbass.lo = 0.0;
	isbass.xlo = 0.0;
	isbass.qphase = 0.0;
	isbass.qinc = (DP / SR) * lfo * 1.5;
	memcpy(&dstbass, &isbass, sizeof(ISAW));

	ISAW ismel;
	ISAW dstmel;
	set_saw(&ismel, 40, 0.0);
	ismel.lo = 0.0;
	ismel.xlo = 0.0;
	ismel.qphase = 0.0;
	ismel.qinc = (DP / SR) * lfo * 8.0;
	memcpy(&dstmel, &ismel, sizeof(ISAW));

	IHOSC hat;
	memset(&hat, 0, sizeof(IHOSC));

	int maintune;
	maintune = 0;

	for (i = 0; i < songlen; i++) {
		sam = ((double)(rand()%101) / 10000.0) - 0.005;
		if (rand()%3000==0)
			sam += ((double)(rand()%101) / 100.0) - 0.5;

		if (ontick == 1) {
		ontick = 0;

		if (c % 8 == 0 && c >= 256)
			start_bdrum(&bd);

		if (c >= 252 && c % 8 == 4)
			set_hat(&hat);

		if (c % 256 >= 240)
			set_saw(&dstmel, 40, 0.5);
		if (c%256==0)
			set_saw(&dstmel, 40, 0.0);

		if (c % 32 == 0)
			set_saw(&dstbass, 4, 0.2);
		if (c%32==4)
			set_saw(&dstbass, 4, 0.7);
		if (c%32==16)
			set_saw(&dstbass, 4, 0.1);
		if (c % 32 == 20)
			set_saw(&dstbass, 9, 0.5);
		if (c % 32 == 28)
			set_saw(&dstbass, 0, 0.6);
		}

		sam += hihat(&hat);
		sam += bassor(&ismel, &dstmel);
		sam += bassor(&isbass, &dstbass);

		e = echo[0] * 0.7;
		memcpy(echo, &echo[1], sizeof(double) * (elen-1));
		echo[elen] = (sam * 0.4) + e;
		sam += echo[0];

		sam = bdrum(&bd, sam);

		buf[250] = sam;
		memcpy(buf, &buf[1], sizeof(double) * 250);

		lft = 0.0; rgt = 0.0;
		for (x = 20; x < 30; x++) lft += buf[x];
		for (x = 200; x < 210; x++) rgt += buf[x];
		lft /= 10.0; lft = sam - lft; lft = run_fx(&lfx, lft);
		rgt /= 10.0; rgt = sam - rgt; rgt = run_fx(&rfx, rgt);

		xsam = (int)(lft * 5000000); fwrite(&xsam, 3, 1, outf);
		xsam = (int)(rgt * 5000000); fwrite(&xsam, 3, 1, outf);

		tln = i * 70;
		tln /= songlen;
		if (tln != oln) {
			oln = tln;
			memset(txt, 0, 80);
			for (x = 0; x < 70; x++) {
				if (x < tln) txt[x] = '#';
				else txt[x] = '_';
			}
			txt[x] = '\r';
			fputs(txt, stdout);
			fflush(stdout);
		}

		if (i > next) {
			pos += ln32;
			next = (int)pos;
			c++;
			ontick = 1;
		}
	}

	fclose(outf);
	puts("");

	return 0;
}
