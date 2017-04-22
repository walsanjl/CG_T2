#ifndef _MYGL_H_
#define _MYGL_H_

#include<math.h>
#include "definitions.h"
#include <iostream>
#include <stdio.h>

#define DIM 4

void PipeLine();    //prototipo da funcao PipeLine()

double angulo = 180;

double M_modelviewprojection[4][4] =  { { 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 0.0},
                                            { 0.0, 0.0, 0.0, 0.0}
                                            };

double M_viewport[4][4] =  { { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0}
                                };

//// declaração de tipos ////
typedef struct {
    double R, G, B, A;   //cor
} tCor;

typedef struct {
    double x, y;   //posicao
    tCor cor;
} tPixel;

typedef struct {
    double x, y, z, w;   //coordenadas
} tVertice;

typedef struct Ponto{
	float x;
	float y;
	float z;
	Ponto *proximo;
}ponto;

typedef struct PtrPonto{
	Ponto *ponto;
	PtrPonto *proximo;
}ptrPonto;

typedef struct Face{
	PtrPonto *ptrP;
	Face *proximo;
}face;

Ponto *pontos = NULL;
Face *faces = NULL;

void imprimeArray(int n, double *a) { // ok
    int i;
    for(i = 0; i < n; i++) {
        printf("%.1f  ", a[i]);
    }
    printf("---------------------------\n");
}

void imprimeMatriz(double m[][DIM]) { // ok
    int i, j;
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            printf("%.1f  ", m[i][j]);
        }
        printf("\n");
    }
    printf("-----------------------------\n");
}

//*****************************************************************************
// Defina aqui as suas funções gráficas
//*****************************************************************************
void PutPixel( double x, double y, tCor cor) {
    unsigned int offset = (4*x)+(4*y*IMAGE_WIDTH);
    FBptr[offset+0]=cor.R;
    FBptr[offset+1]=cor.G;
    FBptr[offset+2]=cor.B;
    FBptr[offset+3]=cor.A;
}

void DrawLine(double x0, double y0, double x1,double y1, tCor cor0, tCor cor1 ) {
    double dx = x1-x0;
    double dy = y1-y0;
    double m = (double)dy/dx;

    if((m>=0.0) and (m<=1.0) and (x1>x0)) {     // 1 configuracao
        //interpolacao linear em x
        double vR = (cor1.R-cor0.R)/dx;
        double vG = (cor1.G-cor0.G)/dx;
        double vB = (cor1.B-cor0.B)/dx;
        double vA = (cor1.A-cor0.A)/dx;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(x_aux < x1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Nordeste
                x_aux++;
                y_aux++;
            } else {        //a esccolha atual serah Leste
                x_aux++;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    } //fim do if
    else if((m>=0.0) and (m<=1.0) and (x1<x0)) {     // 2 configuracao
        dx = dx*(-1);
        dy = dy*(-1);

        //interpolacao linear em x
        double vR = (cor1.R-cor0.R)/dx;
        double vG = (cor1.G-cor0.G)/dx;
        double vB = (cor1.B-cor0.B)/dx;
        double vA = (cor1.A-cor0.A)/dx;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(x_aux > x1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Sudoeste
                x_aux--;
                y_aux--;
            } else {        //a esccolha atual serah Oeste
                x_aux--;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    } //fim do if
    else if((m>1.0) and (y1>y0)) {      // 3 configuracao
        //interpolacao linear em y
        double vR = (cor1.R-cor0.R)/dy;
        double vG = (cor1.G-cor0.G)/dy;
        double vB = (cor1.B-cor0.B)/dy;
        double vA = (cor1.A-cor0.A)/dy;

        //simulando o espelhamento da reta em torno da reta x=y
        double aux = dx;
        dx = dy;
        dy = aux;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(y_aux < y1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Nordeste
                x_aux++;
                y_aux++;
            } else {        //a esccolha atual serah Norte
                y_aux++;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    }
    else if((m>1.0) and (y1<y0)) {  // 4 configuracao
        dx = dx*(-1);
        dy = dy*(-1);

        //interpolacao linear em y
        double vR = (cor1.R-cor0.R)/dy;
        double vG = (cor1.G-cor0.G)/dy;
        double vB = (cor1.B-cor0.B)/dy;
        double vA = (cor1.A-cor0.A)/dy;

        //simulando o espelhamento da reta em torno da reta x=y
        double aux = dx;
        dx = dy;
        dy = aux;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(y_aux > y1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Sudoeste
                x_aux--;
                y_aux--;
            } else {        //a esccolha atual serah Sul
                y_aux--;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    }
    else if((m<0.0) and (m>=-1.0) and (x0<x1)) {    // 5 configuracao
        dy = dy*(-1);

        //interpolacao linear em x
        double vR = (cor1.R-cor0.R)/dx;
        double vG = (cor1.G-cor0.G)/dx;
        double vB = (cor1.B-cor0.B)/dx;
        double vA = (cor1.A-cor0.A)/dx;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(x_aux < x1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Sudeste
                x_aux++;
                y_aux--;
            } else {        //a esccolha atual serah Leste
                x_aux++;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    }
    else if((m<0.0) and (m>=-1.0) and (x0>x1)) {    // 6 configuracao
        dx = dx*(-1);

        //interpolacao linear em x
        double vR = (cor1.R-cor0.R)/dx;
        double vG = (cor1.G-cor0.G)/dx;
        double vB = (cor1.B-cor0.B)/dx;
        double vA = (cor1.A-cor0.A)/dx;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(x_aux > x1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Noroeste
                x_aux--;
                y_aux++;
            } else {        //a esccolha atual serah Oeste
                x_aux--;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    }
    else if((m<-1.0) and (y1<y0)) {     // 7 configuracao
        dy = dy*(-1);

        //interpolacao linear em y
        double vR = (cor1.R-cor0.R)/dy;
        double vG = (cor1.G-cor0.G)/dy;
        double vB = (cor1.B-cor0.B)/dy;
        double vA = (cor1.A-cor0.A)/dy;

        //simulando o espelhamento da reta em torno da reta x=y
        double aux = dx;
        dx = dy;
        dy = aux;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(y_aux > y1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Sudeste
                x_aux++;
                y_aux--;
            } else {        //a esccolha atual serah Sul
                y_aux--;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    }
    else if((m<-1.0) and (y1>y0)) {     // 8 configuracao
        dx = dx*(-1);

        //interpolacao linear em y
        double vR = (cor1.R-cor0.R)/dy;
        double vG = (cor1.G-cor0.G)/dy;
        double vB = (cor1.B-cor0.B)/dy;
        double vA = (cor1.A-cor0.A)/dy;

        //simulando o espelhamento da reta em torno da reta x=y
        double aux = dx;
        dx = dy;
        dy = aux;

        double alpha = dy;
        double beta = -dx;
        double d = 2*dy-dx;
        double x_aux = x0;
        double y_aux = y0;
        PutPixel(x_aux,y_aux,cor0);
        double d_old = d;
        double d_new=0;
        while(y_aux < y1) {
            if(d_old > 0) { //a esccolha anterior foi Nordeste
                d_new = d_old+alpha+beta;
            } else {     //a esccolha anterior foi Leste
                d_new = d_old+alpha;
            }
            if(d_new > 0) { //a esccolha atual serah Noroeste
                x_aux--;
                y_aux++;
            } else {        //a esccolha atual serah Norte
                y_aux++;
            }
            cor0.R+=vR;
            cor0.G+=vG;
            cor0.B+=vB;
            cor0.A+=vA;
            PutPixel(x_aux,y_aux,cor0);
            d_old = d_new;
        } // fim do while
    }
}

void DrawTriangle(double x0, double y0, double x1, double y1, double x2, double y2,
                   tCor cor0, tCor cor1, tCor cor2 ) {
    DrawLine(x0, y0, x1, y1, cor0, cor1 );
    DrawLine(x1, y1, x2, y2, cor1, cor2 );
    DrawLine(x2, y2, x0, y0, cor2, cor0 );
}

//// funções do trabalho 2 ////

void DesenhaMacaco() {
    tCor cor = {255.0,0.0,0.0,0.0};
    double xs[3], ys[3];

    Face *ptrFace;
    PtrPonto *ptrP;
    for(ptrFace = faces ; ptrFace != NULL ; ptrFace = ptrFace->proximo){
        int i = 0;
        for(ptrP = ptrFace->ptrP ; ptrP != NULL ; ptrP = ptrP->proximo){
            xs[i] = ptrP->ponto->x;
            ys[i] = ptrP->ponto->y;
            i++;
        }
        DrawTriangle(xs[0], ys[0], xs[1], ys[1], xs[2], ys[2],
                     cor, cor, cor);
    }
}

void limpaTela() {
    int i,j;
    for(i=0; i<511; i++) {
        for(j=0; j<511; j++) {
            PutPixel( i, j, {50.0, 50.0, 50.0, 0.0});
        }
    }
}

void PipeLine() {
    //*************************************************************************
	// Matriz model
	//*************************************************************************

    // matriz de rotação em Y, com um ângulo de th = pi/5 radianos
    //double th = M_PI / 5;
    double R2[4][4] = { { cos(angulo*rad), 0.0, sin(angulo*rad), 0.0},
                        { 0.0, 1.0, 0.0, 0.0},
                        { -sin(angulo*rad), 0, cos(angulo*rad), 0.0},
                        { 0.0, 0.0, 0.0, 1.0}
                    };
    angulo += 2;	// incrementando o angulo

    // matrix Model composta de duas rotações
    double M_model[4][4] =  {   { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0}
                            };
    
    // codigo para multiplicar duas matrizes de ordem DIM x DIM
    int i, j, k;
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            M_model[i][j] = R2[i][j];
        }
    }

    //*************************************************************************
	// Parametros da camera
	//*************************************************************************
	double camera_pos[] =       {0.0,0.0,4.0}; // Posicao da camera no universo.
	double camera_lookat[] =    {0.0,0.0,0.0}; // Ponto para onde a camera esta olhando.
	double camera_up[] =        {0.0,1.0,0.0}; // 'up' da camera no espaco do universo.

    //*************************************************************************
	// Calculo do sistema ortonormal gerado a partir dos parametros da camera
	//*************************************************************************
	double z_camera[] = {0.0, 0.0, 0.0};
	double aux1[] = {0.0, 0.0, 0.0};
	//int i;
	for(i = 0; i < 3; i++) {
	    aux1[i] = -(camera_lookat[i] - camera_pos[i] );
	}
	double aux2[] = {0.0, 0.0, 0.0};
	for(i = 0; i < 3; i++) {
	    aux2[i] = camera_lookat[i] - camera_pos[i];
	}
	double norma = sqrt( aux2[0]*aux2[0] + aux2[1]*aux2[1] + aux2[2]*aux2[2] );
	for(i = 0; i < 3; i++) {
	    z_camera[i] = aux1[i] / norma;
	}

    double x_camera[] = {0.0, 0.0, 0.0};
    double aux3[] = {0.0, 0.0, 0.0};
    aux3[0] = camera_up[1]*z_camera[2] - camera_up[2]*z_camera[1];
    aux3[1] = -camera_up[0]*z_camera[2] + camera_up[2]*z_camera[0];
    aux3[2] = camera_up[0]*z_camera[1] - camera_up[1]*z_camera[0];
    norma = sqrt( aux3[0]*aux3[0] + aux3[1]*aux3[1] + aux3[2]*aux3[2] );
    for(i = 0; i < 3; i++) {
	    x_camera[i] = aux3[i] / norma;
	}

    double y_camera[] = {0.0, 0.0, 0.0};
    y_camera[0] = z_camera[1]*x_camera[2] - z_camera[2]*x_camera[1];
    y_camera[1] = -z_camera[0]*x_camera[2] + z_camera[2]*x_camera[0];
    y_camera[2] = z_camera[0]*x_camera[1] - z_camera[1]*x_camera[0];


    //*************************************************************************
	// Construcao da matriz view
	//*************************************************************************
	double Bt[4][4] = { {x_camera[0], x_camera[1], x_camera[2], 0.0},
                        {y_camera[0], y_camera[1], y_camera[2], 0.0},
                        {z_camera[0], z_camera[1], z_camera[2], 0.0},
                        {0.0, 0.0, 0.0, 1.0}
                        };

    double T[4][4] = {  {1.0, 0.0, 0.0, -camera_pos[0]},
                        {0.0, 1.0, 0.0, -camera_pos[1]},
                        {0.0, 0.0, 1.0, -camera_pos[2]},
                        {0.0, 0.0, 0.0, 1.0}
                        };

    double M_view[4][4] =  {   { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0},
                                { 0.0, 0.0, 0.0, 0.0}
                            };
    // M_view = Bt * T
    //int j,k;
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            for(k = 0; k < DIM; k++) {
                M_view[i][j] = M_view[i][j] + (Bt[i][k] * T[k][j]);
            }
        }
    }

    //*************************************************************************
	// Construcao da matriz Modelview
	//*************************************************************************
	double M_modelview[4][4] =  {   { 0.0, 0.0, 0.0, 0.0},
                                    { 0.0, 0.0, 0.0, 0.0},
                                    { 0.0, 0.0, 0.0, 0.0},
                                    { 0.0, 0.0, 0.0, 0.0}
                                };

	// M_modelview = M_view * M_model
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            for(k = 0; k < DIM; k++) {
                M_modelview[i][j] = M_modelview[i][j] + (M_view[i][k] * M_model[k][j]);
            }
        }
    }

    //*************************************************************************
	// Construcao da matriz de Projecao
	//*************************************************************************
	double d = 2.0;
	double M_projection[4][4] =  {  { 1.0, 0.0, 0.0, 0.0},
                                    { 0.0, 1.0, 0.0, 0.0},
                                    { 0.0, 0.0, 1.0, d},
                                    { 0.0, 0.0, -1.0/d, 0.0}
                                };

    //*************************************************************************
	// Construcao da matriz Modelviewprojection
	//*************************************************************************
	
	//M_modelviewprojection = M_projection * M_modelview;
	for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            for(k = 0; k < DIM; k++) {
                M_modelviewprojection[i][j] = M_modelviewprojection[i][j] + (M_projection[i][k] * M_modelview[k][j]);
            }
        }
    }

    //*************************************************************************
	// Conversão de coordenadas do espaço canônico para o espaço de tela.
	//*************************************************************************
    double S1[4][4] =  {    { 1.0, 0.0, 0.0, 0.0},
                            { 0.0, -1.0, 0.0, 0.0},
                            { 0.0, 0.0, 1.0, 0.0},
                            { 0.0, 0.0, 0.0, 1.0}
                        };

    double T2[4][4] =  {    { 1.0, 0.0, 0.0, 1.0},
                            { 0.0, 1.0, 0.0, 1.0},
                            { 0.0, 0.0, 1.0, 0.0},
                            { 0.0, 0.0, 0.0, 1.0}
                        };

    double w = 512.0;
    double h = 512.0;

    double S2[4][4] =  {    { (w-1)/2, 0.0, 0.0, 0.0},
                            { 0.0, (h-1)/2, 0.0, 0.0},
                            { 0.0, 0.0, 1.0, 0.0},
                            { 0.0, 0.0, 0.0, 1.0}
                        };

    double aux[4][4] =  {   { 0.0, 0.0, 0.0, 0.0},
                            { 0.0, 0.0, 0.0, 0.0},
                            { 0.0, 0.0, 0.0, 0.0},
                            { 0.0, 0.0, 0.0, 0.0}
                            };
    // M_viewport = S2 * T2 * S1
        // aux = S2 * T2
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            for(k = 0; k < DIM; k++) {
                aux[i][j] = aux[i][j] + (S2[i][k] * T2[k][j]);
            }
        }
    }
        // M_viewport = aux * S1
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            for(k = 0; k < DIM; k++) {
                M_viewport[i][j] = M_viewport[i][j] + (aux[i][k] * S1[k][j]);
            }
        }
    }
}

Ponto* buscarPonto(Ponto **p, int f){
	Ponto *aux;
	int cout = 1;
	for(aux = (*p);!(aux == NULL); aux = aux->proximo){
		if(f == cout){
			return aux;
		}
		cout++;
	}
}

void inserirFace(Face **ptrFace, PtrPonto *ptrPonto){
	Face *novo = NULL;
    novo = (Face*) malloc(sizeof(Face));
    if(novo == NULL){
    	printf("ERRO. Ponteiro não alocado\n");
        return;
    }
    novo->ptrP = ptrPonto;
    novo->proximo = NULL;
    if((*ptrFace) == NULL){
    	(*ptrFace) = novo;
    }else{
    	Face *aux;
    	for(aux = (*ptrFace);!(aux->proximo == NULL); aux = aux->proximo);
        aux->proximo = novo;
    }
}

void inserirPonto(Ponto **p, double x, double y, double z){
	Ponto *novo = NULL;
	double vetor[4] = {0,0,0,1};
    novo = (Ponto*) malloc(sizeof(Ponto));
    if(novo == NULL){
    	printf("ERRO. Ponteiro não alocado\n");
        return;
    }

    vetor[0] = x;
    vetor[1] = y;
    vetor[2] = z;

    // levando pra o espaco de recorte
    double v_clip[4] = {0,0,0,0};
    // codigo para multiplicar uma matriz de ordem DIM x DIM por um vertice de ordem DIM x 1
    int i, j;
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            v_clip[i] = v_clip[i] + (M_modelviewprojection[i][j] * vetor[j]);
        }
    }

    // levando pra o espaco canonico
    double v_canonic[4] = {0,0,0,0};
    for(i = 0; i < DIM; i++) {
	    v_canonic[i] = v_clip[i] / v_clip[3];
	}

    // levando pra o espaco de tela
    double v_screen[4] = {0,0,0,0};
    for(i = 0; i < DIM; i++) {
        for(j = 0; j < DIM; j++) {
            v_screen[i] = v_screen[i] + (M_viewport[i][j] * v_canonic[j]);
        }
    }

    // arredondando os valores no espaco de tela
    for(i = 0; i < DIM; i++) {
        v_screen[i] = round(v_screen[i]);
    }

    novo->x = v_screen[0];
    novo->y = v_screen[1];
    novo->z = v_screen[2];
    novo->proximo = NULL;

    if((*p) == NULL){
    	(*p) = novo;
    }else{
    	Ponto *aux;
    	for(aux = (*p);!(aux->proximo == NULL); aux = aux->proximo);
        aux->proximo = novo;
    }
}

void inserirPtrPontos(PtrPonto **ptrPonto, Ponto *p){
	PtrPonto *novo = NULL;
    novo = (PtrPonto*) malloc(sizeof(PtrPonto));
    if(novo == NULL){
    	printf("ERRO. Ponteiro não alocado\n");
        return;
    }

    novo->ponto = p;
    novo->proximo = NULL;

    if((*ptrPonto) == NULL){
    	(*ptrPonto) = novo;
    }else{
    	PtrPonto *aux;
    	for(aux = (*ptrPonto);!(aux->proximo == NULL); aux = aux->proximo);
        aux->proximo = novo;
    }
}

void lerArquivo(Ponto **pontos, Face **faces){
	FILE *arq = NULL;
    char c, aux[5];
    float x,y,z;
    char temp[10];
    int id = 0;
    PtrPonto *ptrPonto = NULL;

    arq = fopen("monkey_head2.obj","r");
    if(arq == NULL){
    	printf("ERRO. Falha ao abrir aquivo.\n");
    }

    while((c = fgetc(arq))!=EOF){
    	//printf("%c",c);
    	if(c == 'v'){
			c = fgetc(arq);
			if(c == ' '){
				fscanf(arq,"%f %f %f", &x, &y, &z);
				inserirPonto(pontos,x,y,z);
			}
    	}
    	if(c == 'f'){
    		c = fgetc(arq);
    		if(c == ' '){
    			while((c = fgetc(arq))!= EOF && c != '\n'){
    				if((c >= '0') && (c <= '9')){
    					temp[id] = c;
    					id++;
    				}
    				else{
    					temp[id] = '\0';
                        inserirPtrPontos(&ptrPonto,buscarPonto(pontos,atoi(temp)));
                        id = 0;
    				}
    			}
    			//-----garante que a ultima string seja inserida na estrutura-----
    			temp[id] = '\0';
                inserirPtrPontos(&ptrPonto,buscarPonto(pontos,atoi(temp)));
                id = 0;
                inserirFace(faces, ptrPonto);
                ptrPonto = NULL;
    		}
    	}
    }
}

#endif // _MYGL_H_
